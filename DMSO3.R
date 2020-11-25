install.packages("tidyverse")
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Hmisc)
install.packages("Hmisc")
install.packages("lmerTest")
library(lme4)

######   wrangle raw data   ######
rm(DMSOdrop)
DMSOdrop <- as_tibble(DMSOdrop)

#rename sac to treatment
DMSOdrop <- DMSOdrop %>% 
  rename(treatment = sac)

print(DMSOdrop, n= 50)

#making % change column
DMSOdrop.pct <- DMSOdrop %>%
  group_by(treatment, antpost, indvd) %>%
  mutate(
    area.pct.change = ((area - area[1]) / area[1]
    )*100) %>%
  ungroup() #why the ungroup?

print(DMSOdrop.pct, n=100)
str(DMSOdrop.pct)

#Change the treatment category treatment to baf
DMSOdrop.pct$treatment[DMSOdrop.pct$treatment == 'treatment'] <- 'DMSO'

print(DMSOdrop.pct, n= 48)

#dataframe with only the experimental (non-zero) %changes
exp <- subset(DMSOdrop.pct, exposure == "experimental", 
              select = c(treatment, antpost, area, indvd, area.pct.change))
print(exp, n=24)

## compute means and sds
means.sd <-
  DMSOdrop.pct %>% #still has baseline areas
  select(-indvd) %>% 
  group_by(exposure, treatment, antpost) %>% 
  ## now compute mean and sd:
  summarize(across(everything(), #how does the across everything work? is it that chr columns get created based on num columns that are summarized?
                   tibble::lst(mean = mean, sd = sd)))

print(means.sd, n= 50)

#take only esperimental means/ sd 
exp.mean.sd <- subset(means.sd, exposure == "experimental", 
                      select = c(treatment, antpost, area_mean, area_sd, area.pct.change_mean, area.pct.change_sd))

print(exp.mean.sd, n=24)


##### makes strip charts with mean and stdv (from percents)
ggplot(exp, aes(y = area.pct.change, x = treatment, colour= antpost)) +
  geom_jitter(size = 2, pch = 1, position = position_dodge(width = 0.7)) +
  labs(x = "exposure", y = "% change") + #labels axes
  theme_classic() +  #takes out background
  stat_summary(
    fun.data = mean_sdl, position = position_dodge(width = 0.5), geom = "errorbar", width = 0.1, fun.args = list(mult=1)) +
  stat_summary(
    fun = mean, geom = "point", position = position_dodge(width = 0.5),
    size = 3)


########### stats ###########

stats_data <-
  exp %>%
  mutate(treatment = as_factor(treatment)) %>%
  mutate(antpost = as_factor(antpost))

print(stats_data, n= 24)


# doesn't account for random effects
mod1 <-
  lm(area.pct.change ~ treatment * indvd, 
     data = stats_data)
summary(mod1)

mod1.aov <- aov(mod1)
TukeyHSD(mod1.aov)

####

stats_data$indvd <- as.factor(stats_data$indvd)

mcmod <-
  MCMCglmm::MCMCglmm(
    area.pct.change ~ treatment, random = ~indvd,
    data = stats_data, scale = FALSE,
    nitt = 1300000, thin = 1000, burnin = 300000, 
    verbose = FALSE
  )

summary(mcmod)

mean(mcmod$VCV[,1]/(mcmod$VCV[,1] + mcmod$VCV[,2]))

#with interaction
mcmod.inter <-
  MCMCglmm::MCMCglmm(
    area.pct.change ~ treatment * antpost, random = ~indvd,
    data = stats_data, scale = FALSE,
    nitt = 1300000, thin = 1000, burnin = 300000, 
    verbose = FALSE
  )

summary(mcmod.inter)

mean(mcmod.inter$VCV[,1]/(mcmod.inter$VCV[,1] + mcmod.inter$VCV[,2]))

print(stats_data)

z <- lmer(area.pct.change ~ treatment + (1|indvd), data = stats_data)
summary(z)
VarCorr(z)
