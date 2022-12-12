### Code for the analysis of longitudinal diving data (foraging dive frequency, max dive depth, descent rate) collected over 2016-2019 
### on known-age birds equipped with GDR tags
### By Amélie Lescroël
### June 1, 2022
### Revised Nov 24, 2022

### R version 4.1.2 (2021-11-01) -- "Bird Hippie"
### Copyright (C) 2021 The R Foundation for Statistical Computing
### Platform: x86_64-w64-mingw32/x64 (64-bit)

### Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggeffects)
library(performance)

### set WD
setwd(".../")

### Load dive data summarized by period of the annual cycle
by_p5 <- read.csv(".../by_period_dive_data.csv")

## Assign factor levels
by_p5$period <- factor(by_p5$period, levels = c("Incubation", "Chick rearing", "Moult", "Outgoing migration", "Migration max", "Returning migration"))
by_p5$age_class <- factor(by_p5$age_class, levels=c("Young", "Middle-age", "Old"))
by_p5$season <- factor(by_p5$season, levels=c("2016", "2017", "2018"))
by_p5$br_st <- factor(by_p5$br_st, levels=c("Non breeder", "Failed breeder", "Successful breeder"))
by_p5$sex <- factor(by_p5$sex, levels=c("F", "M"))
by_p5$colony <- factor(by_p5$colony, levels=c("ROYD", "CROZ"))

### Variation in foraging dive frequency
## Calculating year-centered foraging dive frequency (F_dives_PH = number of foraging dives per hour, 
## F_dives_PH_an = year-centered number of foraging dives per hour)
t1 <- by_p5 %>%                                       
  group_by(season, period) %>%                         
  summarise_at(vars(F_dives_PH),              
               list(av_F_dives_PH = mean))
t1 %>%
  filter(period!="Incubation"&period!="Chick rearing") %>%
  group_by(season) %>%
  summarize(mean_F_dives_PH=mean(av_F_dives_PH))
# A tibble: 3 x 2
#      season mean_F_dives_PH
# <int>           <dbl>
#  1   2016            2.93
#  2   2017            2.42
#  3   2018            2.80
#t1 %>%
#  filter(period!="Incubation"&period!="Chick rearing") %>%
#  group_by(season) %>%
#  summarize(sd_F_dives_PH=sd(av_F_dives_PH))
# A tibble: 3 x 2
# season sd_F_dives_PH
# <int>         <dbl>
#  1   2016          1.02
#  2   2017          1.32
#  3   2018          1.92
by_p5$F_dives_PH_an <- ifelse(by_p5$season=="2016", (by_p5$F_dives_PH-2.93),
                              ifelse(by_p5$season=="2017", (by_p5$F_dives_PH-2.42),
                                     ifelse(by_p5$season=="2018", (by_p5$F_dives_PH-2.80),
                                            NA)))

## Model selection on overall model with lm (n=18)
o1 <- lm(F_dives_PH_an ~ period*br_st + age_class*deploy_num + sex + colony, data=by_p5, na.action=na.fail)
#plot(o1)
select1 <- dredge(o1, rank="AICc")
subset(select1, delta <2) # show the top models within 2 AICc (6 models)
coefTable(select1)[1] # get the summary for the best model 
b1 <- lm(F_dives_PH_an ~ period*br_st + age_class*deploy_num + sex, data=by_p5, na.action=na.fail)
summary(b1) # Adj R²=0.61
-2*logLik(b1) # 718.721 (df=25)
# Plot predictions with a colorblind-friendly palette:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pr <- ggpredict(b1, c("period","br_st"))
plot(pr) + 
  labs(x = "Period", y = "Centred number of foraging dives per hour", title=NULL, colour="Breeding status") + 
  scale_colour_manual(values=cbbPalette) +
  theme_classic()
pr <- ggpredict(b1, c("age_class","deploy_num"))
plot(pr) + 
  labs(x = "Age class", y = "Centred number of foraging dives per hour", title=NULL, colour="Time since deployment") + 
  scale_colour_manual(values=cbbPalette) +
  theme_classic()
confint(b1) 
coefTable(select1)[2]
b2 <- lm(F_dives_PH_an ~ period*br_st + age_class*deploy_num + sex + colony, data=by_p5, na.action=na.fail)
summary(b2) # Adj R²=0.61
-2*logLik(b2) # 716.2552 (df=26)
confint(b2) 
b3 <- lm(F_dives_PH_an ~ period*br_st + sex, data=by_p5, na.action=na.fail)
summary(b3) # Adj R²=0.59
-2*logLik(b3) # 731.302 (df=20)
confint(b3) 
b4 <- lm(F_dives_PH_an ~ period*br_st + age_class + sex, data=by_p5, na.action=na.fail)
summary(b4) # Adj R²=0.60
-2*logLik(b4) #  726.7226 (df=22)
confint(b4) 
b5 <- lm(F_dives_PH_an ~ period*br_st + colony + sex, data=by_p5, na.action=na.fail)
summary(b5) # Adj R²=0.60
-2*logLik(b5) # 729.2059 (df=21)
confint(b5)
b6 <- lm(F_dives_PH_an ~ period*br_st + age_class + colony + sex, data=by_p5, na.action=na.fail)
summary(b6) # Adj R²=0.60
-2*logLik(b6) # 724.4806 (df=23)
confint(b6) 
b7 <- lm(F_dives_PH_an ~ period*br_st + deploy_num + sex, data=by_p5, na.action=na.fail)
summary(b7) # Adj R²=0.59
-2*logLik(b7) # 731.1984 (df=21)
b8 <- lm(F_dives_PH_an ~ period*br_st + age_class + deploy_num + sex, data=by_p5, na.action=na.fail)
summary(b8) # Adj R²=0.60
-2*logLik(b8) # 726.4116 (df=23)
b9 <- lm(F_dives_PH_an ~ period*br_st + age_class + deploy_num + sex + colony, data=by_p5, na.action=na.fail)
summary(b9) # Adj R²=0.60
-2*logLik(b9) # 724.0742 (df=24)
b10 <- lm(F_dives_PH_an ~ period*br_st + deploy_num + sex + colony, data=by_p5, na.action=na.fail)
summary(b10) # Adj R²=0.59
-2*logLik(b10) # 729.0515 (df=22)
b11 <- lm(F_dives_PH_an ~ 1, data=by_p5, na.action=na.fail)
summary(b11) # Adj R²=0.00
-2*logLik(b11) # 965.3953 (df=2)

### Variation in maximum foraging depth
## Calculating year-centered maximum foraging depth (av_maxdep = maximum foraging depth averaged over periods of the annual cycle in m, 
## av_maxdep_an = year-centered nav_maxdep)
t3 <- by_p5 %>%                                       
  group_by(season, period) %>%                         
  summarise_at(vars(av_maxdep),              
               list(av_maxdep = mean))
t3 %>%
  filter(period!="Incubation" & period!="Chick rearing") %>%
  group_by(season) %>%
  summarize(mean_av_maxdep=mean(av_maxdep))
# A tibble: 3 x 2
#      season mean_F_dives_PH
# <int>           <dbl>
#  1   2016            42.3
#  2   2017            36.5
#  3   2018            39.5
#t3 %>%
#  filter(period!="Incubation"&period!="Chick rearing") %>%
#  group_by(season) %>%
#  summarize(sd_av_maxdep=sd(av_maxdep))
by_p5$av_maxdep_an <- ifelse(by_p5$season=="2016", (by_p5$av_maxdep-42.3),
                             ifelse(by_p5$season=="2017", (by_p5$av_maxdep-36.5),
                                    ifelse(by_p5$season=="2018", (by_p5$av_maxdep-39.5),
                                           NA)))

## Model selection on overall model with lm (n=18)
o3 <- lm(av_maxdep_an ~ period*br_st + age_class*deploy_num + colony + sex, data=by_p5, na.action=na.fail)
selecto3 <- dredge(o3, rank="AICc")
subset(selecto3, delta <2) # show the top models within 2 AICc (3 models)
coefTable(selecto3)[1] 
o3_1 <- lm(av_maxdep_an ~ period*br_st + colony, data=by_p5, na.action=na.fail)
summary(o3_1) # Adj R²=0.55
-2*logLik(o3_1) # 1754.616 (df=20)
# Plot predictions with a colorblind-friendly palette:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pr <- ggpredict(o3_1, c("period","br_st"))
plot(pr) + 
  labs(x = "Period", y = "Centred maximum foraging depth (m)", title=NULL, colour="Breeding status") + 
  scale_colour_manual(values=cbbPalette) +
  theme_classic()
confint(o3_1)
coefTable(selecto3)[2] 
o3_2 <- lm(av_maxdep_an ~ period*br_st + colony + sex, data=by_p5, na.action=na.fail)
summary(o3_2) # Adj R²=0.55
-2*logLik(o3_2) # 1753.732 (df=21)
confint(o3_2)
coefTable(selecto3)[3]
o3_3 <- lm(av_maxdep_an ~ period*br_st, data=by_p5, na.action=na.fail)
summary(o3_3) # Adj R²=0.54
-2*logLik(o3_3) # 1758.596 (df=19)
confint(o3_3)
coefTable(selecto3)[4] 
o3_4 <- lm(av_maxdep_an ~ period*br_st + age_class + colony, data=by_p5, na.action=na.fail)
summary(o3_4) # Adj R²=0.55
-2*logLik(o3_4) # 1751.886 (df=22)
o3_5 <- lm(av_maxdep_an ~ period*br_st + deploy_num + colony, data=by_p5, na.action=na.fail)
summary(o3_5) # Adj R²=0.55
-2*logLik(o3_5) # 1754.604 (df=21)
o3_6 <- lm(av_maxdep_an ~ period*br_st + sex, data=by_p5, na.action=na.fail)
summary(o3_6) # Adj R²=0.54
-2*logLik(o3_6) # 1757.951 (df=20)
o3_7 <- lm(av_maxdep_an ~ period*br_st + age_class, data=by_p5, na.action=na.fail)
summary(o3_7) # Adj R²=0.54
-2*logLik(o3_7) # 1755.865 (df=21)
o3_8 <- lm(av_maxdep_an ~ period*br_st + age_class + sex + colony, data=by_p5, na.action=na.fail)
summary(o3_8) # Adj R²=0.55
-2*logLik(o3_8) # 1750.999 (df=23)
o3_9 <- lm(av_maxdep_an ~ period*br_st + deploy_num + sex + colony, data=by_p5, na.action=na.fail)
summary(o3_9) # Adj R²=0.54
-2*logLik(o3_9) # 1753.721 (df=22)
o3_10 <- lm(av_maxdep_an ~ period*br_st + deploy_num, data=by_p5, na.action=na.fail)
summary(o3_10) # Adj R²=0.54
-2*logLik(o3_10) # 1758.596 (df=20)
o3_11 <- lm(av_maxdep_an ~ period*br_st + age_class + deploy_num + colony, data=by_p5, na.action=na.fail)
summary(o3_11) # Adj R²=0.55
-2*logLik(o3_11) # 1751.885 (df=23)
o3_12 <- lm(av_maxdep_an ~ period*br_st + age_class + sex, data=by_p5, na.action=na.fail)
summary(o3_12) # Adj R²=0.54
-2*logLik(o3_12) # 1755.176 (df=22)
o3_13 <- lm(av_maxdep_an ~ period*br_st + deploy_num + sex, data=by_p5, na.action=na.fail)
summary(o3_13) # Adj R²=0.54
-2*logLik(o3_13) # 1757.951 (df=21)
o3_14 <- lm(av_maxdep_an ~ period*br_st + age_class + deploy_num, data=by_p5, na.action=na.fail)
summary(o3_14) # Adj R²=0.54
-2*logLik(o3_14) # 1755.859 (df=22)
o3_15 <- lm(av_maxdep_an ~ period*br_st + age_class + deploy_num + sex + colony, data=by_p5, na.action=na.fail)
summary(o3_15) # Adj R²=0.55
-2*logLik(o3_15) # 1750.998 (df=24)
o3_16 <- lm(av_maxdep_an ~ 1, data=by_p5, na.action=na.fail)
summary(o3_16) # Adj R²=0.55
-2*logLik(o3_16) # 1962.745 (df=2)


### Load data about all individual foraging dives
### Download from: https://drive.google.com/file/d/106n19QGIYdPLdzH9btHToEc6P1hKgj_o/view?usp=sharing
dsfC_F <- read.csv (".../dsfC_F.csv")

## Assign factor levels
dsfC_F$age_class <- factor(dsfC_F$age_class, levels=c("Young", "Middle-age", "Old"))
dsfC_F$period <- factor(dsfC_F$period, levels = c("Incubation", "Chick rearing", "Moult", "Outgoing migration", "Migration max", "Returning migration"))
levels(dsfC_F$period) <- c("Incubation", "Chick rearing", "Moult", "Outgoing migration", "Migration max", "Returning migration")
dsfC_F$br_st <- factor(dsfC_F$br_st, levels=c("Non breeder", "Failed breeder", "Successful breeder"))
dsfC_F$sex <- factor(dsfC_F$sex, levels=c("F", "M"))

### Variation in descent rate of individual foraging dives
## Calculating year-centered descent rate (descr = descent rate in m/s, 
## descr_an = year-centered descent rate)
t2 <- dsfC_F %>%
  filter(period!="Incubation" & period!="Chick rearing" & descr!="NA") %>%
  group_by(season) %>%
  summarize(mean_descr=mean(descr), sd_descr=sd(descr))
t2
#    season mean_descr sd_descr
# <int>      <dbl>    <dbl>
#  1   2016       1.39    0.426
#  2   2017       1.41    0.428
#  3   2018       1.44    0.440
dsfC_F$ctr.descr <- ifelse(dsfC_F$season=="2016", (dsfC_F$descr-1.39),
                           ifelse(dsfC_F$season=="2017", (dsfC_F$descr-1.41),
                                  (dsfC_F$descr-1.44)))

## Model selection on overall model with lm (n=18)
m4 <- lmer(ctr.descr~maxdep + sex + age_class*deploy_num + period*br_st + (1|bird_id), data=dsfC_F, REML=FALSE, na.action=na.fail)
selectm4 <- dredge(m4, rank="AICc")
subset(selectm4, delta <2) ## 1 top model
# age_class*deploy_num + maxdep + sex + period*br_st
library(performance)
r2(m4) # Conditional 0.25, marginal 0.23
m4_reml <- lmer(ctr.descr~maxdep + sex + age_class*deploy_num + period*br_st + (1|bird_id), data=dsfC_F, REML=TRUE, na.action=na.fail)
r2(m4_reml) # Conditional 0.25, marginal 0.23
summary(m4_reml)
-2*logLik(m4_reml) # 969025 (df=27)
pr <- ggpredict(m4_reml, c("age_class", "deploy_num"))
# Plot predictions with a colorblind-friendly palette:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot(pr) + 
  labs(x = "Age class", y = "Centred descent rate (m/s)", title=NULL, colour="Time since attachment") + 
  scale_colour_manual(values=cbbPalette) +
  theme_classic()
pr <- ggpredict(m4_reml, c("period","br_st"))
plot(pr) + 
  labs(x = "Period", y = "Centred descent rate (m/s)", title=NULL, colour="Breeding status") + 
  scale_colour_manual(values=cbbPalette) +
  theme_classic()
confint(m4_reml)
m4_2_reml <- lmer(ctr.descr~maxdep + age_class*deploy_num + period*br_st + (1|bird_id), data=dsfC_F, REML=TRUE, na.action=na.fail)
r2(m4_2_reml) # Conditional 0.25, marginal 0.21
-2*logLik(m4_2_reml) # 969026.2 (df=26)
m4_i_reml <- lmer(ctr.descr~1 + (1|bird_id), data=dsfC_F, REML=TRUE, na.action=na.fail)
-2*logLik(m4_i_reml) # 1233016 (df=3)

### End of code