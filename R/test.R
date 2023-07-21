# reading datasets

library(tidyverse)
library(here)
library(orchaRd)
library(metafor)

# get data
dat_prey <- read_csv(here("data/prey_19072023.csv"))
dat_pred <- read_csv(here("data/predator_19072023.csv"))
dim(dat_prey)
dim(dat_pred)

##########
# prey
##########

# turn all character strings to factor
dat_prey <- dat_prey %>%
  mutate_if(is.character, as.factor)

summary(dat_prey)

# check whether "T_proportion" and "C_proportion" have 0 or not
dat_prey %>%
  filter(T_proportion == 0 | C_proportion == 0) %>%
  select(T_proportion, C_proportion)
# # A tibble: 3 × 2
#   T_proportion C_proportion
#          <dbl>        <dbl>
# 1        0.184            0
# 2        0.318            0
# 3        0.660            0

# check the samllest values in "C_proportion" apart from 0
dat_prey %>%
  filter(C_proportion != 0) %>%
  arrange(C_proportion) %>%
  head(3) %>%
  select(C_proportion)

#   # A tibble: 3 × 1
#   C_proportion
#          <dbl>
# 1       0.0107
# 2       0.0107
# 3       0.0183

# given this we should repliace 0 in "C_proportion" with 0.01
dat_prey <- dat_prey %>%
  mutate(C_proportion = ifelse(C_proportion == 0, 0.01, C_proportion))

# check whether "T_meam" and "C_mean" have 0 or not 
dat_prey %>%
  filter(T_mean == 0 | C_mean == 0) %>%
  select(T_mean, C_mean)

# check whether "T_sd" and "C_sd" have 0 or not
dat_prey %>%
  filter(T_sd == 0 | C_sd == 0) %>%
  select(T_sd, C_sd)


dat1 <- effect_lnRR_prey(dat_prey)

hist(dat1$lnRR) 

# TODO something
hist(dat1$lnRR_var)

# this looks strange
# TODO check 83rd row
which(dat1$lnRR == min(dat$lnRR))

hist(dat1$lnRR[-83])
hist(log(dat1$lnRR_var[-83]))

# meta-analysis
# TODO - get phylo, VCV and other things
dat1$Obs_ID <- 1:nrow(dat1)
dat1 <- dat1[-83, ]

ma_prey <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       random = list(~1 | Study_ID,
                     #~1 | Cohort_ID, #TODO  some missing values
                     ~1 | Obs_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat1)

summary(ma_prey)
i2_ml(ma_prey)


orchard_plot(ma_prey,
             group = "Study_ID",
             xlab = "log response ratio (lnRR)")

# meta-regression
mr_prey <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Treatment_stimulus,
       random = list(~1 | Study_ID,
                     #~1 | Cohort_ID, #TODO  some missing values
                     ~1 | Obs_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat1)

summary(mr_prey)


orchard_plot(mr_prey,
             mod = "Treatment_stimulus",
             group = "Study_ID",
             xlab = "log response ratio (lnRR)")

# eyespot size

mr_prey1 <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Diameter_eyespot,
       random = list(~1 | Study_ID,
                     #~1 | Cohort_ID, #TODO  some missing values
                     ~1 | Obs_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat1)

summary(mr_prey1)


bubble_plot(mr_prey1,
             mod = "Diameter_eyespot",
             group = "Study_ID",
             xlab = "Eyespots size (mm)")

# number of eyespots
mr_prey2 <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Number_eyespot,
       random = list(~1 | Study_ID,
                     #~1 | Cohort_ID, #TODO  some missing values
                     ~1 | Obs_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat1)

summary(mr_prey2)

bubble_plot(mr_prey2,
             mod = "Number_eyespot",
             group = "Study_ID",
             xlab = "Eyespot number")

# Type_of_prey
mr_prey3 <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Type_of_prey,
       random = list(~1 | Study_ID,
                     #~1 | Cohort_ID, #TODO  some missing values
                     ~1 | Obs_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat1)
       
summary(mr_prey3)

orchard_plot(mr_prey3,
             mod = "Type_of_prey",
             group = "Study_ID",
             xlab = "Type of prey")

##########
# predator
##########
# turn all character strings to factor
dat_pred <- dat_pred %>%
  mutate_if(is.character, as.factor)

summary(dat_pred)

# check whether "T_mean" and "C_mean" have 0 or not
dat_pred %>%
  filter(T_mean == 0 | C_mean == 0) %>%
  select(T_mean, C_mean)

# repliace 0 in "C_mean" with 0.01
# TODO - check whether this is correct
dat_pred <- dat_pred %>%
  mutate(C_mean = ifelse(C_mean == 0, 0.01, C_mean))

dat2 <- effect_lnRR_predator(dat_pred)

hist(dat2$lnRR) 
hist(dat2$lnRR_var)

# meta-analysis
dat2$Obs_ID <- 1:nrow(dat2)

ma_pred <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       random = list(~1 | Study_ID,
                     #~1 | Cohort_ID, #TODO  some missing values
                     ~1 | Obs_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat2)

summary(ma_pred)
i2_ml(ma_pred)


orchard_plot(ma_pred,
             group = "Study_ID",
             xlab = "log response ratio (lnRR)")

# meta-regression
# TODO Treatment_stimulus - mistake in there
# fixing Treatment_stimulus - merge "eye spot" and "eyespot"
dat2$Treatment_stimulus <- 
  ifelse(dat2$Treatment_stimulus == "eye spots", "eyespots", dat2$Treatment_stimulus)

mr_pred <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Treatment_stimulus,
       random = list(~1 | Study_ID,
                     #~1 | Cohort_ID, #TODO  some missing values
                     ~1 | Obs_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat2)

summary(mr_pred)


orchard_plot(mr_pred,
             mod = "Treatment_stimulus",
             group = "Study_ID",
             xlab = "log response ratio (lnRR)")