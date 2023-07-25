# reading datasets

library(tidyverse)
library(here)
library(orchaRd)
library(metafor)
library(ape)

# get data
dat_prey <- read_csv(here("data/prey_22072023.csv"))
dat_pred <- read_csv(here("data/predator_22072023.csv"))
dat_prey <- dat_prey[1:146,] # exclude unclear report
dim(dat_prey)
dim(dat_pred)

# get phylogeny
# TODO how to get phylo and vcv from multiple phylogeny
tree <- read.nexus(here("data/bird_phy.nex"))
tree_1 <-  tree[[1]]
plot(tree[[1]])
tree_1 <- compute.brlen(tree_1)
cor <- vcv(tree_1, cor = T)

##########
# prey
##########

# turn all character strings to factor
dat_prey <- dat_prey %>%
  mutate_if(is.character, as.factor)

summary(dat_prey)

dat1 <- effect_lnRR(dat_prey)

hist(dat1$lnRR) 
hist(dat1$lnRR_var)

#which(dat1$lnRR == min(dat1$lnRR))

# meta-analysis
# TODO - get phylo, VCV and other things
dat1$Shared_control_ID <- 1:nrow(dat1)

ma_prey <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       random = list(~1 | Study_ID,
                     ~1 | Cohort_ID,
                     ~1 | Shared_control_ID),
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
# eyespot or conspicuous?
mr_prey <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Treatment_stimulus,
       random = list(~1 | Study_ID,
                     ~1 | Cohort_ID,
                     ~1 | Shared_control_ID),
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

# size of pattern
mr_prey1 <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Diameter_pattern,
       random = list(~1 | Study_ID,
                     ~1 | Cohort_ID,
                     ~1 | Shared_control_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat1)

summary(mr_prey1)

bubble_prey1 <-  orchaRd::mod_results(mr_prey1, mod = "Diameter_pattern", group = "Study_ID")
bubble_plot(mr_prey1,
             mod = "Diameter_pattern",
             group = "Study_ID",
             xlab = "size (mm)")

#FIXME - Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
#replacement has 0 rows, data has 146

# area of pattern
mr_prey2 <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Area_pattern,
       random = list(~1 | Study_ID,
                     ~1 | Cohort_ID, 
                     ~1 | Shared_control_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat1)

summary(mr_prey2)

bubble_plot(mr_prey2,
             mod = "Area_pattern",
             group = "Study_ID",
             xlab = "Area (mm2)")

# number of pattern
mr_prey3 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Number_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat1)

summary(mr_prey3)

bubble_plot(mr_prey3,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number")

# type of prey
mr_prey4 <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Type_prey,
       random = list(~1 | Study_ID,
                     ~1 | Cohort_ID, 
                     ~1 | Shared_control_ID),
       #R = list(Phylo = cor_tree),
       test = "t",
       method = "REML", 
       sparse = TRUE,
       data = dat1)
       
summary(mr_prey4)

orchard_plot(mr_prey4,
             mod = "Type_prey",
             group = "Study_ID",
             xlab = "Type of prey")

# shape of prey
mr_prey5 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_prey,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat1)

summary(mr_prey5)

orchard_plot(mr_prey5,
             mod = "Shape_prey",
             group = "Study_ID",
             xlab = "Shape of prey")

# shape of pattern
mr_prey6 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat1)

summary(mr_prey6)

orchard_plot(mr_prey6,
             mod = "Shape_pattern",
             group = "Study_ID",
             xlab = "Shape of pattern")

##########
# predator
##########
# turn all character strings to factor
dat_pred <- dat_pred %>%
  mutate_if(is.character, as.factor)

summary(dat_pred)

dat2 <- effect_lnRR(dat_pred)

hist(dat2$lnRR) 
hist(dat2$lnRR_var)

# meta-analysis
dat2$Shared_control_ID <- 1:nrow(dat2)

ma_pred <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       random = list(~1 | Study_ID,
                     ~1 | Cohort_ID, 
                     ~1 | Shared_control_ID),
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
# eyespot or conspicuous?
mr_pred <- rma.mv(yi = lnRR,
       V = lnRR_var, 
       mods = ~ Treatment_stimulus,
       random = list(~1 | Study_ID,
                     ~1 | Cohort_ID, 
                     ~1 | Shared_control_ID),
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

# size of pattern
mr_pred1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Diameter_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID,
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred1)


bubble_plot(mr_pred1,
            mod = "Diameter_pattern",
            group = "Study_ID",
            xlab = "size (mm)")

# area of pattern
mr_pred2 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred2)

bubble_plot(mr_pred2,
            mod = "Area_pattern",
            group = "Study_ID",
            xlab = "Area (mm2)")


# number of pattern
mr_pred3 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Number_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred3)

bubble_plot(mr_pred3,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number")

# type of prey
mr_pred4 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Type_prey,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred4)

orchard_plot(mr_pred4,
             mod = "Type_prey",
             group = "Study_ID",
             xlab = "Type of prey")

# shape of prey
mr_pred5 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_prey,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred5)

orchard_plot(mr_pred5,
             mod = "Shape_prey",
             group = "Study_ID",
             xlab = "Shape of prey")

# shape of pattern
mr_pred6 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   #R = list(Phylo = cor_tree),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred6)

orchard_plot(mr_pred6,
             mod = "Shape_pattern",
             group = "Study_ID",
             xlab = "Shape of pattern")

