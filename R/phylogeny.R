# read libraries
library(ape)
library(clubSandwich)
library(here)
library(metafor)
library(orchaRd)
library(phangorn)
library(tidyverse)

# get data
dat_pred <- read_csv(here("data/predator_22072023.csv"))
dim(dat_pred)

# calculate lnRR and lnRR_var
source("R/function_2.R")
dat <- effect_lnRR(dat_pred)
dat$Obs_ID <- 1:nrow(dat)

# phylogenetic tree
trees <- read.nexus(here("data/bird_phy.nex"))
phylo_vcv <- vector("list", 50)

# Loop over each tree in the list
for (i in 1:50) {
  # Calculate the VCV matrix for the current tree
  vcv_matrix <- vcv(trees[[i]])
  
  # Store the result in the list
  phylo_vcv[[i]] <- vcv_matrix
}

phylo_vcv

# test 1 - loop through all trees and find minimum AIC value tree
ma_pred_test <- lapply(phylo_vcv, function(phylo_vcv) {
  rma.mv(yi = lnRR,
         V = lnRR_var,
         random = list(~1 | Study_ID,
                       ~1 | Cohort_ID, 
                       ~1 | Shared_control_ID,
                       ~1 | Bird_species,
                       ~1 | Phylogeny,
                       ~1 |Obs_ID),
         R = list(Bird_species = phylo_vcv), 
         test = "t",
         method = "REML",
         sparse = TRUE,
         data = dat)
})

aic_values <- sapply(ma_pred_test, AIC)
min_aic <- min(aic_values)

# get the index and result of the model with the smallest AIC value
best_model_index <- which(aic_values == min_aic)
best_model_result <- ma_pred_test[[best_model_index]]
print(best_model_result)

# Multivariate Meta-Analysis Model (k = 117; method: REML)
# Variance Components:
#           estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.1209  0.3477     33     no          Cohort_ID   no 
# sigma^2.3  0.0830  0.2882     29     no  Shared_control_ID   no 
# sigma^2.4  0.0000  0.0000      7     no       Bird_species  yes 
# sigma^2.5  0.0000  0.0000      7     no          Phylogeny   no 
#sigma^2.6  0.5344  0.7310    117     no             Obs_ID    no 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139    

i2_ml(ma_pred_test[[37]])
#        I2_Total             I2_Study_ID         I2_Cohort_ID 
#        9.960611e+01         4.286317e-07          1.631045e+01 
#    I2_Shared_control_ID    I2_Bird_species      I2_Phylogeny 
#       1.120321e+01         1.054767e-09         8.383898e-10 
#           I2_Obs_ID 
#        7.209245e+01 

# test 2 - average the phylogenetic covariance matrix and use it in the model
phylo_vcv1 <- Reduce("+", phylo_vcv) / length(phylo_vcv)

#                     Gallus_gallus Cyanocitta_cristata Sturnus_vulgaris
# Gallus_gallus                   1           0.0000000        0.0000000
# Cyanocitta_cristata             0           1.0000000        0.4958565
# Sturnus_vulgaris                0           0.4958565        1.0000000
# Ficedula_hypoleuca              0           0.4958565        0.7022586
# Emberiza_sulphurata             0           0.4958565        0.5620193
# Parus_major                     0           0.4958565        0.5533103
# Parus_caeruleus                 0           0.4958565        0.5533103
#                     Ficedula_hypoleuca Emberiza_sulphurata Parus_major
# Gallus_gallus                0.0000000           0.0000000   0.0000000
# Cyanocitta_cristata          0.4958565           0.4958565   0.4958565
# Sturnus_vulgaris             0.7022586           0.5620193   0.5533103
# Ficedula_hypoleuca           1.0000000           0.5620193   0.5533103
# Emberiza_sulphurata          0.5620193           1.0000000   0.5533103
# Parus_major                  0.5533103           0.5533103   1.0000000
# Parus_caeruleus              0.5533103           0.5533103   0.8124203
#                    Parus_caeruleus
# Gallus_gallus             0.0000000
# Cyanocitta_cristata       0.4958565
# Sturnus_vulgaris          0.5533103
# Ficedula_hypoleuca        0.5533103
# Emberiza_sulphurata       0.5533103
# Parus_major               0.8124203
# Parus_caeruleus           1.0000000

ma_pred_test1 <- rma.mv(yi = lnRR,
                         V = lnRR_var,
                         random = list(~1 | Study_ID,
                                       ~1 | Cohort_ID, 
                                       ~1 | Shared_control_ID,
                                       ~1 | Bird_species,
                                       ~1 | Phylogeny,
                                       ~1 | Obs_ID),
                         R = list(Bird_species = phylo_vcv1), 
                         test = "t",
                         method = "REML",
                         sparse = TRUE,
                         data = dat)

summary(ma_pred_test1)

# Multivariate Meta-Analysis Model (k = 117; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -155.1741   310.3482   324.3482   343.6233   325.3852   

# Variance Components:
#             estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.1209  0.3477     33     no          Cohort_ID   no 
# sigma^2.3  0.0830  0.2882     29     no  Shared_control_ID   no 
# sigma^2.4  0.0000  0.0000      7     no       Bird_species  yes 
# sigma^2.5  0.0000  0.0000      7     no          Phylogeny  yes 
# sigma^2.6  0.5344  0.7310    117     no             Obs_ID   no 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139    

i2_ml(ma_pred_test1)
#             I2_Total          I2_Study_ID         I2_Cohort_ID 
#         9.960611e+01         5.619035e-07         1.631044e+01 
# I2_Shared_control_ID      I2_Bird_species         I2_Phylogeny 
#         1.120322e+01         2.087399e-09         2.087399e-09 
#            I2_Obs_ID 
#         7.209245e+01 

# test 3 - use the maximum clade credibility tree
# Use 1,000 trees downloaded from www.birdtree.org based on predator dataset, 
# compute the maximum clade credibility tree, compute branch lengths, compute the correlation matrix
tips <- dat_pred$Phylogeny
pruned.trees <- lapply(trees, keep.tip, tip=tips)
class(pruned.trees) <- "multiPhylo" 

mcc_pruned <- maxCladeCred(pruned.trees)
mcc_pruned <- ladderize(mcc_pruned)
mcc_pruned$tip.label
# [1] "Gallus_gallus"       "Cyanocitta_cristata" "Sturnus_vulgaris"   
# [4] "Ficedula_hypoleuca"  "Emberiza_sulphurata" "Parus_major"        
# [7] "Parus_caeruleus" 

mcc_ult <- compute.brlen(mcc_pruned, power = 1)
phylo_vcv2 <- vcv(mcc_ult, cor=T)

# meta-analysis
ma_pred_test2 <- rma.mv(yi = lnRR,
                         V = lnRR_var,
                         random = list(~1 | Study_ID,
                                       ~1 | Cohort_ID, 
                                       ~1 | Shared_control_ID,
                                       ~1 | Bird_species,
                                       ~1 | Phylogeny,
                                       ~1 | Obs_ID),
                         R = list(Bird_species = phylo_vcv2), 
                         test = "t",
                         method = "REML",
                         sparse = TRUE,
                         data = dat)

summary(ma_pred_test2)
#    logLik   Deviance        AIC        BIC       AICc   
# -155.1741   310.3482   324.3482   343.6233   325.3852   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.1209  0.3477     33     no          Cohort_ID   no 
# sigma^2.3  0.0830  0.2882     29     no  Shared_control_ID   no 
# sigma^2.4  0.0000  0.0000      7     no       Bird_species  yes 
# sigma^2.5  0.0000  0.0000      7     no          Phylogeny   no 
# sigma^2.6  0.5344  0.7310    117     no             Obs_ID   no 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139    

i2_ml(ma_pred_test2)
#             I2_Total          I2_Study_ID         I2_Cohort_ID 
#         9.960611e+01         4.658025e-07         1.631064e+01 
# I2_Shared_control_ID      I2_Bird_species         I2_Phylogeny 
#         1.120305e+01         3.760417e-10         1.534219e-09 
#            I2_Obs_ID 
#         7.209242e+01 