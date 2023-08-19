# read libraries
pacman::p_load(ape, here, metafor, tidyverse)
source("R/function_2.R")

# get data
dat_pred <- read_csv(here("data/predator_22072023.csv"))
dim(dat_pred)

# calculate lnRR and lnRR_var
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

# test 1 - loop through all trees and find minimum AIC value tree
ma_pred_test <- lapply(phylo_vcv, function(phylo_vcv) {
  rma.mv(yi = lnRR,
                         V = lnRR_var,
                         random = list(~1 | Study_ID,
                                       ~1 | Cohort_ID, 
                                       ~1 | Shared_control_ID,
                                       ~1 | Obs_ID,
                                       ~1 | Bird_species,
                                       ~1 | Phylogeny
                                       ),
                         R = list(Phylogeny = phylo_vcv), 
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

#    logLik   Deviance        AIC        BIC       AICc   
# -155.1741   310.3482   324.3482   343.6233   325.3852   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.1209  0.3477     33     no          Cohort_ID   no 
# sigma^2.3  0.0830  0.2882     29     no  Shared_control_ID   no 
# sigma^2.4  0.5344  0.7310    117     no             Obs_ID   no 
# sigma^2.5  0.0000  0.0000      7     no       Bird_species   no 
# sigma^2.6  0.0000  0.0000      7     no          Phylogeny  yes 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139  

i2_ml(best_model_result)
#             I2_Total          I2_Study_ID         I2_Cohort_ID 
#         9.960611e+01         4.279216e-07         1.631067e+01 
# I2_Shared_control_ID            I2_Obs_ID      I2_Bird_species 
#         1.120302e+01         7.209242e+01         8.590309e-10 
#         I2_Phylogeny 
#         9.732004e-10 

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
                                       ~1 | Obs_ID,
                                       ~1 | Bird_species,
                                       ~1 | Phylogeny
                                       ),
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
# sigma^2.4  0.5344  0.7310    117     no             Obs_ID   no 
# sigma^2.5  0.0000  0.0000      7     no       Bird_species   no 
# sigma^2.6  0.0000  0.0000      7     no          Phylogeny  yes 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139    

i2_ml(ma_pred_test1)
#             I2_Total          I2_Study_ID         I2_Cohort_ID 
#         9.960611e+01         4.356581e-07         1.631066e+01 
# I2_Shared_control_ID            I2_Obs_ID      I2_Bird_species 
#         1.120303e+01         7.209242e+01         9.850812e-10 
#         I2_Phylogeny 
#         8.926356e-10 


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
                                       ~1 | Obs_ID,
                                       ~1 | Bird_species,
                                       ~1 | Phylogeny
                                       ),
                         R = list(Phylogeny = phylo_vcv2), 
                         test = "t",
                         method = "REML",
                         sparse = TRUE,
                         data = dat)

summary(ma_pred_test2)
# Multivariate Meta-Analysis Model (k = 117; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -155.1741   310.3482   324.3482   343.6233   325.3852   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.1209  0.3477     33     no          Cohort_ID   no 
# sigma^2.3  0.0830  0.2882     29     no  Shared_control_ID   no 
# sigma^2.4  0.5344  0.7310    117     no             Obs_ID   no 
# sigma^2.5  0.0000  0.0000      7     no       Bird_species   no 
# sigma^2.6  0.0000  0.0000      7     no          Phylogeny  yes 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139  

i2_ml(ma_pred_test2)
#             I2_Total          I2_Study_ID         I2_Cohort_ID 
#         9.960611e+01         4.658025e-07         1.631064e+01 
# I2_Shared_control_ID            I2_Obs_ID      I2_Bird_species  
#         1.120305e+01         7.209241e+01         3.743351e-10 
#         I2_Phylogeny 
#         1.527690e-09  

# test 4 - use vcalc()
# ASK Shinichi

model_phy <- NULL
phylo_vcv3 <- NULL

class(model_phy) <- "list"

V <- vcalc(lnRR_var, cluster = Phylogeny, subgroup = Obs_ID, data = dat)

for (i in 1:50) {
  # Calculate the VCV matrix for the current tree
  vcv_matrix <- vcv(trees[[i]]) 
  # Store the result in the list
  phylo_vcv3[[i]] <- vcv_matrix

# meta-analysis
 ma_pred_test3 <- rma.mv(yi = lnRR,
                          V = V,
                          random = list(~1 | Study_ID,
                                        ~1 | Cohort_ID, 
                                        ~1 | Shared_control_ID,
                                        ~1 | Obs_ID,
                                        ~1 | Bird_species,
                                        ~1 | Phylogeny
                                        ),
                          R = list(Phylogeny = phylo_vcv3[[i]]), 
                          test = "t",
                          method = "REML",
                          sparse = TRUE,
                          data = dat)
   
 model_phy[[i]] <- ma_pred_test3
 print(i)

}

summary(model_phy[[1]])
# Multivariate Meta-Analysis Model (k = 117; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -155.1741   310.3482   324.3482   343.6233   325.3852   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.1209  0.3477     33     no          Cohort_ID   no 
# sigma^2.3  0.0830  0.2882     29     no  Shared_control_ID   no 
# sigma^2.4  0.5344  0.7310    117     no             Obs_ID   no 
# sigma^2.5  0.0000  0.0000      7     no       Bird_species   no 
# sigma^2.6  0.0000  0.0000      7     no          Phylogeny  yes 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139  

i2_ml(model_phy[[1]])
#         I2_Total          I2_Study_ID         I2_Cohort_ID 
#         9.960611e+01         4.296546e-07         1.631066e+01 
# I2_Shared_control_ID            I2_Obs_ID      I2_Bird_species 
#         1.120302e+01         7.209242e+01         8.558411e-10 
#         I2_Phylogeny 
#         1.029748e-09 

all.m.phy <- model_phy[[1]]$sigma2

for (i in 2:50) {
  all.m.phy <- rbind(all.m.phy, model_phy[[i]]$sigma2)
}
all.m.phy <- data.frame(all.m.phy, byrow = TRUE)
col_names <- c("sigma^2.1_Study_ID", "sigma^2.2_Cohort_ID", "sigma^2.3_SharedControl_ID", "sigma^2.4_Obs_ID", "sigma^2.5_BirdSpecies", "sigma^2.6_Phylogeny")
colnames(all.m.phy) <- col_names
print(all.m.phy)