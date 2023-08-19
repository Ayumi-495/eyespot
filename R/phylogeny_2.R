# read libraries
pacman::p_load(ape, here, tidyverse)
source("R/function_2.R")

# data and phylogeny
dat_pred <- read_csv(here("data/predator_22072023.csv"))
trees <- read.nexus(here("data/bird_phy.nex"))
is.ultrametric(trees[[1]])
# [1] TRUE

vcv(trees[[1]], cor = T) # example of vcv matrix
#                     Gallus_gallus Cyanocitta_cristata Sturnus_vulgaris
# Gallus_gallus                   1           0.0000000        0.0000000
# Cyanocitta_cristata             0           1.0000000        0.4437522
# Sturnus_vulgaris                0           0.4437522        1.0000000
# Ficedula_hypoleuca              0           0.4437522        0.6390457
# Emberiza_sulphurata             0           0.4437522        0.5328950
# Parus_major                     0           0.4437522        0.5288567
# Parus_caeruleus                 0           0.4437522        0.5288567
#                     Ficedula_hypoleuca Emberiza_sulphurata Parus_major
# Gallus_gallus                0.0000000           0.0000000   0.0000000
# Cyanocitta_cristata          0.4437522           0.4437522   0.4437522
# Sturnus_vulgaris             0.6390457           0.5328950   0.5288567
# Ficedula_hypoleuca           1.0000000           0.5328950   0.5288567
# Emberiza_sulphurata          0.5328950           1.0000000   0.5288567
# Parus_major                  0.5288567           0.5288567   1.0000000
# Parus_caeruleus              0.5288567           0.5288567   0.8176725
#                     Parus_caeruleus
# Gallus_gallus             0.0000000
# Cyanocitta_cristata       0.4437522
# Sturnus_vulgaris          0.5288567
# Ficedula_hypoleuca        0.5288567
# Emberiza_sulphurata       0.5288567
# Parus_major               0.8176725
# Parus_caeruleus           1.0000000

# calculate lnRR and lnRR_var
dat <- effect_lnRR(dat_pred)
dat$Obs_ID <- 1:nrow(dat)

# test whether phylogeny should be included in our analyses
model_phy <- NULL
phylo_cor <- NULL
class(model_phy) <- "list"
class(phylo_cor) <- "list"

V <- vcalc(lnRR_var, cluster = Cohort_ID, subgroup = Obs_ID, data = dat)

for (i in 1:50) {
  # Calculate the VCV matrix for each trees 
  vcv_matrix <- vcv(trees[[i]], cor = T) 
  # Store the results in the list
  phylo_cor[[i]] <- vcv_matrix

# run meta-analysis
phylo_model <- rma.mv(yi = lnRR,
                          V = V,
                          random = list(~1 | Study_ID,
                                        ~1 | Shared_control_ID,
                                        ~1 | Obs_ID,
                                        ~1 | Bird_species,
                                        ~1 | Phylogeny
                                        ),
                          R = list(Phylogeny = phylo_cor[[i]]), 
                          test = "t",
                          method = "REML",
                          sparse = TRUE,
                          data = dat)

 model_phy[[i]] <- phylo_model
 print(i)

}

summary(model_phy[[1]]) # example of summary
# Multivariate Meta-Analysis Model (k = 117; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -155.4106   310.8212   322.8212   339.3428   323.5919   

# Variance Components:

#            estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.1725  0.4153     29     no  Shared_control_ID   no 
# sigma^2.3  0.5565  0.7460    117     no             Obs_ID   no 
# sigma^2.4  0.0000  0.0000      7     no       Bird_species   no 
# sigma^2.5  0.0000  0.0000      7     no          Phylogeny  yes 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0802  0.1186  0.6762  116  0.5003  -0.1547  0.3150     
