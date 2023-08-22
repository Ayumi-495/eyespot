# test whether phylogeny should be included in our analyses

# read libraries
pacman::p_load(ape, here, tidyverse, miWQS, miceadds)
source("R/function_2.R")

# data and phylogeny
dat_predator <- read_csv(here("data/predator_22072023.csv"))
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
dat_pred <- effect_lnRR(dat_predator)
dat_pred$Obs_ID <- 1:nrow(dat_pred)

# check heterogeneity of phylogeny
model_phy <- NULL
phylo_cor <- NULL
class(model_phy) <- "list"
class(phylo_cor) <- "list"

V <- vcalc(lnRR_var, cluster = Cohort_ID, subgroup = Obs_ID, data = dat_pred)

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
                                        ~1 | Cohort_ID,
                                        ~1 | Obs_ID,
                                        ~1 | Bird_species,
                                        ~1 | Phylogeny
                                        ),
                          R = list(Phylogeny = phylo_cor[[i]]), 
                          test = "t",
                          method = "REML",
                          sparse = TRUE,
                          data = dat_pred)

  model_phy[[i]] <- phylo_model
  print(i)

}

summary(model_phy[[1]]) # example of summary
# Multivariate Meta-Analysis Model (k = 117; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -154.6149   309.2298   323.2298   342.5049   324.2668    

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.0923  0.3037     29     no  Shared_control_ID   no 
# sigma^2.3  0.1009  0.3177     33     no          Cohort_ID   no 
# sigma^2.4  0.5323  0.7296    117     no             Obs_ID   no 
# sigma^2.5  0.0000  0.0000      7     no       Bird_species   no 
# sigma^2.6  0.0000  0.0000      7     no          Phylogeny  yes 

# Test for Heterogeneity:
# Q(df = 116) = 5316.7677, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.1394  0.1192  1.1690  116  0.2448  -0.0968  0.3755       

# ASK - extract sigma^2　- I am not sure whether this is correct ＊ 間違っている and いらないかも
sigma2_mod <- do.call(rbind, lapply(model_phy, function(x) x$sigma2))
sigma2_mod <- data.frame(sigma2_mod)
colnames(sigma2_mod) <- c("sigma^2.1_Study_ID", "sigma^2.2_SharedControl_ID", "sigma^2.3_Cohort_ID",
                          "sigma^2.4_Obs_ID", "sigma^2.5_BirdSpecies", "sigma^2.6_Phylogeny")

print(sigma2_mod)

colMeans(sigma2_mod)
#         sigma^2.1_Study_ID sigma^2.2_SharedControl_ID 
#               6.604261e-10               2.888791e-10  
#        sigma^2.3_Cohort_ID           sigma^2.4_Obs_ID 
#               1.009451e-01               5.323456e-01 
#      sigma^2.5_BirdSpecies        sigma^2.6_Phylogeny 
#               6.604261e-10               2.888791e-10

i2_ml(model_phy[[1]]) # example of i2_ml
#             I2_Total          I2_Study_ID I2_Shared_control_ID 
#         9.959916e+01         4.743726e-07         1.266457e+01 
#         I2_Cohort_ID            I2_Obs_ID      I2_Bird_species 
#         1.385718e+01         7.307740e+01         8.611110e-08 
#         I2_Phylogeny 
#         4.416624e-08


# BUG
# When I want to use the pool.mi(), I need SD values of each model.  

# Extract the mean and standard deviation parameters from the rma.mv results
# to.pool <- array(NA, dim = c(2, 2, 50))
# for (i in 1:50) {
#  to.pool[1, 1, i] <- model_phy[[i]]$b
#  to.pool[1, 2, i] <- sqrt(model_phy[[i]]$V) - Do I need to extract SD here?
# }

# Then, I tried to use pool() function (mice package) and pool_mi() function (miceadds package)
# but I got an error message in pool_mi().

# use pool function
summary(pool(model_phy))
#      term  estimate std.error statistic      df   p.value
# 1 overall 0.1393776 0.1192319  1.168962 114.039 0.2448566

# use pool_mi function
# extract estimates and standard errors from each model
estimates <- lapply(model_phy, function(x) x$b)
se <- lapply(model_phy, function(x) x$se)


# use pool_mi function
results <- pool_mi(qhat = estimates, se = se)
# Error in u[ii, , ] <- u0[[ii]] : replacement has length zero