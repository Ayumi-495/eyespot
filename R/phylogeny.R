# read libraries
library(ape)
library(here)
library(metafor)
library(orchaRd)
library(phangorn)
library(tidyverse)

effect_lnRR <- function(dt) {
  
  # replace 0 in "C_mean", "T_sd", "C_sd", "C_proportion" with each minimum values
  # if proportion too extreme value, replace minimum value (only "T_proportion")
  
  dt <- dt %>% 
    mutate(C_mean = ifelse(C_mean == 0, 0.04, C_mean))
  dt <- dt %>% 
    mutate(T_sd = ifelse(T_sd == 0, 0.01, T_sd))
  dt <- dt %>% 
    mutate(C_sd = ifelse(C_sd == 0, 0.05, C_sd))
  dt <- dt %>% 
    mutate(C_proportion = ifelse(C_proportion == 0, 0.01, C_proportion))
  dt <- dt %>% 
    mutate(T_proportion = ifelse(T_proportion < 0.01, 0.01, T_proportion))
  
  # copy dataset for adding effect size and its variation (lnRR /lnRR_var) column
  dt1 <- dt %>% 
    mutate(lnRR     = NA,
           lnRR_var = NA)
  
  for(i in seq_len(nrow(dt1))) {
    
    Tn <- dt1$Tn[i]
    Cn <- dt1$Cn[i]
    T_mean <- dt1$T_mean[i]
    C_mean <- dt1$C_mean[i]
    T_proportion <- dt1$T_proportion[i]
    C_proportion <- dt1$C_proportion[i]
    T_sd <- dt1$T_sd[i]
    C_sd <- dt1$C_sd[i]
    Response <- dt1$Response[i]
    Measure <- dt1$Measure[i]
    Study_design <- dt1$Study_design[i]
    Direction <- dt1$Direction[i]
    
    # continuous data - using escalc() (metafor package)
    if (Response == "continuous" & Study_design == "independent") {
      
      # reverse means - response expect decrease when birds present eyespots
         
      effect <- escalc(measure = "ROM", 
                       n1i = Tn, n2i = Cn, 
                       m1i = T_mean, m2 = C_mean, 
                       sd1i = T_sd, sd2i = C_sd,
                       var.names = c("lnRR", "lnRR_var")
      )
      
      dt1$lnRR[i] <- effect$lnRR
      dt1$lnRR_var[i] <- effect$lnRR_var
      
    }
    else if (Response == "continuous" & Study_design == "dependent") {
      
      # reverse means - response expect decrease when birds present eyespots
        
      effect <- escalc(measure = "ROMC", 
                       ni = (Tn + Cn)/2, 
                       m1i = T_mean, m2 = C_mean, 
                       sd1i = T_sd, sd2i = C_sd,
                       ri = 0.5,
                       var.names = c("lnRR", "lnRR_var")
      )
      
      dt1$lnRR[i] <- effect$lnRR
      dt1$lnRR_var[i] <- effect$lnRR_var
      
    }
    
    # proportion data (no sd values)
    else if (Response == "proportion1" & Study_design == "independent") {
      
      T_proportion <- replace(T_proportion, Direction == "decrease",
                        (1-T_mean[Direction == "decrease"]))
      C_proportion <- replace(C_proportion, Direction == "decrease",
                        (1-C_mean[Direction == "decrease"]))
      
      # transform proportion value
      asin_trans <- function(proportion) {
        trans <- asin(sqrt(proportion)) 
        trans
      }
      
      T_proportion <- asin_trans(T_proportion)
      C_proportion <- asin_trans(C_proportion)
      
      # calculate lnRR and lnRR variance   
      lnRR_pro1 <- log(T_proportion / C_proportion)
      lnRR_var_pro1 <- (1 / sqrt(8))^2 * (1 / (T_proportion^2 * Tn) +
                                            1 / (C_proportion^2 * Cn))
      
      dt1$lnRR[i] <- lnRR_pro1
      dt1$lnRR_var[i] <- lnRR_var_pro1
      
    }
    else if (Response == "proportion1" & Study_design == "dependent") {
      
      T_proportion <- replace(T_proportion, Direction == "decrease",
                              (1-T_mean[Direction == "decrease"]))
      C_proportion <- replace(C_proportion, Direction == "decrease",
                              (1-C_mean[Direction == "decrease"]))
      
      # transform proportion value
      asin_trans <- function(proportion) {
        trans <- asin(sqrt(proportion)) 
        trans
      }
      
      T_proportion <- asin_trans(T_proportion)
      C_proportion <- asin_trans(C_proportion)
      
      # calculate lnRR and lnRR variance   
      lnRR_pro1 <- log(T_proportion / C_proportion)
      lnRR_var_pro1 <- (1 / sqrt(8))^2 * (1 / (T_proportion^2 * Tn) +
                                            1 / (C_proportion^2 * Cn)) - 
        2 * 0.5 * sqrt((1 / sqrt(8))^2 * (1 / (T_proportion^2 * Tn))) * 
        sqrt((1 / sqrt(8))^2 * (1 / (C_proportion^2 * Cn)))
      
      dt1$lnRR[i] <- lnRR_pro1
      dt1$lnRR_var[i] <- lnRR_var_pro1
      
    }
    # proportion data (exist sd values) 
    else if (Response == "proportion2" & Study_design == "independent") {
      
      # transform proportion mean value
      asin_trans <- function(proportion) {
        trans <- asin(sqrt(proportion)) 
        trans
      }
      
      T_SD <- T_sd^2/(4*(T_proportion)*(1-(T_proportion)))
      C_SD <- C_sd^2/(4*(C_proportion)*(1-(C_proportion)))
      
      T_proportion <- asin_trans(T_proportion)
      C_proportion <- asin_trans(C_proportion)
      
      # calculate lnRR and lnRR variance 
      lnRR_pro2 <- log(T_proportion / C_proportion)
      lnRR_var_pro2 <- (T_SD)^2 * (1 / (T_proportion^2 * Tn)) +
        (C_SD)^2 * (1 / (C_proportion^2 * Cn))
      
      dt1$lnRR[i] <- lnRR_pro2
      dt1$lnRR_var[i] <- lnRR_var_pro2
      
    }
    
    else if (Response == "proportion2" & Study_design == "dependent") {
      
      # transform proportion mean value
      asin_trans <- function(proportion) {
        trans <- asin(sqrt(proportion)) 
        trans
      }
      
      T_SD <- T_sd^2/(4*(T_proportion)*(1-(T_proportion)))
      C_SD <- C_sd^2/(4*(C_proportion)*(1-(C_proportion)))
      
      T_proportion <- asin_trans(T_proportion)
      C_proportion <- asin_trans(C_proportion)
      
      # calculate lnRR and lnRR variance 
      lnRR_pro2 <- log(T_proportion / C_proportion)
      lnRR_var_pro2 <- (T_SD)^2 * (1 / (T_proportion^2 * Tn)) +
        (C_SD)^2 * (1 / (C_proportion^2 * Cn)) -
        2 * 0.5 * sqrt((T_SD)^2 * (1 / (T_proportion^2 * Tn)) ) *
        sqrt((C_SD)^2 * (1 / (C_proportion^2 * Cn)))
      
      dt1$lnRR[i] <- lnRR_pro2
      dt1$lnRR_var[i] <- lnRR_var_pro2
      
    }
  }
  
  return(dt1)
  
}



# get data
dat_pred <- read_csv(here("data/predator_22072023.csv"))
dim(dat_pred)

dat2 <- effect_lnRR(dat_pred)
dat2$Obs_ID <- 1:nrow(dat2)

# phylogenetic tree
trees <- read.nexus(here("data/bird_phy.nex"))
t <- 50 # number of trees
phylo_vcv <- lapply(1:t, function(i) vcv(trees[[i]], corr = TRUE))
phylo_vcv

# test 1 - loop through all trees and find minimum AIC value tree
ma_pred_test <- lapply(phylo_vcv, function(phylo_vcv) {
  rma.mv(yi = lnRR,
         V = lnRR_var,
         random = list(~1 | Study_ID,
                       ~1 | Cohort_ID, 
                       ~1 | Shared_control_ID,
                       ~1 | Bird_species,
                       ~1 | Bird_species,
                       ~1 |Obs_ID),
         R = list(Bird_species = phylo_vcv), 
         test = "t",
         method = "REML",
         sparse = TRUE,
         data = dat2)
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
# sigma^2.5  0.0000  0.0000      7     no       Bird_species  yes 
#sigma^2.6  0.5344  0.7310    117     no             Obs_ID   no 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139    

i2_ml(ma_pred_test[[37]])
#        I2_Total             I2_Study_ID         I2_Cohort_ID 
#        9.960611e+01         5.535444e-07         1.631045e+01 
#    I2_Shared_control_ID    I2_Bird_species      I2_Bird_species 
#       1.120321e+01         2.046888e-09         2.046888e-09 
#           I2_Obs_ID 
#        7.209245e+01 
### 

# test 2 - average the phylogenetic covariance matrix and use it in the model
cov_comb <- Reduce("+", phylo_vcv) / length(phylo_vcv)

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

ma_pred_test_1 <- rma.mv(yi = lnRR,
                         V = lnRR_var,
                         random = list(~1 | Study_ID,
                                       ~1 | Cohort_ID, 
                                       ~1 | Shared_control_ID,
                                       ~1 | Bird_species,
                                       ~1 | Bird_species,
                                       ~1 |Obs_ID),
                         R = list(Bird_species = cov_comb), 
                         test = "t",
                         method = "REML",
                         sparse = TRUE,
                         data = dat2)

summary(ma_pred_test_1)

# Multivariate Meta-Analysis Model (k = 117; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -155.1741   310.3482   324.3482   343.6233   325.3852   

# Variance Components:
#             estim    sqrt  nlvls  fixed             factor    R 
# sigma^2.1  0.0000  0.0001     18     no           Study_ID   no 
# sigma^2.2  0.1209  0.3477     33     no          Cohort_ID   no 
# sigma^2.3  0.0830  0.2882     29     no  Shared_control_ID   no 
# sigma^2.4  0.0000  0.0000      7     no       Bird_species  yes 
# sigma^2.5  0.0000  0.0000      7     no       Bird_species  yes 
# sigma^2.6  0.5344  0.7310    117     no             Obs_ID   no 

# Test for Heterogeneity:
# Q(df = 116) = 5282.9638, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval    ci.lb   ci.ub    
#   0.0747  0.1208  0.6187  116  0.5373  -0.1645  0.3139    

i2_ml(ma_pred_test_1)
#             I2_Total          I2_Study_ID         I2_Cohort_ID 
#         9.960611e+01         5.619035e-07         1.631044e+01 
# I2_Shared_control_ID      I2_Bird_species      I2_Bird_species 
#         1.120322e+01         2.087399e-09         2.087399e-09 
#            I2_Obs_ID 
#         7.209245e+01 