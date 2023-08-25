# read libraries
if (!require(MetBrewer)) {
  install.packages("MetBrewer")
} # 
pacman::p_load(here, MetBrewer, orchaRd, metafor) #delete MetBrewer parts in code, but I will use this package
source("R/function_2.R")

# get data
dat_all <-  read_csv(here("data/all_15082023.csv"))
dim(dat_all)

# add log-transform diameter and area
# ASK - I think I should use log-transformed background area in our analyses rather than natural background area
dat_all$Log_diameter <- log(dat_all$Diameter_pattern)
dat_all$Log_area <- log(dat_all$Area_pattern)
dat_all$Log_background <- log(dat_all$Area_background)

# turn all character strings to factor
dat_all <- dat_all %>%
  mutate_if(is.character, as.factor)
summary(dat_all)

# calculate lnRR and lnRR variance
dat <- effect_lnRR(dat_all)
dat$Obs_ID <- 1:nrow(dat)
hist(dat$lnRR) 
hist(dat$lnRR_var)

# a. meta-analysis
##################
# overall effect #
##################
# I exclude cohort_ID because sigma^2.2 = 0 and I2 = 0

# use vcalc to calculate variance-covariance matrix
VCV <- vcalc(vi = lnRR_var, 
             cluster = Cohort_ID,
             obs = Obs_ID,
             rho = 0.5,
             data = dat)

# now results are a bit different 
ma_all <- rma.mv(yi = lnRR,
                  V = VCV, 
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat)

summary(ma_all)
#    logLik   Deviance        AIC        BIC       AICc   
# -256.8776   513.7551   521.7551   536.0285   521.9108   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0850  0.2916     32     no           Study_ID 
# sigma^2.2  0.0163  0.1279     88     no  Shared_control_ID 
# sigma^2.3  0.2371  0.4870    263     no             Obs_ID 

# Test for Heterogeneity:
# Q(df = 262) = 6425.6391, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval   ci.lb   ci.ub     
#   0.2298  0.0713  3.2225  262  0.0014  0.0894  0.3702  ** 

i2_ml(ma_all)
#              I2_Total          I2_Study_ID I2_Shared_control_ID 
#             98.472649            24.736691             4.755604 
#            I2_Obs_ID 
#             68.980353 

orchard_plot(ma_all,
              group = "Study_ID",
              xlab = "log response ratio (lnRR)", angle = 45) +
              scale_x_discrete(labels = c("Overall effect"))

caterpillars(ma_all, group = "Study_ID", xlab = "log response ratio (lnRR)")

# b. meta-regression
##############################
# 0. eyespot or conspicuous? #
##############################
## simple model remove intercept
# TODO use VCV

mr_eyespot <- rma.mv(yi = lnRR,
                    V = lnRR_var, 
                    mods = ~ Treatment_stimulus -1,
                    random = list(~1 | Study_ID,
                                  ~1 | Shared_control_ID,
                                  ~1 | Obs_ID),
                    test = "t",
                    method = "REML", 
                    sparse = TRUE,
                    data = dat)

summary(mr_eyespot)
# Multivariate Meta-Analysis Model (k = 263; method: REML)
#     logLik   Deviance        AIC        BIC       AICc      
# -255.8370   511.6740   521.6740   539.4966   521.9093   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0893  0.2989     32     no           Study_ID 
# sigma^2.2  0.0164  0.1282     88     no  Shared_control_ID 
# sigma^2.3  0.2372  0.4870    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6354.1948, p-val < .0001

# Test of Moderators (coefficients 1:2):
# F(df1 = 2, df2 = 261) = 5.2968, p-val = 0.0056

# Model Results:

#                                 estimate      se    tval   df    pval    ci.lb 
# Treatment_stimulus conspicuous    0.1709  0.1127  1.5168  261  0.1305  -0.0510 
# Treatment_stimulus eyespot        0.2545  0.0804  3.1647  261  0.0017   0.0961 
#                                  ci.ub     
# Treatment_stimulus conspicuous  0.3927     
# Treatment_stimulus eyespot      0.4128  ** 

mr_eyespot1 <- rma.mv(yi = lnRR,
                      V = lnRR_var, 
                      mods = ~ Treatment_stimulus,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = dat)

summary(mr_eyespot1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -255.8370   511.6740   521.6740   539.4966   521.9093   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0893  0.2989     32     no           Study_ID 
# sigma^2.2  0.0164  0.1282     88     no  Shared_control_ID 
# sigma^2.3  0.2372  0.4870    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6354.1948, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 0.4755, p-val = 0.4911

# Model Results:

#                             estimate      se    tval   df    pval    ci.lb 
# intrcpt                       0.1709  0.1127  1.5168  261  0.1305  -0.0510 
# Treatment_stimulus eyespot    0.0836  0.1212  0.6896  261  0.4911  -0.1551 
#                              ci.ub    
# intrcpt                     0.3927    
# Treatment_stimulus eyespot  0.3223    

orchard_plot(mr_eyespot,
            mod = "Treatment_stimulus",
            group = "Study_ID",
            xlab = "log response ratio (lnRR)",
            angle = 45)

##########################
# 1. Diameter of pattern #
##########################
## simple model
mr_diameter <- rma.mv(yi = lnRR,
                      V = lnRR_var, 
                      mods = ~ Diameter_pattern,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = dat)

summary(mr_diameter)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

# logLik   Deviance        AIC        BIC       AICc   
# -255.1009   510.2018   520.2018   538.0244   520.4371   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0909  0.3015     32     no           Study_ID 
# sigma^2.2  0.0144  0.1202     88     no  Shared_control_ID 
# sigma^2.3  0.2373  0.4871    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6161.0355, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 1.1046, p-val = 0.2942

# Model Results:

#                   estimate      se    tval   df    pval    ci.lb   ci.ub    
# intrcpt             0.1366  0.1154  1.1836  261  0.2377  -0.0907  0.3639    
# Diameter_pattern    0.0109  0.0104  1.0510  261  0.2942  -0.0095  0.0314

bubble_plot(mr_diameter,
            mod = "Diameter_pattern",
            group = "Study_ID",
            xlab = "Pattern diameter (mm)")

#CHECK 
## use log-transformed diameter
mr_diameter_log <- rma.mv(yi = lnRR,
                          V = lnRR_var, 
                          mods = ~ Log_diameter,
                          random = list(~1 | Study_ID,
                                        ~1 | Shared_control_ID,
                                        ~1 | Obs_ID),
                          test = "t",
                          method = "REML", 
                          sparse = TRUE,
                          data = dat)

summary(mr_diameter_log)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -253.2494   506.4989   516.4989   534.3215   516.7342   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0824  0.2871     32     no           Study_ID 
# sigma^2.2  0.0139  0.1178     88     no  Shared_control_ID 
# sigma^2.3  0.2330  0.4827    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6268.5425, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 5.3256, p-val = 0.0218

# Model Results:

#               estimate      se     tval   df    pval    ci.lb   ci.ub    
# intrcpt        -0.1742  0.1886  -0.9238  261  0.3564  -0.5456  0.1971    
# Log_diameter    0.2100  0.0910   2.3077  261  0.0218   0.0308  0.3893  *

bubble_plot(mr_diameter_log,
            mod = "Log_diameter",
            group = "Study_ID",
            xlab = "Log-transformed pattern diameter")

## add treatment stimulus type (eyespots or conspicuous)
mr_diameter_log1 <- rma.mv(yi = lnRR,
                        V = lnRR_var, 
                        mods = ~ Log_diameter + Treatment_stimulus,
                        random = list(~1 | Study_ID,
                                      ~1 | Shared_control_ID,
                                      ~1 | Obs_ID),
                        test = "t",
                        method = "REML",
                        sparse = TRUE,
                        data = dat)

summary(mr_diameter_log1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

# logLik   Deviance        AIC        BIC       AICc   
# -251.9584   503.9168   515.9168   537.2809   516.2488   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0894  0.2990     32     no           Study_ID 
# sigma^2.2  0.0145  0.1204     88     no  Shared_control_ID 
# sigma^2.3  0.2314  0.4810    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 6265.8920, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 3.1383, p-val = 0.0450

# Model Results:

#                             estimate      se     tval   df    pval    ci.lb 
# intrcpt                      -0.2873  0.2208  -1.3013  260  0.1943  -0.7222 
# Log_diameter                  0.2245  0.0933   2.4065  260  0.0168   0.0408 
# Treatment_stimulus eyespot    0.1212  0.1210   1.0018  260  0.3174  -0.1171 
#                              ci.ub    
# intrcpt                     0.1475    
# Log_diameter                0.4082  * 
# Treatment_stimulus eyespot  0.3595

## add log-transformed background area
mr_diameter_log2 <- rma.mv(yi = lnRR,
                          V = lnRR_var, 
                          mods = ~ Log_diameter + Log_background,
                          random = list(~1 | Study_ID,
                                        ~1 | Shared_control_ID,
                                        ~1 | Obs_ID),
                          test = "t",
                          method = "REML",
                          sparse = TRUE,
                          data = dat)

summary(mr_diameter_log2)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

# logLik   Deviance        AIC        BIC       AICc   
# -247.8091   495.6182   507.6182   528.9823   507.9502   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0498  0.2232     32     no           Study_ID 
# sigma^2.2  0.0163  0.1275     88     no  Shared_control_ID 
# sigma^2.3  0.2300  0.4796    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 5989.5561, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 7.8896, p-val = 0.0005

# Model Results:

#                 estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt           1.1423  0.4439   2.5733  260  0.0106   0.2682   2.0165    * 
# Log_diameter      0.4218  0.1087   3.8816  260  0.0001   0.2078   0.6358  *** 
# Log_background   -0.2486  0.0784  -3.1714  260  0.0017  -0.4029  -0.0942   **

## add treatment stimulus type and log-transformed background area
mr_diameter_log3 <- rma.mv(yi = lnRR,
                          V = lnRR_var, 
                          mods = ~ Log_diameter + Treatment_stimulus + 
                                   Log_background,
                          random = list(~1 | Study_ID,
                                        ~1 | Shared_control_ID,
                                        ~1 | Obs_ID),
                          test = "t",
                          method = "REML",
                          sparse = TRUE,
                          data = dat)

summary(mr_diameter_log3)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -246.1965   492.3930   506.3930   531.2908   506.8392   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0581  0.2410     32     no           Study_ID 
# sigma^2.2  0.0166  0.1290     88     no  Shared_control_ID 
# sigma^2.3  0.2270  0.4765    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 259) = 5987.3348, p-val < .0001

# Test of Moderators (coefficients 2:4):
# F(df1 = 3, df2 = 259) = 5.6684, p-val = 0.0009

# Model Results:

#                             estimate      se     tval   df    pval    ci.lb 
# intrcpt                       1.1004  0.4645   2.3692  259  0.0186   0.1858 
# Log_diameter                  0.4466  0.1123   3.9774  259  <.0001   0.2255 
# Treatment_stimulus eyespot    0.1537  0.1137   1.3520  259  0.1776  -0.0701 
# Log_background               -0.2645  0.0817  -3.2373  259  0.0014  -0.4253 
#                               ci.ub      
# intrcpt                      2.0149    * 
# Log_diameter                 0.6676  *** 
# Treatment_stimulus eyespot   0.3775      
# Log_background              -0.1036   **

######################
# 2. area of pattern #
######################
## simple model
mr_area <- rma.mv(yi = lnRR,
                  V = lnRR_var,
                  mods = ~ Area_pattern,
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat)

summary(mr_area)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -255.3630   510.7260   520.7260   538.5486   520.9613   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0942  0.3069     32     no           Study_ID 
# sigma^2.2  0.0148  0.1215     88     no  Shared_control_ID 
# sigma^2.3  0.2378  0.4877    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6098.4231, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 0.2576, p-val = 0.6122

# Model Results:

#               estimate      se    tval   df    pval    ci.lb   ci.ub    
# intrcpt         0.2089  0.0858  2.4355  261  0.0155   0.0400  0.3777  * 
# Area_pattern    0.0003  0.0005  0.5076  261  0.6122  -0.0008  0.0013


bubble_plot(mr_area,
            mod = "Area_pattern",
            group = "Study_ID",
            xlab = "Pattern area (mm2)")

# CHECK - use log-transformed pattern area
mr_area_log <- rma.mv(yi = lnRR,
                      V = lnRR_var,
                      mods = ~ Log_area,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML",
                      sparse = TRUE,
                      data = dat)

summary(mr_area_log)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -252.5353   505.0706   515.0706   532.8932   515.3059   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0863  0.2938     32     no           Study_ID 
# sigma^2.2  0.0124  0.1112     88     no  Shared_control_ID 
# sigma^2.3  0.2315  0.4811    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6258.2214, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 6.7146, p-val = 0.0101

# Model Results:

#           estimate      se     tval   df    pval    ci.lb   ci.ub    
# intrcpt    -0.2078  0.1835  -1.1329  261  0.2583  -0.5691  0.1534    
# Log_area    0.1238  0.0478   2.5913  261  0.0101   0.0297  0.2178  *

bubble_plot(mr_area_log,
            mod = "Log_area",
            group = "Study_ID",
            xlab = "Log-transformed pattern area")

## add treatment stimulus type (eyespots or conspicuous)
mr_area_log1 <- rma.mv(yi = lnRR,
                      V = lnRR_var, 
                      mods = ~ Log_area + Treatment_stimulus,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = dat)

summary(mr_area_log1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#   logLik   Deviance        AIC        BIC       AICc   
# -251.5190   503.0379   515.0379   536.4020   515.3699   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0908  0.3014     32     no           Study_ID 
# sigma^2.2  0.0132  0.1147     88     no  Shared_control_ID 
# sigma^2.3  0.2310  0.4806    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 6257.1020, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 3.5728, p-val = 0.0295

# Model Results:

#                             estimate      se     tval   df    pval    ci.lb 
# intrcpt                      -0.2695  0.2042  -1.3197  260  0.1881  -0.6715 
# Log_area                      0.1248  0.0484   2.5803  260  0.0104   0.0296 
# Treatment_stimulus eyespot    0.0823  0.1200   0.6860  260  0.4933  -0.1540 
#                              ci.ub    
# intrcpt                     0.1326    
# Log_area                    0.2200  * 
# Treatment_stimulus eyespot  0.3186

## add log-transformed background area
mr_area_log2 <- rma.mv(yi = lnRR,
                      V = lnRR_var, 
                      mods = ~ Log_area + Log_background,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = dat)

summary(mr_area_log2)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -245.6071   491.2141   503.2141   524.5782   503.5462   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0511  0.2261     32     no           Study_ID 
# sigma^2.2  0.0132  0.1151     88     no  Shared_control_ID 
# sigma^2.3  0.2273  0.4767    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 5998.2405, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 10.1596, p-val < .0001

# Model Results:

#                 estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt           1.3700  0.4556   3.0068  260  0.0029   0.4728   2.2672   ** 
# Log_area          0.2581  0.0583   4.4296  260  <.0001   0.1434   0.3728  *** 
# Log_background   -0.2960  0.0813  -3.6396  260  0.0003  -0.4561  -0.1358  ***

mr_area_log3 <- rma.mv(yi = lnRR,
                      V = lnRR_var, 
                      mods = ~ Log_area + Treatment_stimulus + 
                               Log_background,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = dat)

summary(mr_area_log3)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#   logLik   Deviance        AIC        BIC       AICc   
# -244.6386   489.2771   503.2771   528.1749   503.7233   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0557  0.2359     32     no           Study_ID 
# sigma^2.2  0.0142  0.1190     88     no  Shared_control_ID 
# sigma^2.3  0.2264  0.4758    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 259) = 5997.5157, p-val < .0001

# Test of Moderators (coefficients 2:4):
# F(df1 = 3, df2 = 259) = 6.8174, p-val = 0.0002

# Model Results:

#                             estimate      se     tval   df    pval    ci.lb 
# intrcpt                       1.3350  0.4714   2.8320  259  0.0050   0.4067 
# Log_area                      0.2594  0.0591   4.3903  259  <.0001   0.1430 
# Treatment_stimulus eyespot    0.0822  0.1107   0.7426  259  0.4584  -0.1358 
# Log_background               -0.2996  0.0830  -3.6088  259  0.0004  -0.4632 
#                               ci.ub      
# intrcpt                      2.2633   ** 
# Log_area                     0.3757  *** 
# Treatment_stimulus eyespot   0.3002      
# Log_background              -0.1361  ***

########################
# 3. Number of pattern #
########################
## simple model
mr_num <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Number_pattern,
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat)

summary(mr_num)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#   logLik   Deviance        AIC        BIC       AICc   
# -253.9357   507.8713   517.8713   535.6939   518.1066   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0827  0.2877     32     no           Study_ID 
# sigma^2.2  0.0122  0.1106     88     no  Shared_control_ID 
# sigma^2.3  0.2340  0.4838    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6410.9284, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 5.0996, p-val = 0.0248

# Model Results:

#                 estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt           0.3884  0.0991   3.9201  261  0.0001   0.1933   0.5835  *** 
# Number_pattern   -0.0651  0.0288  -2.2582  261  0.0248  -0.1219  -0.0083    *

bubble_plot(mr_all3,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number of patterns")

## add treatment stimulus type (eyespots or conspicuous)
mr_num1 <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Number_pattern + Treatment_stimulus,
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML",
                  sparse = TRUE,
                  data = dat)

summary(mr_num1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -252.7886   505.5772   517.5772   538.9413   517.9093   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0863  0.2937     32     no           Study_ID 
# sigma^2.2  0.0126  0.1124     88     no  Shared_control_ID 
# sigma^2.3  0.2336  0.4834    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 6339.7133, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 2.8957, p-val = 0.0570

# Model Results:

#                             estimate      se     tval   df    pval    ci.lb 
# intrcpt                       0.3210  0.1283   2.5025  260  0.0129   0.0684 
# Number_pattern               -0.0668  0.0290  -2.3046  260  0.0220  -0.1239 
# Treatment_stimulus eyespot    0.1014  0.1196   0.8478  260  0.3973  -0.1341 
#                               ci.ub    
# intrcpt                      0.5736  * 
# Number_pattern              -0.0097  * 
# Treatment_stimulus eyespot   0.3369

## add log-transformed background area
mr_num2 <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Number_pattern + Log_background,
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML",
                  sparse = TRUE,
                  data = dat)

summary(mr_num2)
# logLik   Deviance        AIC        BIC       AICc   
# -252.2622   504.5243   516.5243   537.8884   516.8563   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0887  0.2978     32     no           Study_ID 
# sigma^2.2  0.0137  0.1169     88     no  Shared_control_ID 
# sigma^2.3  0.2324  0.4821    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 6072.0986, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 2.9741, p-val = 0.0528

# Model Results:

#                 estimate      se     tval   df    pval    ci.lb    ci.ub    
# intrcpt           0.8576  0.5050   1.6981  260  0.0907  -0.1369   1.8521  . 
# Number_pattern   -0.0667  0.0290  -2.2987  260  0.0223  -0.1238  -0.0096  * 
# Log_background   -0.0665  0.0703  -0.9461  260  0.3450  -0.2050   0.0719

## add treatment stimulus type and background area
mr_num3 <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Number_pattern + Log_background + Treatment_stimulus,
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML",
                  sparse = TRUE,
                  data = dat)

summary(mr_num3)

#################
# 4. background #
#################
mr_background <- rma.mv(yi = lnRR,
                        V = lnRR_var, 
                        mods = ~ Area_background,
                        random = list(~1 | Study_ID,
                                      ~1 | Shared_control_ID,
                                      ~1 | Obs_ID),
                        test = "t",
                        method = "REML",
                        sparse = TRUE,
                        data = dat)

summary(mr_background)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -251.9129   503.8258   513.8258   531.6484   514.0611   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0800  0.2828     32     no           Study_ID 
# sigma^2.2  0.0131  0.1144     88     no  Shared_control_ID 
# sigma^2.3  0.2330  0.4827    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6073.6810, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 7.6165, p-val = 0.0062

# Model Results:

#                  estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt            0.3287  0.0781   4.2062  261  <.0001   0.1748   0.4825  *** 
# Area_background   -0.0000  0.0000  -2.7598  261  0.0062  -0.0001  -0.0000   **

bubble_plot(mr_background,
            mod = "Area_background",
            group = "Study_ID",
            xlab = "Background area (mm²)")

mr_background_log <- rma.mv(yi = lnRR,
                        V = lnRR_var, 
                        mods = ~ Log_background,
                        random = list(~1 | Study_ID,
                                      ~1 | Shared_control_ID,
                                      ~1 | Obs_ID),
                        test = "t",
                        method = "REML",
                        sparse = TRUE,
                        data = dat)

summary(mr_background_log)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

# logLik   Deviance        AIC        BIC       AICc   
# -255.2981   510.5961   520.5961   538.4187   520.8314   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0912  0.3020     32     no           Study_ID 
# sigma^2.2  0.0172  0.1313     88     no  Shared_control_ID 
# sigma^2.3  0.2362  0.4860    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6114.0099, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 0.6646, p-val = 0.4157

# Model Results:

#                 estimate      se     tval   df    pval    ci.lb   ci.ub    
# intrcpt           0.6367  0.5036   1.2641  261  0.2073  -0.3550  1.6284    
# Log_background   -0.0582  0.0714  -0.8152  261  0.4157  -0.1988  0.0824

bubble_plot(mr_background_log,
            mod = "Log_background",
            group = "Study_ID",
            xlab = "Log_background")


###########################################
# 5. type of prey (Artificial or natural) #
###########################################
## simple model remove intercept
mr_prey_type <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Type_prey -1,
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat)

summary(mr_prey_type)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -255.8353   511.6706   521.6706   539.4933   521.9059   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0938  0.3063     32     no           Study_ID 
# sigma^2.2  0.0161  0.1267     88     no  Shared_control_ID 
# sigma^2.3  0.2369  0.4867    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6337.4755, p-val < .0001

# Test of Moderators (coefficients 1:2):
# F(df1 = 2, df2 = 261) = 5.0461, p-val = 0.0071

# Model Results:

#                      estimate      se    tval   df    pval   ci.lb   ci.ub    
# Type_preyartificial    0.2057  0.0918  2.2401  261  0.0259  0.0249  0.3866  * 
# Type_preyreal          0.2761  0.1226  2.2525  261  0.0251  0.0347  0.5174  *

orchard_plot(mr_prey_type,
              mod = "Type_prey",
              group = "Study_ID",
              xlab = "Type of prey",
              angle = 45)

mr_prey_type1 <- rma.mv(yi = lnRR,
                        V = lnRR_var, 
                        mods = ~ Type_prey,
                        random = list(~1 | Study_ID,
                                      ~1 | Shared_control_ID,
                                      ~1 | Obs_ID),
                        test = "t",
                        method = "REML", 
                        sparse = TRUE,
                        data = dat)

summary(mr_prey_type1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -255.8353   511.6706   521.6706   539.4933   521.9059   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0938  0.3063     32     no           Study_ID 
# sigma^2.2  0.0161  0.1267     88     no  Shared_control_ID 
# sigma^2.3  0.2369  0.4867    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6337.4755, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 0.2109, p-val = 0.6465

# Model Results:

#                estimate      se    tval   df    pval    ci.lb   ci.ub    
# intrcpt          0.2057  0.0918  2.2401  261  0.0259   0.0249  0.3866  * 
# Type_preyreal    0.0703  0.1531  0.4592  261  0.6465  -0.2312  0.3719

####################
# 6. shape of prey #
####################
## simple model remove intercept
mr_prey_shape <- rma.mv(yi = lnRR,
                        V = lnRR_var, 
                        mods = ~ Shape_prey - 1,
                        random = list(~1 | Study_ID,
                                      ~1 | Shared_control_ID,
                                      ~1 | Obs_ID),
                                      test = "t",
                                      method = "REML", 
                                      sparse = TRUE,
                                      data = dat)

summary(mr_prey_shape)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -252.4490   504.8980   518.8980   543.7958   519.3442   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.1039  0.3223     32     no           Study_ID 
# sigma^2.2  0.0206  0.1436     88     no  Shared_control_ID 
# sigma^2.3  0.2323  0.4820    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 259) = 5984.8859, p-val < .0001

# Test of Moderators (coefficients 1:4):
# F(df1 = 4, df2 = 259) = 2.7782, p-val = 0.0274

# Model Results:

#                               estimate      se    tval   df    pval    ci.lb 
# Shape_preyabstract_butterfly    0.3274  0.1386  2.3620  259  0.0189   0.0545 
# Shape_preyabstract_stimuli      0.0138  0.2339  0.0591  259  0.9529  -0.4467 
# Shape_preybutterfly             0.2794  0.1261  2.2157  259  0.0276   0.0311 
# Shape_preycaterpillar           0.1265  0.1606  0.7878  259  0.4315  -0.1897 
#                                ci.ub    
# Shape_preyabstract_butterfly  0.6003  * 
# Shape_preyabstract_stimuli    0.4743    
# Shape_preybutterfly           0.5277  * 
# Shape_preycaterpillar         0.4428

orchard_plot(mr_prey_shape,
              mod = "Shape_prey",
              group = "Study_ID",
              xlab = "Shape of prey", angle = 45)

##############################
# 7. all potential moderator # 
##############################
#  mods = ~ Treatment_stimulus + Log_diameter + Log_background +
#            Log_area + Number_pattern + Type_prey + Shape_prey,
# Warning message: Redundant predictors dropped from the model.

mr_all <- rma.mv(yi = lnRR,
                V = lnRR_var, 
                mods = ~ Treatment_stimulus + Log_diameter + Log_background +
                          Log_area + Number_pattern + Type_prey + Shape_prey,
                random = list(~1 | Study_ID,
                              ~1 | Shared_control_ID,
                              ~1 | Obs_ID),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dat)

summary(mr_all)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -238.7127   477.4254   501.4254   543.8734   502.7200   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0532  0.2306     32     no           Study_ID 
# sigma^2.2  0.0202  0.1420     88     no  Shared_control_ID 
# sigma^2.3  0.2208  0.4699    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 254) = 5622.9922, p-val < .0001

# Test of Moderators (coefficients 2:9):
# F(df1 = 8, df2 = 254) = 3.2167, p-val = 0.0017

# Model Results:

#                              estimate      se     tval   df    pval    ci.lb 
# intrcpt                        0.5475  0.6470   0.8462  254  0.3982  -0.7267 
# Treatment_stimulus eyespot     0.1418  0.1464   0.9685  254  0.3337  -0.1465 
# Log_diameter                  -0.1597  0.3942  -0.4052  254  0.6857  -0.9360 
# Log_background                -0.1807  0.1049  -1.7234  254  0.0860  -0.3872 
# Log_area                      0.3628  0.2150   1.6873  254  0.0928  -0.0606 
# Number_pattern                -0.0189  0.0326  -0.5787  254  0.5633  -0.0831 
# Type_preyreal                 -0.0485  0.1787  -0.2715  254  0.7862  -0.4005 
# Shape_prey abstract_stimuli   -0.6749  0.3391  -1.9903  254  0.0476  -1.3427 
# Shape_prey caterpillar        -0.0317  0.1908  -0.1660  254  0.8683  -0.4073 
#                                ci.ub    
# intrcpt                       1.8216    
# Treatment_stimulus eyespot    0.4300    
# Log_diameter                  0.6166    
# Log_background                0.0258  . 
# Log_area                      0.7862  . 
# Number_pattern                0.0454    
# Type_preyreal                 0.3034    
# Shape_prey abstract_stimuli  -0.0071  * 
# Shape_prey caterpillar        0.3440

# ASK now doing
####################
# publication bias #
####################
# traditional funnel plot
funnel(ma_all)

df_bias <- dat %>% mutate(sqrt_inv_e_n = sqrt((Cn + Tn)/(Cn * Tn)))
data.frame(df_bias)
bias_model <- rma.mv(yi = lnRR,
                      V = lnRR_var, 
                      mods = ~sqrt_inv_e_n,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = df_bias)

summary(bias_model)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

# logLik   Deviance        AIC        BIC       AICc   
# -255.7498   511.4995   521.4995   539.3221   521.7348   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0751  0.2741     32     no           Study_ID 
# sigma^2.2  0.0215  0.1467     88     no  Shared_control_ID 
# sigma^2.3  0.2369  0.4857    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6081.2099, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 1, df2 = 261) = 0.5047, p-val = 0.4781

# Model Results:

#               estimate      se     tval   df    pval    ci.lb   ci.ub    
# intrcpt         0.3137  0.1401   2.2400  261  0.0259   0.0379  0.5895  * 
# sqrt_inv_e_n   -0.2866  0.4034  -0.7104  261  0.4781  -1.0810  0.5078    

bubble_plot(bias_model,
            mod = "sqrt_inv_e_n",
            group = "Study_ID",
            xlab = "sqrt_inv_e_n")

year_model <- rma.mv(yi = lnRR,
                      V = lnRR_var, 
                      mods = ~ Year,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = df_bias)

summary(year_model)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -255.8732   511.7464   521.7464   539.5690   521.9817   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0877  0.2961     32     no           Study_ID 
# sigma^2.2  0.0167  0.1291     88     no  Shared_control_ID 
# sigma^2.3  0.2377  0.4876    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6411.1908, p-val < .0001

# Test of Moderators (coefficient 2):
#F(df1 = 1, df2 = 261) = 0.0529, p-val = 0.8182

#Model Results:

#          estimate       se     tval   df    pval     ci.lb    ci.ub    
# intrcpt    3.5299  14.3444   0.2461  261  0.8058  -24.7156  31.7753    
# Year      -0.0016   0.0071  -0.2300  261  0.8182   -0.0157   0.0124    

bubble_plot(year_model,
            mod = "Year",
            group = "Study_ID",
            xlab = "Year")

###########
# dataset #
###########
mr_dataset <- rma.mv(yi = lnRR,
                      V = lnRR_var,
                      mods = ~ Dataset,
                      random = list(~1 | Study_ID,
                                    ~1 | Shared_control_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML",
                      sparse = TRUE,
                      data = dat)
summary(mr_dataset)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -255.4718   510.9436   520.9436   538.7662   521.1789   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0901  0.3002     32     no           Study_ID 
# sigma^2.2  0.0169  0.1300     88     no  Shared_control_ID 
# sigma^2.3  0.2362  0.4860    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6189.5292, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 0.7658, p-val = 0.3823

# Model Results:

#              estimate      se    tval   df    pval    ci.lb   ci.ub    
# intrcpt        0.1690  0.1009  1.6744  261  0.0952  -0.0297  0.3678  . 
# Dataset prey    0.1213  0.1387  0.8751  261  0.3823  -0.1517  0.3944

orchard_plot(mr_dataset,
              mod = "Dataset",
              group = "Study_ID",
              xlab = "Dataset",
              angle = 45)