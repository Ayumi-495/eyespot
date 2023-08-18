# read libraries
pacman::p_load(here, MetBrewer, orchaRd)
source("R/function_2.R")

# get data
dat_all <-  read_csv(here("data/all_15082023.csv"))
dim(dat_all)


# turn all character strings to factor
dat_all <- dat_all %>%
  mutate_if(is.character, as.factor)

summary(dat_all)

# calculate lnRR and lnRR variance
dat <- effect_lnRR(dat_all)
dat$Obs_ID <- 1:nrow(dat)

hist(dat$lnRR) 
hist(dat$lnRR_var)


# meta-analysis 
# I exclude cohort_ID because sigma^2.2 = 0 and I2 = 0

ma_all <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat)

summary(ma_all)

#    logLik   Deviance        AIC        BIC       AICc   
# -259.7358   519.4716   527.4716   541.7450   527.6273   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0785  0.2802     32     no           Study_ID 
# sigma^2.2  0.0235  0.1534     88     no  Shared_control_ID 
# sigma^2.3  0.2429  0.4928    263     no             Obs_ID 

# Test for Heterogeneity:
# Q(df = 262) = 6465.9171, p-val < .0001

# Model Results:
# estimate      se    tval   df    pval   ci.lb   ci.ub     
#   0.2056  0.0707  2.9078  262  0.0040  0.0664  0.3448  ** 

i2_all <- i2_ml(ma_all)
i2_all
#             I2_Total          I2_Study_ID         I2_Shared_control_ID 
#         98.502998            22.418191             6.722398 
#            I2_Obs_ID 
#         69.362409 

p1_all <-  orchard_plot(ma_all,
                        group = "Study_ID",
                        xlab = "log response ratio (lnRR)", angle = 45) +
                        scale_x_discrete(labels = c("Overall effect")) +
                        scale_fill_manual(values = met.brewer("Homer2")) +
                        scale_colour_manual(values = met.brewer("Homer2"))

p1_all
ggsave("overall_all_11Aug.pdf", dpi = 450)

p1_all_cat <- caterpillars(ma_all, group = "Study_ID", xlab = "log response ratio (lnRR)")

p1_all_cat
ggsave("overall_cat_all.pdf", dpi = 450)

# meta-regression 

##############################
# 1. eyespot or conspicuous? #
##############################
mr_all <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Treatment_stimulus,
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat)

summary(mr_all)
#     logLik   Deviance        AIC        BIC       AICc      
#  -258.6788   517.3576   527.3576   545.1802   527.5929     

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0831  0.2882     32     no           Study_ID  
# sigma^2.2  0.0234  0.1529     88     no  Shared_control_ID  
# sigma^2.3  0.2430  0.4929    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6388.1567, p-val < .0001

# Test of Moderators (coefficients 1:2):
# F(df1 = 2, df2 = 261) = 4.3523, p-val = 0.0138

# Model Results:

#                                 estimate      se    tval   df    pval    ci.lb 
# Treatment_stimulus conspicuous    0.1449  0.1125  1.2881  261  0.1989  -0.0766 
# Treatment_stimulus eyespot        0.2308  0.0801  2.8819  261  0.0043   0.0731
#                                   ci.ub     
# Treatment_stimulus conspicuous   0.3664   
# Treatment_stimulus eyespot       0.3885 ** 

r2_ml(mr_all)
#   R2_marginal R2_conditional 
#   0.004373912    0.307613500   

p2_all <- orchard_plot(mr_all,
                       mod = "Treatment_stimulus",
                       group = "Study_ID",
                       xlab = "log response ratio (lnRR)",
                       angle = 45) +
                       scale_fill_manual(values = met.brewer("Homer1", 2)) +
                       scale_colour_manual(values = met.brewer("Homer1", 2))

p2_all
ggsave("Treatment_all_11Aug.pdf", dpi = 450)

##########################
# 2. Diameter of pattern #
##########################
mr_all1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Diameter_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#   logLik   Deviance        AIC        BIC       AICc   
# -258.0192   516.0385   526.0385   543.8611   526.2738

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0858  0.2929     32     no           Study_ID 
# sigma^2.2  0.0207  0.1439     88     no  Shared_control_ID 
# sigma^2.3  0.2433  0.4933    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6296.6831, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 0.9909, p-val = 0.3204

# Model Results:

#                   estimate      se    tval   df    pval    ci.lb   ci.ub    
# intrcpt             0.1174  0.1152  1.0189  261  0.3092  -0.1094  0.3442    
# Diameter_pattern    0.0103  0.0104  0.9954  261  0.3204  -0.0101  0.0308 

r2_ml(mr_all1)
# R2_marginal R2_conditional 
# 0.02391305     0.77083869

bubble_plot(mr_all1,
            mod = "Diameter_pattern",
            group = "Study_ID",
            k = TRUE, g = TRUE,
            xlab = "Diameter (mm)")
# plotしたいんだけどなぁ
# FIXME - Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
# replacement has 0 rows, data has 261

mr_all1_1 <- rma.mv(yi = lnRR,
                    V = lnRR_var, 
                    mods = ~ Diameter_pattern + Treatment_stimulus,
                    random = list(~1 | Study_ID,
                                  ~1 | Shared_control_ID,
                                  ~1 | Obs_ID),
                    test = "t",
                    method = "REML",
                    sparse = TRUE,
                    data = dat)

summary(mr_all1_1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#   logLik   Deviance        AIC        BIC       AICc   
# -256.9127   513.8255   525.8255   547.1896   526.1575   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0923  0.3038     32     no           Study_ID 
# sigma^2.2  0.0208  0.1442     88     no  Shared_control_ID 
# sigma^2.3  0.2428  0.4928    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 6295.6998, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 0.7978, p-val = 0.4514


# Model Results:
# 
#                            estimate      se    tval   df    pval    ci.lb 
# intrcpt                      0.0420  0.1511  0.2776  260  0.7815  -0.2556 
# Diameter_pattern             0.0111  0.0106  1.0470  260  0.2961  -0.0098 
# Treatment_stimuluseyespot    0.0970  0.1237  0.7836  260  0.4340  -0.1467 
#                             ci.ub    
# intrcpt                    0.3395    
# Diameter_pattern           0.0321    
# Treatment_stimuluseyespot  0.3407 

mr_all1_2 <- rma.mv(yi = lnRR,
                    V = lnRR_var, 
                    mods = ~ Treatment_stimulus + Diameter_pattern * Area_background,
                    random = list(~1 | Study_ID,
                                  ~1 | Shared_control_ID,
                                  ~1 | Obs_ID),
                    test = "t",
                    method = "REML",
                    sparse = TRUE,
                    data = dat)

summary(mr_all1_2)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -243.8199   487.6397   503.6397   532.0634   504.2180  

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0174  0.1319     32     no           Study_ID 
# sigma^2.2  0.0112  0.1059     88     no  Shared_control_ID 
# sigma^2.3  0.2445  0.4944    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 258) = 6123.9332, p-val < .0001

# Test of Moderators (coefficients 2:5):
# F(df1 = 4, df2 = 258) = 8.2091, p-val < .0001

# Model Results:

#                                   estimate      se     tval   df    pval 
# intrcpt                            -0.1830  0.1517  -1.2065  258  0.2287 
# Diameter_pattern                    0.0506  0.0118   4.2867  258  <.0001 
# Treatment_stimuluseyespot           0.1784  0.0990   1.8013  258  0.0728 
# Area_background                    -0.0000  0.0001  -0.5012  258  0.6167 
# Diameter_pattern:Area_background   -0.0000  0.0000  -0.6007  258  0.5486 
#                                     ci.lb   ci.ub      
# intrcpt                           -0.4816  0.1157      
# Diameter_pattern                   0.0274  0.0739  *** 
# Treatment_stimuluseyespot         -0.0166  0.3734    . 
# Area_background                   -0.0002  0.0001      
# Diameter_pattern:Area_background  -0.0000  0.0000

SizeModel <- orchaRd::mod_results(mr_all1_2, group = "Study_ID",
                                  mod = "Treatment_stimulus",
                                  at = list(Diameter_pattern = c(2, 6, 10, 14)), by = "Diameter_pattern")

orchaRd::orchard_plot(SizeModel, xlab = "lnRR", angle = 45, g = FALSE,
                      condition.lab = "Diameter difference (mm)") + 
                      theme(legend.direction = "vertical") +
                      scale_fill_manual(values = met.brewer("Klimt")) +
                      scale_colour_manual(values = met.brewer("Klimt"))
ggsave("Size_all_18Aug.pdf", dpi = 450)

#######################
# 3. area of pattern #
######################
mr_all2 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all2)

# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -258.2197   516.4394   526.4394   544.2620   526.6747   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0896  0.2994     32     no           Study_ID 
# sigma^2.2  0.0210  0.1449     88     no  Shared_control_ID 
# sigma^2.3  0.2436  0.4935    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6254.4754, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 0.2757, p-val = 0.6000

# Model Results:
#               estimate      se    tval   df    pval    ci.lb   ci.ub    
# intrcpt         0.1836  0.0856  2.1463  261  0.0328   0.0152  0.3521  * 
# Area_pattern    0.0003  0.0005  0.5251  261  0.6000  -0.0008  0.0013    

r2_ml(mr_all2)
# R2_marginal R2_conditional 
# 0.007728165    0.317609328 

bubble_plot(mr_all2,
            mod = "Area_pattern",
            group = "Study_ID",
            k = TRUE, g = TRUE,
            xlab = "Area (mm2)")
# plotしたいんだけどなぁ
# FIXME - Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
# replacement has 0 rows, data has 263

mr_all2_1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_pattern + Treatment_stimulus,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all2_1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -257.2296   514.4592   526.4592   547.8232   526.7912   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0941  0.3067     32     no           Study_ID 
# sigma^2.2  0.0215  0.1466     88     no  Shared_control_ID 
# sigma^2.3  0.2435  0.4934    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 6253.9101, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 0.3425, p-val = 0.7103

# Model Results:

#                           estimate      se    tval   df    pval    ci.lb 
# intrcpt                      0.1302  0.1202  1.0829  260  0.2799  -0.1065 
# Area_pattern                 0.0002  0.0005  0.4318  260  0.6663  -0.0008 
# Treatment_stimuluseyespot    0.0803  0.1251  0.6423  260  0.5213  -0.1660 
#                             ci.ub    
# intrcpt                    0.3669    
# Area_pattern               0.0013    
# Treatment_stimuluseyespot  0.3266    

mr_all2_2 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Treatment_stimulus + Area_pattern * Area_background,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all2_2)
#    Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -244.9117   489.8233   505.8233   534.2470   506.4017  

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0371  0.1925     32     no           Study_ID 
# sigma^2.2  0.0151  0.1228     88     no  Shared_control_ID 
# sigma^2.3  0.2414  0.4913    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 258) = 6192.3803, p-val < .0001

# Test of Moderators (coefficients 2:5):
# F(df1 = 4, df2 = 258) = 6.2146, p-val < .0001

# Model Results:

#                                estimate      se     tval   df    pval    ci.lb 
# intrcpt                          0.0854  0.1201   0.7112  258  0.4776  -0.1510 
# Area_pattern                     0.0028  0.0008   3.5611  258  0.0004   0.0012 
# Treatment_stimulus eyespot       0.0948  0.1081   0.8770  258  0.3813  -0.1180 
# Area_background                 -0.0000  0.0001  -0.5968  258  0.5512  -0.0002 
# Area_pattern:Area_background    -0.0000  0.0000  -1.1223  258  0.2628  -0.0000 
#                                 ci.ub      
# intrcpt                        0.3218      
# Area_pattern                   0.0043  *** 
# Treatment_stimulus eyespot     0.3076      
# Area_background                0.0001      
# Area_pattern:Area_background   0.0000

AreaModel <- orchaRd::mod_results(mr_all2_2, group = "Study_ID",
                                  mod = "Treatment_stimulus", 
                                  at = list(Area_pattern = c(40, 80, 120, 160)), by = "Area_pattern")

orchaRd::orchard_plot(AreaModel, xlab = "lnRR", angle = 45, g = FALSE,
                      condition.lab = "Area difference (mm²)") + 
                      theme(legend.direction = "vertical") +
                      scale_fill_manual(values = met.brewer("Redon")) +
                      scale_colour_manual(values = met.brewer("Redon"))

ggsave("Area_all_18Aug.pdf", dpi = 450)


########################
# 4. Number of pattern #
########################
mr_all3 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Number_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all3)

#    logLik   Deviance        AIC        BIC       AICc   
# -257.2227   514.4455   526.4455   547.8326   526.7762   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0772  0.2778     32     no           Study_ID 
# sigma^2.2  0.0195  0.1395     88     no  Shared_control_ID 
# sigma^2.3  0.2408  0.4907    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6440.4382, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 4.1797, p-val = 0.0419

# Model Results:

#                 estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt           0.3521  0.0998   3.5284  261  0.0005   0.1556   0.5486  *** 
# Number_pattern   -0.0601  0.0294  -2.0444  261  0.0419  -0.1180  -0.0022    * 

r2_ml(mr_all3)
# R2_marginal R2_conditional 
#  0.01793765     0.29920050 

bubble_plot(mr_all3,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number")
# Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
# replacement has 0 rows, data has 263

mr_all3_1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Number_pattern + Treatment_stimulus,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -256.0600   512.1201   526.1201   551.0448   526.5645   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0805  0.2838     32     no           Study_ID 
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0197  0.1405     88     no  Shared_control_ID 
# sigma^2.4  0.2405  0.4905    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 6363.3087, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 2.4450, p-val = 0.0887

# Model Results:

#                             estimate      se     tval   df    pval    ci.lb 
# intrcpt                       0.2836  0.1289   2.2001  260  0.0287   0.0298 
# Number_pattern               -0.0619  0.0295  -2.0949  260  0.0371  -0.1201 
# Treatment_stimulus eyespot    0.1032  0.1207   0.8552  260  0.3932  -0.1344 
#                               ci.ub    
# intrcpt                      0.5375  * 
# Number_pattern              -0.0037  * 
# Treatment_stimulus eyespot   0.3408

NumModel <- orchaRd::mod_results(mr_all3_1, group = "Study_ID",
                                  mod = "Treatment_stimulus", 
                                  at = list(Number_pattern = c(1, 2, 3)), by = "Number_pattern")

orchaRd::orchard_plot(NumModel, xlab = "lnRR", angle = 45, g = FALSE,
                      condition.lab = "Number difference ") + 
                      theme(legend.direction = "vertical") +
                      scale_fill_manual(values = met.brewer("Cross", direction = -1)) +
                      scale_colour_manual(values = met.brewer("Cross", direction = -1))

ggsave("Number_all_11Aug.pdf", dpi = 450)

# interaction
mr_all3_2 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Treatment_stimulus + Number_pattern * Area_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all3_2)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -253.0013   506.0025   522.0025   550.4262   522.5808   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.1042  0.3228     32     no           Study_ID 
# sigma^2.2  0.0211  0.1451     88     no  Shared_control_ID 
# sigma^2.3  0.2354  0.4852    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 258) = 6181.2881, p-val < .0001

# Test of Moderators (coefficients 2:5):
# F(df1 = 4, df2 = 258) = 1.7700, p-val = 0.1353

# Model Results:

#                              estimate      se     tval   df    pval    ci.lb 
# intrcpt                         0.2872  0.1429   2.0099  258  0.0455   0.0058 
# Number_pattern                 -0.0723  0.0310  -2.3307  258  0.0205  -0.1334 
# Treatment_stimulus eyespot      0.1048  0.1261   0.8310  258  0.4068  -0.1435 
# Area_pattern                   -0.0015  0.0012  -1.2498  258  0.2125  -0.0038 
# Number_pattern:Area_pattern     0.0009  0.0006   1.4991  258  0.1351  -0.0003 
#                                 ci.ub    
# intrcpt                        0.5685  * 
# Number_pattern                -0.0112  * 
# Treatment_stimulus eyespot     0.3531    
# Area_pattern                   0.0008    
# Number_pattern:Area_pattern    0.00208

###################
# 5. type of prey #
###################
mr_all4 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Type_prey -1,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all4)
#     logLik   Deviance        AIC        BIC       AICc   
# -258.7741   517.5482   527.5482   545.3708   527.7835   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0875  0.2959     32     no           Study_ID 
# sigma^2.2  0.0228  0.1511     88     no  Shared_control_ID 
# sigma^2.3  0.2429  0.4929    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6465.5593, p-val < .0001

# Test of Moderators (coefficients 1:2):
# F(df1 = 2, df2 = 261) = 4.0268, p-val = 0.0189

# Model Results:

#                      estimate      se    tval   df    pval    ci.lb   ci.ub    
# Type_preyartificial    0.1932  0.0909  2.1250  261  0.0345   0.0142  0.3723  * 
# Type_preyreal          0.2298  0.1222  1.8810  261  0.0611  -0.0108  0.4704  . 

orchard_plot(mr_all4,
             mod = "Type_prey",
             group = "Study_ID",
             xlab = "Type of prey",
             angle = 45)

ggsave("type_prey_all.pdf", dpi = 450)

# 6. shape of prey
mr_all5 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_prey -1,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all5)
#    logLik   Deviance        AIC        BIC       AICc   
# -255.5402   511.0804   525.0804   549.9782   525.5266  

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.1009  0.3177     32     no           Study_ID 
# sigma^2.2  0.0268  0.1636     88     no  Shared_control_ID 
# sigma^2.3  0.2386  0.4884    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 259) = 6128.2304, p-val < .0001

# Test of Moderators (coefficients 1:4):
# F(df1 = 4, df2 = 259) = 2.1759, p-val = 0.0721

# Model Results:
#                                estimate      se    tval   df    pval    ci.lb 
# Shape_prey abstract_butterfly    0.3025  0.1390  2.1766  259  0.0304   0.0288 
# Shape_prey abstract_stimuli      0.0135  0.2346  0.0577  259  0.9541  -0.4485 
# Shape_prey butterfly             0.2324  0.1268  1.8336  259  0.0679  -0.0172 
#Shape_prey caterpillar            0.1245  0.1606  0.7750  259  0.4391  -0.1918 
#                                ci.ub    
#                               0.5762  * 
#                               0.4756    
#                               0.4820  . 
#                               0.4407    
           
orchard_plot(mr_all5,
             mod = "Shape_prey",
             group = "Study_ID",
             xlab = "Shape of prey", angle = 45)
ggsave("Shape_prey_all.pdf", dpi = 450)

#############################################
# 7. background and pattern characteristics #
#############################################
mr_all6 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_ratio,
                   random = list(~1 | Study_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all6)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -251.8259   502.4560   512.4560   530.2786   512.6913   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0324  0.1800     32     no           Study_ID 
# sigma^2.2  0.0175  0.1323     88     no  Shared_control_ID 
# sigma^2.3  0.2417  0.4916    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6317.1260, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 18.0485, p-val < .0001

# Model Results:

#             estimate      se     tval   df    pval    ci.lb   ci.ub      
# intrcpt      -0.0689  0.0845  -0.8157  261  0.4154  -0.2353  0.0975      
# Area_ratio    3.0507  0.7181   4.2483  261  <.0001   1.6367  4.4648  *** 

# CHECK
####################
# publication bias #
####################
# traditional funnel plot
funnel(ma_all)

df_bias <- dat %>% mutate(sqrt_inv_e_n = sqrt((Cn + Tn)/(Cn * Tn))) 

bias_model <- rma.mv(yi = lnRR,
                     V = lnRR_var, 
                     mods = ~1 + scale(sqrt_inv_e_n) + scale(Year),
                     random = list(~1 | Study_ID,
                                   ~1 | Shared_control_ID,
                                   ~1 | Obs_ID),
                     test = "t",
                     method = "REML", 
                     sparse = TRUE,
                     data = df_bias)

summary(bias_model)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -257.6678   515.3355   527.3355   548.6996   527.6675  

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0741  0.2722     32     no           Study_ID 
# sigma^2.2  0.0274  0.1654     88     no  Shared_control_ID 
# sigma^2.3  0.2436  0.4936    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 260) = 6260.2656, p-val < .0001

# Test of Moderators (coefficients 2:3):
# F(df1 = 2, df2 = 260) = 0.1597, p-val = 0.8525

# Model Results:

#                      estimate      se     tval   df    pval    ci.lb   ci.ub 
# intrcpt                0.2026  0.0703   2.8820  260  0.0043   0.0642  0.3411 ** 
# scale(sqrt_inv_e_n)   -0.0256  0.0643  -0.3976  260  0.6913  -0.1521  0.1010 
# scale(Year)           -0.0333  0.0676  -0.4923  260  0.6229  -0.1663  0.0998 