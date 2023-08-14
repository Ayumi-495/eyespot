# read libraries
library(here)
library(MetBrewer)
library(phangorn)
library(orchaRd)

# get data
dat_all <-  read_csv(here("data/all_31072023.csv"))
dim(dat_all)


# turn all character strings to factor
dat_all <- dat_all %>%
  mutate_if(is.character, as.factor)

summary(dat_all)

# calculate lnRR and lnRR variance
source("R/function_2.R")
dat <- effect_lnRR(dat_all)
dat$Obs_ID <- 1:nrow(dat)

hist(dat$lnRR) 
hist(dat$lnRR_var)

###############
# meta-analysis#
###############
# I may exclude cohort_ID because sigma^2.2 = 0 and I2 = 0
ma_all <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  random = list(~1 | Study_ID,
                                ~1 | Cohort_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat)

summary(ma_all)

#    logLik   Deviance        AIC        BIC       AICc   
# -259.7358   519.4716   529.4716   547.3134   529.7060   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0785  0.2802     32     no           Study_ID 
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0235  0.1534     88     no  Shared_control_ID 
# sigma^2.4  0.2429  0.4928    263     no             Obs_ID 

# Test for Heterogeneity:
# Q(df = 262) = 6465.9171, p-val < .0001

# Model Results:
# estimate      se    tval   df    pval   ci.lb   ci.ub     
#   0.2056  0.0707  2.9078  262  0.0040  0.0664  0.3448  ** 

i2_all <- i2_ml(ma_all)
i2_all
#             I2_Total          I2_Study_ID         I2_Cohort_ID 
#         9.850300e+01         2.241821e+01         1.187997e-07 
# I2_Shared_control_ID            I2_Obs_ID 
#        6.722317e+00         6.936247e+01 

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

# check publication bias
funnel(ma_all)
dat$inv_n_tilda <-  with(dat, (Cn + Tn)/(Cn*Tn))
dat$sqrt_inv_n_tilda <-  with(dat, sqrt(inv_n_tilda))

publication_bias <- rma.mv(yi = lnRR,
                            V = lnRR_var, 
                            mods = ~ 1 + sqrt_inv_n_tilda,
                            random = list(~1 | Study_ID,
                                          ~1 | Cohort_ID,
                                          ~1 | Shared_control_ID,
                                          ~1 | Obs_ID),
                            test = "t",
                            method = "REML", 
                            sparse = TRUE,
                            data = dat)

summary(publication_bias)
bubble_plot(publication_bias,
            mod = "sqrt_inv_n_tilda",
            group = "Study_ID",
            xlab = "sqrt(inv_n_tilda)")
###################
# meta-regression #
###################
# 1. eyespot or conspicuous?
mr_all <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Treatment_stimulus -1,
                  random = list(~1 | Study_ID,
                                ~1 | Cohort_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat)

summary(mr_all)

#     logLik   Deviance        AIC        BIC       AICc      
#  -258.6788   517.3576   529.3576   550.7447   529.6883   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0831  0.2882     32     no           Study_ID  
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID  
# sigma^2.3  0.0234  0.1529     88     no  Shared_control_ID  
# sigma^2.4  0.2430  0.4929    263     no             Obs_ID 

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


# 2. size
mr_all1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Diameter_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all1)
#                   estimate      se    tval   df    pval    ci.lb   ci.ub    
# intrcpt             0.0725  0.1054  0.6878  259  0.4922  -0.1351  0.2801    
# Diameter_pattern    0.0091  0.0096  0.9496  259  0.3432  -0.0098  0.0280  

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
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

SizeModel <- orchaRd::mod_results(mr_all1_1, group = "Study_ID",
                                  mod = "Treatment_stimulus", 
                                  at = list(Diameter_pattern = c(2, 6, 10, 14)), by = "Diameter_pattern")

orchaRd::orchard_plot(SizeModel, xlab = "lnRR", angle = 45, g = FALSE,
                      condition.lab = "Diameter difference (mm)") + 
                      theme(legend.direction = "vertical") +
                      scale_fill_manual(values = met.brewer("Klimt")) +
                      scale_colour_manual(values = met.brewer("Klimt"))
ggsave("Size_all_11Aug.pdf", dpi = 450)


# 3. area of pattern
mr_all2 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all2)
#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0896  0.2994     32     no           Study_ID 
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0210  0.1449     88     no  Shared_control_ID 
# sigma^2.4  0.2436  0.4935    263     no             Obs_ID 

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
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all2_1)
AreaModel <- orchaRd::mod_results(mr_all2_1, group = "Study_ID",
                                  mod = "Treatment_stimulus", 
                                  at = list(Area_pattern = c(40, 60, 80, 100, 200)), by = "Area_pattern")

orchaRd::orchard_plot(AreaModel, xlab = "lnRR", angle = 45, g = FALSE,
                      condition.lab = "Area difference (mm²)") + 
                      theme(legend.direction = "vertical") +
                      scale_fill_manual(values = met.brewer("Redon")) +
                      scale_colour_manual(values = met.brewer("Redon"))

ggsave("Area_all_11Aug.pdf", dpi = 450)

# 4. number of pattern
mr_all3 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Number_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
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
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0195  0.1395     88     no  Shared_control_ID 
# sigma^2.4  0.2408  0.4907    263     no             Obs_ID 

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
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all3_1)
NumModel <- orchaRd::mod_results(mr_all3_1, group = "Study_ID",
                                  mod = "Treatment_stimulus", 
                                  at = list(Number_pattern = c(1, 2, 3)), by = "Number_pattern")

orchaRd::orchard_plot(NumModel, xlab = "lnRR", angle = 45, g = FALSE,
                      condition.lab = "Number difference (mm²)") + 
                      theme(legend.direction = "vertical") +
                      scale_fill_manual(values = met.brewer("Cross", direction = -1)) +
                      scale_colour_manual(values = met.brewer("Cross", direction = -1))

ggsave("Number_all_11Aug.pdf", dpi = 450)


# 5. type of prey
mr_all4 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Type_prey -1,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all4)
#     logLik   Deviance        AIC        BIC       AICc   
# -258.7741   517.5482   529.5482   550.9353   529.8789   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0875  0.2959     32     no           Study_ID 
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0228  0.1511     88     no  Shared_control_ID 
# sigma^2.4  0.2429  0.4929    263     no             Obs_ID 

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
             xlab = "Type of prey")

ggsave("type_prey_all.pdf", dpi = 450)

# 6. shape of prey
mr_all5 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_prey -1,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all5)
#    logLik   Deviance        AIC        BIC       AICc   
# -255.5402   511.0804   527.0804   555.5350   527.6564   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.1009  0.3177     32     no           Study_ID 
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0268  0.1636     88     no  Shared_control_ID 
# sigma^2.4  0.2386  0.4884    263     no             Obs_ID 

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

# 7. ratio of background and pattern
mr_all6 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_ratio,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all6)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -251.8259   503.6518   515.6518   537.0389   515.9825   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0303  0.1740     32     no           Study_ID 
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0216  0.1470     88     no  Shared_control_ID 
# sigma^2.4  0.2411  0.4911    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 261) = 6305.9635, p-val < .0001

# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 261) = 16.8195, p-val < .0001

# Model Results:

#             estimate      se     tval   df    pval    ci.lb   ci.ub      
# intrcpt      -0.0706  0.0861  -0.8203  261  0.4128  -0.2401  0.0989      
# Area_ratio    2.9897  0.7290   4.1012  261  <.0001   1.5542  4.4251  *** 

# interaction between area of pattern and background
mr_all6_1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_pattern * Area_background,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat)

summary(mr_all6_1)

#    logLik   Deviance        AIC        BIC       AICc   
# -246.0685   492.1369   508.1369   536.5916   508.7129   

# Variance Components:

#             estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0343  0.1853     32     no           Study_ID 
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0144  0.1201     88     no  Shared_control_ID 
# sigma^2.4  0.2424  0.4923    263     no             Obs_ID 

# Test for Residual Heterogeneity:
# QE(df = 259) = 6202.0338, p-val < .0001

# Test of Moderators (coefficients 2:4):
# F(df1 = 3, df2 = 259) = 8.0952, p-val < .0001

# Model Results:

#                               estimate      se     tval   df    pval    ci.lb 
# intrcpt                         0.1361  0.0905   1.5037  259  0.1339  -0.0421 
# Area_pattern                    0.0028  0.0008   3.7314  259  0.0002   0.0013 
# Area_background                -0.0000  0.0001  -0.4501  259  0.6530  -0.0001 
# Area_pattern:Area_background   -0.0000  0.0000  -1.2589  259  0.2092  -0.0000 
 #                               ci.ub      
# intrcpt                       0.3142      
# Area_pattern                  0.0043  *** 
# Area_background               0.0001      
# Area_pattern:Area_background  0.0000      

