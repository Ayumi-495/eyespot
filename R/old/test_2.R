# read libraries
library(ape)
library(here)
library(metafor)
library(MetBrewer)
library(orchaRd)
library(phangorn)
library(tidyverse)

# get data
dat_prey <- read_csv(here("data/prey_22072023.csv"))
dat_pred <- read_csv(here("data/predator_22072023.csv"))
dat_all <-  read_csv(here("data/all_31072023.csv"))
dim(dat_prey)
dim(dat_pred)
dim(dat_all)

trees <- read.nexus(here("data/bird_phy.nex"))

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


# meta-analysis
dat1$Obs_ID <- 1:nrow(dat1)

ma_prey <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  random = list(~1 | Study_ID,
                                ~1 | Cohort_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat1)

summary(ma_prey)
# estimate      se    tval   df    pval   ci.lb   ci.ub      
# 0.3023  0.0775  3.8996  145  0.0001  0.1491  0.4555  *** 
  
i2_prey <- i2_ml(ma_prey)
i2_prey
# I2_Total          I2_Study_ID         I2_Cohort_ID I2_Shared_control_ID 
# 88.64257             46.61096             22.47216             19.55945 

p1_prey <-  orchard_plot(ma_prey,
             group = "Study_ID",
             xlab = "log response ratio (lnRR)", angle = 45) +
             scale_x_discrete(labels = c("Overall effect")) +
             scale_y_continuous(limit = c(-1, 2.5), breaks = seq(-1, 2.5, 0.5)) +
             scale_fill_manual(values = "paleturquoise3") +
             scale_colour_manual(values = "paleturquoise3")
p1_prey
ggsave("overall_prey_2.pdf", dpi = 450)

p1_prey_cat <- caterpillars(ma_prey, group = "Study_ID", xlab = "log response ratio (lnRR)")
ggsave("overall_cat_prey.pdf", p1_prey_cat, dpi = 450)

# meta-regression
# eyespot or conspicuous?
mr_prey <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Treatment_stimulus -1,
                  random = list(~1 | Study_ID,
                                ~1 | Cohort_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat1)

summary(mr_prey)
#                                 estimate      se    tval   df    pval   ci.lb   ci.ub      
#  Treatment_stimulusconspicuous    0.2309  0.0884  2.6109  144  0.0100  0.0561  0.4056   ** 
#  Treatment_stimuluseyespot        0.3514  0.0833  4.2192  144  <.0001  0.1868  0.5160  *** 

r2_ml(mr_prey)
# R2_marginal R2_conditional 
# 0.03003801     0.77586177 

p2_prey <- orchard_plot(mr_prey,
             mod = "Treatment_stimulus",
             group = "Study_ID",
             xlab = "log response ratio (lnRR)",
             angle = 45) +
             scale_y_continuous(limit = c(-1, 2.5), breaks = seq(-1, 2.5, 0.5)) +
             scale_fill_manual(values = met.brewer("Homer2", 2)) +
             scale_colour_manual(values = met.brewer("Homer2", 2))
ggsave("treatment_prey.pdf", p2_prey, dpi = 450)

# size
mr_prey1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Diameter_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID,
                                 ~1 | Obs_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat1)

summary(mr_prey1)
#                   estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt            -0.2060  0.0897  -2.2958  144  0.0231  -0.3833  -0.0286    * 
# Diameter_pattern    0.0709  0.0120   5.9232  144  <.0001   0.0473   0.0946  *** 

r2_ml(mr_prey1)
# R2_marginal R2_conditional 
# 0.3967114      0.6990439 
res <- mod_results(mr_prey1, group = "Study_ID")
bubble_plot(res,
            mod = "Diameter_pattern",
            group = "Study_ID",
            xlab = "Diameter (mm)")

install.packages("pacman")
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)

devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)
# plotしたいんだけどなぁ
# FIXME - Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
# replacement has 0 rows, data has 146

# area of pattern
mr_prey2 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat1)

summary(mr_prey2)
#               estimate      se     tval   df    pval    ci.lb   ci.ub      
# intrcpt        -0.0789  0.0592  -1.3322  144  0.1849  -0.1959  0.0382      
# Area_pattern    0.0094  0.0013   7.4592  144  <.0001   0.0069  0.0119  *** 

r2_ml(mr_prey2)
# R2_marginal R2_conditional 
# 0.5409066      0.7809412 

bubble_plot(mr_prey2,
             mod = "Area_pattern",
             group = "Study_ID",
             k = TRUE, g = TRUE,
             xlab = "Area (mm2)")
# plotしたいんだけどなぁ
# FIXME - Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
# replacement has 0 rows, data has 146

# number of pattern
mr_prey3 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Number_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat1)

summary(mr_prey3)
#                 estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt           0.5103  0.0896   5.6976  144  <.0001   0.3333   0.6873  *** 
# Number_pattern   -0.0672  0.0135  -4.9838  144  <.0001  -0.0938  -0.0405  *** 
  
r2_ml(mr_prey3)
#R2_marginal R2_conditional 
#0.09169542     0.78937013 

bubble_plot(mr_prey3,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number")

# type of prey
mr_prey4 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Type_prey-1,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat1)
       
summary(mr_prey4)
#                      estimate      se    tval   df    pval   ci.lb   ci.ub     
# Type_preyartificial    0.3149  0.0998  3.1560  144  0.0019  0.1177  0.5121  ** 
# Type_preyreal          0.2849  0.1408  2.0241  144  0.0448  0.0067  0.5631   * 

r2_ml(mr_prey4)
# R2_marginal R2_conditional 
# 0.001176003    0.764598307

orchard_plot(mr_prey4,
             mod = "Type_prey",
             group = "Study_ID",
             xlab = "Type of prey", 
             angle = 45)

ggsave("type_prey_prey.pdf", dpi = 450)

# shape of prey
mr_prey5 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_prey,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat1)

summary(mr_prey5)

#                        estimate      se     tval   df    pval    ci.lb   ci.ub      
# intrcpt                  0.4052  0.1200   3.3758  143  0.0009   0.1679  0.6425  *** 
# Shape_prey butterfly     -0.1219  0.1820  -0.6698  143  0.5041  -0.4817  0.2379      
# Shape_prey caterpillar   -0.2591  0.2009  -1.2899  143  0.1992  -0.6562  0.1380  

r2_ml(mr_prey5)
# R2_marginal R2_conditional 
# 0.09978477     0.77878872 

orchard_plot(mr_prey5,
             mod = "Shape_prey",
             group = "Study_ID",
             xlab = "Shape of prey")

ggsave("shape_prey_prey.pdf", dpi = 450)

##########
# predator
##########

# turn all character strings to factor
dat_pred <- dat_pred %>%
  mutate_if(is.character, as.factor)

summary(dat_pred)

dat2 <- effect_lnRR(dat_pred)
dat2$Obs_ID <- 1:nrow(dat2)

hist(dat2$lnRR) 
hist(dat2$lnRR_var)

# meta-analysis
ma_pred <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  random = list(~1 | Study_ID,
                                ~1 | Cohort_ID, 
                                ~1 | Shared_control_ID,
                                ~1 | Bird_species,
                                ~1 | Bird_species,
                                ~1 | Obs_ID),
                  R = list(Bird_species = phylo_vcv), 
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat2)

summary(ma_pred)

i2_ml(ma_pred)

orchard_plot(ma_pred,
             group = "Study_ID",
             xlab = "log response ratio (lnRR)", angle = 45) +
             scale_x_discrete(labels = c("Overall effect")) +
             scale_fill_manual(values = "darkolivegreen3") +
             scale_colour_manual(values = "darkolivegreen3")

ggsave("overall_predator.pdf", dpi = 450)

caterpillars(ma_pred, group = "Study_ID", xlab = "log response ratio (lnRR)")
ggsave("overall_cat_pred.pdf", dpi = 450)

# meta-regression
# eyespot or conspicuous?
mr_pred <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Treatment_stimulus -1,
                  random = list(~1 | Study_ID,
                                ~1 | Cohort_ID, 
                                ~1 | Shared_control_ID,
                                ~1 | Bird_species),
                  R = list(Bird_species = phylo_cor), 
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat2)

summary(mr_pred)
#                                 estimate      se     tval   df    pval    ci.lb   ci.ub    
# Treatment_stimulus conspicuous   -0.3602  0.1850  -1.9475  115  0.0539  -0.7266  0.0062  . 
# Treatment_stimulus eyespot        0.0185  0.1303   0.1420  115  0.8873  -0.2396  0.2766    

r2_ml(mr_pred)
# R2_marginal R2_conditional 
# 0.01721872     0.40342671 

orchard_plot(mr_pred,
             mod = "Treatment_stimulus",
             group = "Study_ID",
             xlab = "log response ratio (lnRR)",
             angle = 45) +
             scale_fill_manual(values = met.brewer("Tara")) +
             scale_colour_manual(values = met.brewer("Tara"))
ggsave("treatment_predator_pred.pdf", dpi = 450)

# size of pattern
mr_pred1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Diameter_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID,
                                 ~1 | Shared_control_ID,
                                 ~1 | Bird_species),
                   R = list(Bird_species = phylo_cor), 
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred1)
#                   estimate      se     tval   df    pval    ci.lb   ci.ub    
# intrcpt            -0.0972  0.1869  -0.5197  113  0.6043  -0.4675  0.2732    
# Diameter_pattern    0.0032  0.0128   0.2459  113  0.8062  -0.0222  0.0286    

r2_ml(mr_pred1)
# R2_marginal R2_conditional 
# 0.002071476    0.385630177 

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
                                 ~1 | Shared_control_ID,
                                 ~1 | Bird_species),
                   R = list(Bird_species = phylo_cor), 
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred2)
# estimate      se     tval   df    pval    ci.lb   ci.ub    
# intrcpt        -0.0376  0.1605  -0.2342  113  0.8152  -0.3555  0.2803    
# Area_pattern   -0.0002  0.0007  -0.2890  113  0.7731  -0.0017  0.0012  

r2_ml(mr_pred2)
# R2_marginal R2_conditional 
# 0.004911973    0.393166538 

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
                                 ~1 | Shared_control_ID,
                                 ~1 | Bird_species),
                   R = list(Bird_species = phylo_cor), 
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred3)
#                 estimate      se     tval   df    pval    ci.lb   ci.ub    
# intrcpt          -0.0206  0.2104  -0.0978  115  0.9222  -0.4374  0.3962    
# Number_pattern   -0.0128  0.0744  -0.1717  115  0.8640  -0.1601  0.1346    

r2_ml(mr_pred3)
# R2_marginal R2_conditional 
# 0.0002400908   0.3822988799 

bubble_plot(mr_pred3,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number")

# type of prey
mr_pred4 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Type_prey -1,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Bird_species),
                   R = list(Bird_species = phylo_cor), 
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)

summary(mr_pred4)
#                      estimate      se     tval   df    pval    ci.lb   ci.ub    
# Type_preyartificial   -0.1604  0.1655  -0.9695  115  0.3343  -0.4882  0.1673    
# Type_preyreal          0.1036  0.1948   0.5318  115  0.5959  -0.2823  0.4895  

r2_ml(mr_pred4)
# R2_marginal R2_conditional 
# 0.02774625     0.38827500 

orchard_plot(mr_pred4,
             mod = "Type_prey",
             group = "Study_ID",
             xlab = "Type of prey")
ggsave("type_prey_pred.pdf", dpi = 450)


# shape of prey
mr_pred5 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_prey -1,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID,
                                 ~1 | Bird_species),
                   R = list(Bird_species = phylo_cor), 
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat2)
# Error: Optimizer (nlminb) did not achieve convergence (convergence = 1). something is wrong

summary(mr_pred5)

r2_ml(mr_pred5)

orchard_plot(mr_pred5,
             mod = "Shape_prey",
             group = "Study_ID",
             xlab = "Shape of prey")

ggsave("shape_predator_pred.pdf", dpi = 450)

##########
# all
##########

# turn all character strings to factor
dat_all <- dat_all %>%
  mutate_if(is.character, as.factor)

summary(dat_all)

dat3$Obs_ID <- 1:nrow(dat3)

dat3 <- effect_lnRR(dat_all)

hist(dat3$lnRR) 
hist(dat3$lnRR_var)

# meta-analysis
dat1$Shared_control_ID <- 1:nrow(dat1)

ma_all <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  random = list(~1 | Study_ID,
                                ~1 | Cohort_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat3)

summary(ma_all)
# estimate      se    tval   df    pval   ci.lb   ci.ub    
#   0.1508  0.0665  2.2660  262  0.0243  0.0198  0.2818  * 

i2_all <- i2_ml(ma_all)
i2_all
#   I2_Total          I2_Study_ID         I2_Cohort_ID I2_Shared_control_ID 
#   96.62613             24.89248             24.19243             47.54122 

p1_all_test2 <-  orchard_plot(ma_all,
                        group = "Study_ID",
                        xlab = "log response ratio (lnRR)", angle = 45) +
                        scale_x_discrete(labels = c("Overall effect")) +
                        scale_fill_manual(values = met.brewer("Homer2")) +
                        scale_colour_manual(values = met.brewer("Homer2"))

p1_all_test2
ggsave("overall_all.pdf", dpi = 450)

p1_all_cat <- caterpillars(ma_all, group = "Study_ID", xlab = "log response ratio (lnRR)")

p1_all_cat
ggsave("overall_cat_all.pdf", dpi = 450)

# meta-regression
# eyespot or conspicuous?
mr_all <- rma.mv(yi = lnRR,
                  V = lnRR_var, 
                  mods = ~ Treatment_stimulus -1,
                  random = list(~1 | Study_ID,
                                ~1 | Cohort_ID,
                                ~1 | Shared_control_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dat3)

summary(mr_all)
#                                  estimate      se    tval   df    pval    ci.lb   ci.ub     
#  Treatment_stimulus conspicuous    0.0419  0.0883  0.4748  261  0.6353  -0.1320  0.2159     
#  Treatment_stimulus eyespot        0.1970  0.0735  2.6804  261  0.0078   0.0523  0.3418  ** 

r2_ml(mr_all)
#   R2_marginal R2_conditional 
#   0.02182893     0.77884570  

p2_all <- orchard_plot(mr_all,
                       mod = "Treatment_stimulus",
                       group = "Study_ID",
                       xlab = "log response ratio (lnRR)",
                       angle = 45) +
                       scale_fill_manual(values = met.brewer("Homer1", 2)) +
                       scale_colour_manual(values = met.brewer("Homer1", 2))

p2_all
ggsave("treatment_all.pdf", dpi = 450)


# size
mr_all1 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Diameter_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID,
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat3)

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

# area of pattern
mr_all2 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Area_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat3)

summary(mr_all2)
# estimate      se     tval   df    pval    ci.lb   ci.ub    
# intrcpt         0.1523  0.0803   1.8977  259  0.0588  -0.0057  0.3104  . 
# Area_pattern   -0.0001  0.0006  -0.1537  259  0.8779  -0.0012  0.0010    

r2_ml(mr_all2)
# R2_marginal R2_conditional 
# 0.001191449    0.758804003 

bubble_plot(mr_all2,
            mod = "Area_pattern",
            group = "Study_ID",
            k = TRUE, g = TRUE,
            xlab = "Area (mm2)")
# plotしたいんだけどなぁ
# FIXME - Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
# replacement has 0 rows, data has 261

# number of pattern
mr_all3 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Number_pattern,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat3)

summary(mr_all3)
# estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt           0.2962  0.0872   3.3947  261  0.0008   0.1244   0.4679  *** 
# Number_pattern   -0.0557  0.0213  -2.6139  261  0.0095  -0.0977  -0.0137   ** 


r2_ml(mr_all3)
# R2_marginal R2_conditional 
# 0.02513448     0.75011044 
bubble_plot(mr_all3,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number")
# Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
# replacement has 0 rows, data has 263

# type of prey
mr_all4 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Type_prey -1,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat3)

summary(mr_all4)
#                       estimate      se    tval   df    pval    ci.lb   ci.ub    
# Type_prey artificial    0.1380  0.0875  1.5775  261  0.1159  -0.0342  0.3102    
# Type_prey real          0.1724  0.1115  1.5470  261  0.1231  -0.0470  0.3919 

orchard_plot(mr_all4,
             mod = "Type_prey",
             group = "Study_ID",
             xlab = "Type of prey")

ggsave("type_prey_all.pdf", dpi = 450)

# shape of prey
mr_all5 <- rma.mv(yi = lnRR,
                   V = lnRR_var, 
                   mods = ~ Shape_prey -1,
                   random = list(~1 | Study_ID,
                                 ~1 | Cohort_ID, 
                                 ~1 | Shared_control_ID),
                   test = "t",
                   method = "REML", 
                   sparse = TRUE,
                   data = dat3)

summary(mr_all5)
#                                estimate      se     tval   df    pval    ci.lb   ci.ub    
# Shape_prey abstract_butterfly    0.3135  0.1215   2.5793  259  0.0105   0.0742  0.5528  * 
# Shape_prey abstract_stimuli     -0.3134  0.2340  -1.3395  259  0.1816  -0.7742  0.1473    
# Shape_prey butterfly             0.1706  0.1077   1.5840  259  0.1144  -0.0415  0.3827    
# Shape_prey caterpillar           0.0763  0.1329   0.5741  259  0.5664  -0.1855  0.3381 

orchard_plot(mr_all5,
             mod = "Shape_prey",
             group = "Study_ID",
             xlab = "Shape of prey")
ggsave("shape_prey_all.pdf", dpi = 450)
