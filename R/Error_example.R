# read libraries
library(here)
library(orchaRd)
source("R/function_2.R")

# get data
dat_all <-  read_csv(here("data/all_31072023.csv"))

# calculate lnRR and lnRR variance
dat <- effect_lnRR(dat_all)
dat$Obs_ID <- 1:nrow(dat)

# meta-analysis
model1 <- rma.mv(yi = lnRR,
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

summary(model1)
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

bubble_plot(model1,
            mod = "Area_pattern",
            group = "Study_ID",
            xlab = "Area (mm2)")
# Error in `$<-.data.frame`(`*tmp*`, "condition", value = integer(0)) : 
#  replacement has 0 rows, data has 263