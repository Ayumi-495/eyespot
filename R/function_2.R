library(metafor)
library(tidyverse)
library(orchaRd)

## calculate effect size and its variation ##

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
    mutate(
      lnRR = NA,
      lnRR_var = NA
    )

  for (i in seq_len(nrow(dt1))) {
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

      effect <- escalc(
        measure = "ROM",
        n1i = Tn, n2i = Cn,
        m1i = T_mean, m2 = C_mean,
        sd1i = T_sd, sd2i = C_sd,
        var.names = c("lnRR", "lnRR_var")
      )

      dt1$lnRR[i] <- effect$lnRR
      dt1$lnRR_var[i] <- effect$lnRR_var
    } else if (Response == "continuous" & Study_design == "dependent") {

      effect <- escalc(
        measure = "ROMC",
        ni = (Tn + Cn) / 2,
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
      T_proportion <- replace(
        T_proportion, Direction == "decrease",
        (1 - T_mean[Direction == "decrease"])
      )
      C_proportion <- replace(
        C_proportion, Direction == "decrease",
        (1 - C_mean[Direction == "decrease"])
      )

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

    } else if (Response == "proportion1" & Study_design == "dependent") {
      T_proportion <- replace(
        T_proportion, Direction == "decrease",
        (1 - T_mean[Direction == "decrease"])
      )
      C_proportion <- replace(
        C_proportion, Direction == "decrease",
        (1 - C_mean[Direction == "decrease"])
      )

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

      T_SD <- T_sd^2 / (4 * (T_proportion) * (1 - (T_proportion)))
      C_SD <- C_sd^2 / (4 * (C_proportion) * (1 - (C_proportion)))

      T_proportion <- asin_trans(T_proportion)
      C_proportion <- asin_trans(C_proportion)

      # calculate lnRR and lnRR variance
      lnRR_pro2 <- log(T_proportion / C_proportion)
      lnRR_var_pro2 <- (T_SD)^2 * (1 / (T_proportion^2 * Tn)) +
        (C_SD)^2 * (1 / (C_proportion^2 * Cn))

      dt1$lnRR[i] <- lnRR_pro2
      dt1$lnRR_var[i] <- lnRR_var_pro2

    } else if (Response == "proportion2" & Study_design == "dependent") {
      # transform proportion mean value
      asin_trans <- function(proportion) {
        trans <- asin(sqrt(proportion))
        trans
      }

      T_SD <- T_sd^2 / (4 * (T_proportion) * (1 - (T_proportion)))
      C_SD <- C_sd^2 / (4 * (C_proportion) * (1 - (C_proportion)))

      T_proportion <- asin_trans(T_proportion)
      C_proportion <- asin_trans(C_proportion)

      # calculate lnRR and lnRR variance
      lnRR_pro2 <- log(T_proportion / C_proportion)
      lnRR_var_pro2 <- (T_SD)^2 * (1 / (T_proportion^2 * Tn)) +
        (C_SD)^2 * (1 / (C_proportion^2 * Cn)) -
        2 * 0.5 * sqrt((T_SD)^2 * (1 / (T_proportion^2 * Tn))) *
          sqrt((C_SD)^2 * (1 / (C_proportion^2 * Cn)))

      dt1$lnRR[i] <- lnRR_pro2
      dt1$lnRR_var[i] <- lnRR_var_pro2
    }
  }

  return(dt1)
}

## test
library(here)
dat <- read.csv(here("data", "all_31072023.csv"))
dat_all <- effect_lnRR(dat)

dat_all$Obs_ID <- 1:nrow(dat_all)
hist(dat_all$lnRR)
hist(dat_all$lnRR_var)

dat_all %>%
  filter(lnRR > 4) %>%
  arrange(lnRR) %>%
  head(3) %>%
  select(lnRR)
#   lnRR
# 1 4.062596
# 2 4.167772

which(dat_all$lnRR == max(dat_all$lnRR))
# 9

model1 <- rma.mv(
  yi = lnRR,
  V = lnRR_var,
  random = list(
    ~ 1 | Study_ID,
    ~ 1 | Cohort_ID,
    ~ 1 | Shared_control_ID,
    ~ 1 | Obs_ID
  ),
  test = "t",
  method = "REML",
  sparse = TRUE,
  data = dat_all
)

summary(model1)
# Multivariate Meta-Analysis Model (k = 263; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -259.7358   519.4716   529.4716   547.3134   529.7060   

# Variance Components:

 #            estim    sqrt  nlvls  fixed             factor 
# sigma^2.1  0.0785  0.2802     32     no           Study_ID 
# sigma^2.2  0.0000  0.0000    157     no          Cohort_ID 
# sigma^2.3  0.0235  0.1534     88     no  Shared_control_ID 
# sigma^2.4  0.2429  0.4928    263     no             Obs_ID 

# Test for Heterogeneity:
# Q(df = 262) = 6465.9171, p-val < .0001

# Model Results:

# estimate      se    tval   df    pval   ci.lb   ci.ub     
#  0.2056  0.0707  2.9078  262  0.0040  0.0664  0.3448  ** 

orchard_plot(model1,
  group = "Study_ID",
  xlab = "log response ratio (lnRR)",
  angle = 45
)
