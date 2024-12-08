library(metafor)
library(tidyverse)

## calculate effect size and its variation ##

effect_lnRR <- function(dt) {
  # replace 0 in "C_mean", "T_sd", "C_sd", "C_proportion" with each minimum values
  # if proportion too extreme value, replace minimum value (only "T_proportion")

  dt <- dt %>%
  mutate(
    C_mean = ifelse(C_mean == 0, 0.04, C_mean),
    T_sd = ifelse(T_sd == 0, 0.01, T_sd),
    C_sd = ifelse(C_sd == 0, 0.05, C_sd),
    C_proportion = ifelse(C_proportion == 0, 0.01, C_proportion),
    T_proportion = ifelse(T_proportion < 0.01, 0.01, T_proportion),
    T_proportion = ifelse(T_proportion == 1, 0.9, T_proportion)
  )
    
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
        T_proportion, Direction == "Decrease",
        (1 - T_proportion[Direction == "Decrease"])
      )
      C_proportion <- replace(
        C_proportion, Direction == "Decrease",
        (1 - C_proportion[Direction == "Decrease"])
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
        T_proportion, Direction == "Decrease",
        (1 - T_proportion[Direction == "Decrease"])
      )
      C_proportion <- replace(
        C_proportion, Direction == "Decrease",
        (1 - C_proportion[Direction == "Decrease"])
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
      
      T_proportion <- replace(
        T_proportion, Direction == "Decrease",
        (1 - T_proportion[Direction == "Decrease"])
      )
      C_proportion <- replace(
        C_proportion, Direction == "Decrease",
        (1 - C_proportion[Direction == "Decrease"])
      )
      
      # transform proportion mean value
      asin_trans <- function(proportion) {
        trans <- asin(sqrt(proportion))
        trans
      }

      T_SD <- sqrt(T_sd^2 / (4 * (T_proportion) * (1 - (T_proportion))))
      C_SD <- sqrt(C_sd^2 / (4 * (C_proportion) * (1 - (C_proportion))))

      T_proportion <- asin_trans(T_proportion)
      C_proportion <- asin_trans(C_proportion)

      # calculate lnRR and lnRR variance
      lnRR_pro2 <- log(T_proportion / C_proportion)
      lnRR_var_pro2 <- (T_SD)^2 * (1 / (T_proportion^2 * Tn)) +
        (C_SD)^2 * (1 / (C_proportion^2 * Cn))

      dt1$lnRR[i] <- lnRR_pro2
      dt1$lnRR_var[i] <- lnRR_var_pro2

    } else if (Response == "proportion2" & Study_design == "dependent") {
      
      T_proportion <- replace(
        T_proportion, Direction == "Decrease",
        (1 - T_proportion[Direction == "Decrease"])
      )
      C_proportion <- replace(
        C_proportion, Direction == "Decrease",
        (1 - C_proportion[Direction == "Decrease"])
      )
      
      # transform proportion mean value
      asin_trans <- function(proportion) {
        trans <- asin(sqrt(proportion))
        trans
      }

      T_SD <- sqrt(T_sd^2 / (4 * (T_proportion) * (1 - (T_proportion))))
      C_SD <- sqrt(C_sd^2 / (4 * (C_proportion) * (1 - (C_proportion))))

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