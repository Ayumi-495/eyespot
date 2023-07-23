
library(metafor)
library(tidyverse)

# TODO insert - when "C_mean = 0, C_proportion = 0, C_sd = 0" replace 0.01 or minimum values using ifelse()


## function for predator dataset ##
effect_lnRR_pred <- function(dt) 
  {
  
# copy dataset for adding effect size and 
# its variation (lnRR /lnRR_var) column
  dt1 <- dt %>% 
    mutate(lnRR     = NA,
           lnRR_var = NA)
  
  for(i in seq_len(nrow(dt1))) {
    
    Tn <- dt1$Tn[i]
    Cn <- dt1$Cn[i]
    T_proportion <- dt1$T_proportion[i]
    C_proportion <- dt1$C_proportion[i]
    T_mean <- dt1$T_mean[i]
    C_mean <- dt1$C_mean[i]
    T_sd <- dt1$T_sd[i]
    C_sd <- dt1$C_sd[i]
    Response <- dt1$Response[i]
    
    # continuous data - using escalc() (metafor package)
    if (Response == "continuous") {
      
      # reverse means when latency was measured
      T_mean <- replace(T_mean, Measure == "latency",
                        T_mean[Measure == "latency"] * (-1))
      C_mean <- replace(C_mean, Measure == "latency",
                        C_mean[Measure == "latency"] * (-1))
      
      effect <- escalc(measure = "ROM", 
                       n1i = Tn, n2i = Cn, 
                       m1i = T_mean, m2 = C_mean, 
                       sd1i = T_sd, sd2i = C_sd,
                       var.names = c("lnRR", "lnRR_var")
      )
      
      dt1$lnRR[i] <- effect$lnRR
      dt1$lnRR_var[i] <- effect$lnRR_var
      
    }
    
    # proportion data (no sd values)
    else if (Response == "proportion1") {
      
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
    
   # proportion data (exist sd values) 
    else if (Response == "proportion2") {
      
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
  }
  
  return(dt1)
  
}

## function for prey dataset ##
effect_lnRR_prey <- function(dt) 
  {
  
    dt1 <- dt %>% 
    mutate(lnRR     = NA,
           lnRR_var = NA)
  
  for(i in seq_len(nrow(dt1))) {
    
    Tn <- dt1$Tn[i]
    Cn <- dt1$Cn[i]
    T_proportion <- dt1$T_proportion[i]
    C_proportion <- dt1$C_proportion[i]
    T_mean <- dt1$T_mean[i]
    C_mean <- dt1$C_mean[i]
    T_sd <- dt1$T_sd[i]
    C_sd <- dt1$C_sd[i]
    Response <- dt1$Response[i]
    
   # continuous data - using escalc() (metafor package)
    if (Response == "continuous") {
      
      effect <- escalc(measure = "ROM", 
                       n1i = Tn, n2i = Cn, 
                       m1i = T_mean, m2 = C_mean, 
                       sd1i = T_sd, sd2i = C_sd,
                       var.names = c("lnRR", "lnRR_var")
      )
      
      dt1$lnRR[i] <- effect$lnRR
      dt1$lnRR_var[i] <- effect$lnRR_var
      
    }
    
    # proportion data (no sd values)
    else if (Response == "proportion1") {
      
      # transform proportion mean value
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
    
    # proportion data (exist sd values) 
    else if (Response == "proportion2") {
      
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
  }
  
  return(dt1)
  
}

