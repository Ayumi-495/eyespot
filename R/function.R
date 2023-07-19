
library(easypackages)
libraries(c("metafor", "tidyverse", "gt"))


##################### create original function for predator dataset #####################
## continuous data, proportion1 (proportion_all), proportion2 (proportion_bird))

effect_lnRR_predator <- function(dt) 
{

  ## copy dataset for adding effect size and its variation (lnRR /lnRR_var) column
  dt1 <- dt %>% 
    mutate(lnRR     = NA,
           lnRR_var = NA)
  
  for(i in seq_len(nrow(dt1))) {
    
    Tn <- dt1$Tn[i]
    Cn <- dt1$Cn[i]
    T_mean <- dt1$T_mean[i]
    C_mean <- dt1$C_mean[i]
    T_sd <- dt1$T_sd[i]
    C_sd <- dt1$C_sd[i]
    Measure <- dt1$Measure[i]
    Response <- dt1$Response[i]
    
    ## continuous data - using escalc() (metafor package)
    if (Response == "continuous") {
      
      ### reverse means when latency was measured
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
    
    ## proportion data (bird treatment group vs. control group)
    else if (Response == "proportion1") {
      
      ### transform proportion value
      asin_trans <- function(proportion) {asin(sqrt(proportion))}
      
      T_mean <- asin_trans(T_mean)
      C_mean <- asin_trans(C_mean)
      
      ### calculate lnRR and lnRR variance   
      lnRR_pro1 <- log(T_mean / C_mean)
      lnRR_var_pro1 <- (1 / sqrt(8))^2 * (1 / (T_mean^2) * Tn + 
                                            1 / (C_mean^2) * Cn)
      
      dt1$lnRR[i] <- lnRR_pro1
      dt1$lnRR_var[i] <- lnRR_var_pro1
      
    }
    
    ## proportion data (within single group) 
    ## For example… 
    ## if the sample size is 30 (n = 40), 10 individuals attacked when they saw the treatment (0.25 (25%)), 
    ## but 30 individuals attacked when they saw the control (0.75 (75%)).
    
    else if (Response == "proportion2") {
      
      ### transform proportion mean value
      asin_trans <- function(proportion) { asin(sqrt(proportion)) }
      
      T_SD <- sqrt((1/sqrt(8))^2/(4*(T_mean)*(1-(T_mean))))
      C_SD <- sqrt((1/sqrt(8))^2/(4*(C_mean)*(1-(C_mean))))
      
      T_mean <- asin_trans(T_mean)
      C_mean <- asin_trans(C_mean)
      
      ### calculate lnRR and lnRR variance 
      lnRR_pro2 <- log(T_mean / C_mean)
      lnRR_var_pro2 <- (T_SD)^2 * (1 / (T_mean^2) * Tn) +
        (C_SD)^2*(1 / (C_mean^2) * Cn)
      
      dt1$lnRR[i] <- lnRR_pro2
      dt1$lnRR_var[i] <- lnRR_var_pro2
      
    }
  }
  
  return(dt1)
  
}

#####################  create original function for prey dataset #####################

effect_lnRR_prey <- function(dt) 
{
  
  ## copy dataset for adding effect size (lnRR /lnRR_var) column
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
    
    ## continuous data - using escalc() (metafor package)
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
    
    ## proportion data (treatment group vs. control group)
    else if (Response == "proportion") {
      
      ### transform proportion value
      asin_trans <- function(proportion) {asin(sqrt(proportion))}
      
      T_proportion <- asin_trans(T_proportion)
      C_proportion <- asin_trans(C_proportion)
      
      ### calculate lnRR and lnRR variance   
      lnRR_pro1 <- log(T_proportion / C_proportion)
      lnRR_var_pro1 <- (1 / sqrt(8))^2 * (1 / (T_proportion^2) * Tn + 
                                            1 / (C_proportion^2) * Cn)
      
      dt1$lnRR[i] <- lnRR_pro1
      dt1$lnRR_var[i] <- lnRR_var_pro1
      
    }
  }
  
  return(dt1)
  
}
##################### check original functions #####################

## predator dataset - import predator_19072023.csv from data folder
test_lnRR_predator  <- effect_lnRR_predator(predator)

## prey dataset - import prey_19072023.csv from data folder
test_lnRR_prey <- effect_lnRR_prey(prey)

## see tables
test_lnRR_predator %>% gt()
test_lnRR_prey     %>% gt()
