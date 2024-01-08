# Load packages and custom function
pacman::p_load(ape,
               GoodmanKruskal,
               DT,
               dtplyr,
               ggokabeito,
               ggcorrplot,
               here, 
               knitr,
               tidyverse,
               patchwork,
               metafor,
               miWQS,
               multcomp,
               orchaRd,
               parallel,
               RColorBrewer)

# function for calculating effect size (lnRR) - the same as function.R
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


# Dataset for analysis

# read raw data()
dat_all <-  read.csv(here("data/all.csv"), header = T, fileEncoding = "CP932")

# calculate effect sizes and their variances
dt_all <- effect_lnRR(dat_all)

# add effect size ID (obs_ID)
dt_all$Obs_ID <- 1:nrow(dt_all)

# use vcalc to calculate variance-covariance matrix
VCV <- vcalc(vi = lnRR_var, 
             cluster = Cohort_ID,
             obs = Obs_ID,
             rho = 0.5,
             data = dt_all)

# convert pattern maximum diameter/length, area total area to logarithm
dt_all$area_pattern_T  <- dt_all$Area_pattern * dt_all$Number_pattern
dt_all$Log_diameter <- log(dt_all$Diameter_pattern)
dt_all$Log_area <- log(dt_all$Area_pattern)
dt_all$Log_background <- log(dt_all$Area_background)
dt_all$Log_area_pattern_T  <- log(dt_all$Area_pattern * dt_all$Number_pattern)

dt_all <- dt_all %>%
  mutate_if(is.character, as.factor)

hist(dt_all$lnRR)
hist(dt_all$lnRR_var)


# Let's meta-analysis and meta-regressions

# Meta-analysis: overall effect size

ma_all <- rma.mv(yi = lnRR,
                  V = VCV, 
                  random = list(~1 | Study_ID,
                                ~1 | Shared_control_ID,
                                ~1 | Cohort_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dt_all)

summary(ma_all)

i2 <- round(i2_ml(ma_all), 4)
r2 <- round(r2_ml(ma_all), 4)
t_i2 <- t(i2)
t_r2 <- t(r2)
knitr::kable(t_i2)
knitr::kable(t_r2)

orchard_plot(ma_all,
              group = "Study_ID",
              xlab = "log response ratio (lnRR)", 
              angle = 45) + 
              scale_x_discrete(labels = c("Overall effect")) + 
              scale_colour_brewer(palette = "Accent") +
              scale_fill_brewer(palette = "Accent") 


# Meta-regressions: uni-moderator

## 1. pattern type

## Eyespot vs. non-eyespot patterns 
### normal model

mr_eyespot <- rma.mv(yi = lnRR,
                    V = VCV, 
                    mods = ~ Treatment_stimulus,
                    random = list(~1 | Study_ID,
                                  ~1 | Obs_ID),
                    test = "t",
                    method = "REML", 
                    sparse = TRUE,
                    data = dt_all)

summary(mr_eyespot)

r2_1 <- round(r2_ml(mr_eyespot), 4)
t_r2_1 <- t(r2_1)
knitr::kable(t_r2_1)  

### intercept-removed model 
mr_eyespot1 <- rma.mv(yi = lnRR,
                      V = VCV, 
                      mods = ~ Treatment_stimulus - 1,
                      random = list(~1 | Study_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = dt_all)

summary(mr_eyespot1)

orchard_plot(mr_eyespot,
            mod = "Treatment_stimulus",
            group = "Study_ID",
            xlab = "log response ratio (lnRR)",
            angle = 45) + 
            scale_colour_brewer(palette = "Set2") +
            scale_fill_brewer(palette = "Set2")


## 2. pattern area (mm²)

mr_area <- rma.mv(yi = lnRR,
                  V = VCV, 
                  mods = ~ Log_area,
                  random = list(~1 | Study_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dt_all)

summary(mr_area)

r2_2 <- round(r2_ml(mr_area), 4)
t_r2_2 <- t(r2_2)
knitr::kable(t_r2_2)  

bubble_plot(mr_area,
            mod = "Log_area",
            group = "Study_ID",
            xlab = "Log-transformed area")


## 3. the number of pattern shapes

mr_num <- rma.mv(yi = lnRR,
                  V = VCV,
                  mods = ~ Number_pattern,
                  random = list(~1 | Study_ID,
                                ~1 | Obs_ID),
                      test = "t",
                      method = "REML",
                      sparse = TRUE,
                      data = dt_all)

summary(mr_num)

r2_3 <- round(r2_ml(mr_num), 4)
t_r2_3 <- t(r2_3)
knitr::kable(t_r2_3)  

bubble_plot(mr_num,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number of patterns") 


## 4. prey material type
### normal model

mr_prey_type <- rma.mv(yi = lnRR,
                  V = VCV, 
                  mods = ~ Type_prey,
                  random = list(~1 | Study_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dt_all)

summary(mr_prey_type)

r2_4 <- round(r2_ml(mr_prey_type), 4)
t_r2_4 <- t(r2_4)
knitr::kable(t_r2_4)

orchard_plot(mr_prey_type,
              mod = "Type_prey",
              group = "Study_ID",
              xlab = "Type of prey",
              angle = 45) +
              scale_colour_brewer(palette = "Set3") +
              scale_fill_brewer(palette = "Set3")


### intercept-removed model

mr_prey_type1 <- rma.mv(yi = lnRR,
                        V = VCV, 
                        mods = ~ Type_prey -1,
                        random = list(~1 | Study_ID,
                                      ~1 | Obs_ID),
                        test = "t",
                        method = "REML", 
                        sparse = TRUE,
                        data = dt_all)

summary(mr_prey_type1)


## 5. diameter/length (mm) 

mr_diameter <- rma.mv(yi = lnRR,
                      V = VCV, 
                      mods = ~ Log_diameter,
                      random = list(~1 | Study_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = dt_all)

summary(mr_diameter)

r2_5 <- round(r2_ml(mr_diameter), 4)
t_r2_5 <- t(r2_5)
knitr::kable(t_r2_5)  

bubble_plot(mr_diameter,
            mod = "Log_diameter",
            group = "Study_ID",
            xlab = "Log-transformed diameter")


## 6. Total pattern area (mm²)

mr_area_T  <- rma.mv(yi = lnRR,
                  V = VCV, 
                  mods = ~ Log_area_pattern_T ,
                  random = list(~1 | Study_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML", 
                  sparse = TRUE,
                  data = dt_all)

summary(mr_area_T)

r2_6 <- round(r2_ml(mr_area_T), 4)
t_r2_6 <- t(r2_6)
knitr::kable(t_r2_6)

bubble_plot(mr_area_T,
            mod = "Log_area_pattern_T",
            group = "Study_ID",
            xlab = "Log-transformed pattern total area")
```

## 7. total area of prey surface (mm²)

mr_background <- rma.mv(yi = lnRR,
                        V = VCV, 
                        mods = ~ Log_background,
                        random = list(~1 | Study_ID,
                                      ~1 | Obs_ID),
                        test = "t",
                        method = "REML",
                        sparse = TRUE,
                        data = dt_all)

summary(mr_background)

r2_7 <- round(r2_ml(mr_background), 4)
t_r2_7 <- t(r2_7)
knitr::kable(t_r2_7)  

bubble_plot(mr_background,
            mod = "Log_background",
            group = "Study_ID",
            xlab = "Log-transformed total prey surface area")


## 8. Prey shape type

mr_prey_shape <- rma.mv(yi = lnRR,
                        V = VCV, 
                        mods = ~ Shape_prey -1,
                        random = list(~1 | Study_ID,
                                      ~1 | Obs_ID),
                                      test = "t",
                                      method = "REML", 
                                      sparse = TRUE,
                                      data = dt_all)

summary(mr_prey_shape)

r2_8 <- round(r2_ml(mr_prey_shape), 4)
t_r2_8 <- t(r2_8)
knitr::kable(t_r2_8) 

dat_all$Shape_prey <- factor(dat_all$Shape_prey)
levels(dat_all$Shape_prey)

mat_ex <- cbind(contrMat(rep(1,length(mr_prey_shape$ci.ub)), type = "Tukey"))
sig_test <- summary(glht(mr_prey_shape, linfct= mat_ex), test = adjusted("none"))

sig_test

orchard_plot(mr_prey_shape,
              mod = "Shape_prey",
              group = "Study_ID",
              xlab = "Shape of prey",
              angle = 45) +
              scale_colour_brewer(palette = "Set3") +
              scale_fill_brewer(palette = "Set3") 

# Correlation visualisation and choosing moderators
# Before we run multi-moderator meta-regressions, we need to consider the correlation between moderators. Area, diameter and background seem to be correlated. Therefore, we visualised the correlation between these variables.

## Correlation between continuous variables

corr_cont <- round(cor(dt_all[, c("Diameter_pattern", "Area_pattern", 
                    "Number_pattern","Log_area_pattern_T", "Area_background")]),2)

p_cont <- ggcorrplot(corr_cont, hc.order = TRUE, lab = TRUE, 
                    outline.col = "white", colors = c("#6D9EC1", "white", "#E46726"), 
                    title = "(a) Continuous variables")

corr_cont_log <- round(cor(dt_all[, c("Log_diameter", "Log_area", 
                    "Number_pattern","Log_area_pattern_T", "Log_background")]),2)

p_cont_log <- ggcorrplot(corr_cont_log, hc.order = TRUE, lab = TRUE, 
                    outline.col = "white", colors = c("#6D9EC1", "white", "#E46726"), 
                    title = "(b) Log-transormed continuous variables")

p_cont + p_cont_log + plot_layout(guides = 'collect')

## Correlation between continuous variables

dat1 <- dt_all %>%
  dplyr::select("Treatment_stimulus", "Type_prey", "Shape_prey")

corr_cat <- GKtauDataframe(dat1)
plot(corr_cat)

## Choose moderators
# We used model R2 values to find better models and moderator VIF values to check multicollinearity between moderators. 
# Higher R2 indicates a better model, and VIF > 2 indicates multicollinearity. 

r2_area <- rma.mv(yi = lnRR,
                  V = VCV,
                  mods = ~ Log_area + Number_pattern,
                  random = list(~1 | Study_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML",
                  sparse = TRUE,
                  data = dt_all)

r2_area_T <- rma.mv(yi = lnRR,
                  V = VCV,
                  mods = ~ Log_area_pattern_T + Number_pattern,
                  random = list(~1 | Study_ID,
                                ~1 | Obs_ID),
                  test = "t",
                  method = "REML",
                  sparse = TRUE,
                  data = dt_all)

r2_diameter <- rma.mv(yi = lnRR,
                      V = VCV,
                      mods = ~ Log_diameter + Number_pattern,
                      random = list(~1 | Study_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML",
                      sparse = TRUE,
                      data = dt_all)

# check r2 values
r2_ml(r2_area)
r2_ml(r2_area_T)
r2_ml(r2_diameter)

# Meta-regressions: multi-moderator
Third, We also ran a multi-moderator meta-regression model, including treatment stimulus pattern types, pattern area, the number of pattern shapes, and prey material type variables due to moderator correlations. 

mr_all <- rma.mv(yi = lnRR,
                V = VCV, 
                mods = ~ Treatment_stimulus + Log_area 
                + Number_pattern + Type_prey,
                random = list(~1 | Study_ID,
                              ~1 | Obs_ID),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dt_all)

summary(mr_all)

r2_9 <- round(r2_ml(mr_all), 4)
t_r2_9 <- t(r2_9)
knitr::kable(t_r2_9)  


# Publication bias

## funnel plot - standard error
funnel(ma_all, yaxis = "sei",
      xlab = "Standarised residuals",
      ylab = "Precision (inverse of SE)" )

# funnel plot - inverse of standard error
funnel(ma_all, yaxis = "seinv",
      xlab = "Standarised residuals",
      ylab = "Precision (inverse of SE)",  col = c(alpha("orange", 0.5)))


## Egger's test

df_bias <- dt_all %>% mutate(sqrt_inv_e_n = sqrt((Cn + Tn)/(Cn * Tn)))

bias_model <- rma.mv(yi = lnRR,
                      V = lnRR_var, 
                      mods = ~ 1 + sqrt_inv_e_n,
                      random = list(~1 | Study_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = df_bias)

summary(bias_model)

bubble_plot(bias_model,
            mod = "sqrt_inv_e_n",
            group = "Study_ID",
            xlab = "Square root of inverse of effective sample size")


## Decline effect

year_model <- rma.mv(yi = lnRR,
                      V = VCV, 
                      mods = ~  1 + Year,
                      random = list(~1 | Study_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = df_bias)

summary(year_model)

bubble_plot(year_model,
            mod = "Year",
            group = "Study_ID",
            xlab = "Year of publication")


## Multi-moderator

multi_bias <- rma.mv(yi = lnRR,
                      V = VCV, 
                      mods = ~ Treatment_stimulus + Log_area 
                            + Number_pattern + Type_prey
                            + sqrt_inv_e_n + Year,
                      random = list(~1 | Study_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML", 
                      sparse = TRUE,
                      data = df_bias)

summary(multi_bias)

bubble_plot(multi_bias,
            mod = "sqrt_inv_e_n",
            group = "Study_ID",
            xlab = "Square root of inverse of effective sample size")

bubble_plot(multi_bias,
            mod = "Year",
            group = "Study_ID",
            xlab = "Year of publication")


# Additional analyses


mr_dataset <- rma.mv(yi = lnRR,
                      V = VCV,
                      mods = ~ Dataset,
                      random = list(~1 | Study_ID,
                                    ~1 | Obs_ID),
                      test = "t",
                      method = "REML",
                      sparse = TRUE,
                      data = dt_all)

summary(mr_dataset)

orchard_plot(mr_dataset,
            mod = "Dataset",
            group = "Study_ID",
            xlab = "Dataset", 
            angle = 45) 


### Total pattern area

mr_all_2 <- rma.mv(yi = lnRR,
                V = VCV, 
                mods = ~ Treatment_stimulus + Log_area_pattern_T 
                + Number_pattern + Type_prey,
                random = list(~1 | Study_ID,
                              ~1 | Obs_ID),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dt_all)
summary(mr_all_2)


## Moderators associated with conspicuousness
### pattern area and number of patterns

mr_cons <- rma.mv(yi = lnRR,
                V = VCV, 
                mods = ~ Log_area + Number_pattern,
                random = list(~1 | Study_ID,
                              ~1 | Obs_ID),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dt_all)

summary(mr_cons)

r2_ml(mr_cons)

bubble_plot(mr_cons,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Log-transformed area")


### total pattern area and number of patterns

mr_cons1 <- rma.mv(yi = lnRR,
                V = VCV, 
                mods = ~ Log_area_pattern_T + Number_pattern,
                random = list(~1 | Study_ID,
                              ~1 | Obs_ID),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dt_all)

summary(mr_cons1)

r2_ml(mr_cons1)

bubble_plot(mr_cons1,
            mod = "Log_area_pattern_T",
            group = "Study_ID",
            xlab = "Log-transformed area")


### pattern area, number of patterns, and treatment stimulus

mr_cons2 <- rma.mv(yi = lnRR,
                V = VCV, 
                mods = ~ Log_area + Number_pattern + Treatment_stimulus,
                random = list(~1 | Study_ID,
                              ~1 | Obs_ID),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dt_all)

summary(mr_cons2)

r2_ml(mr_cons2)

### total pattern area, number of patterns, and treatment stimulus

mr_cons3 <- rma.mv(yi = lnRR,
                V = VCV, 
                mods = ~ Log_area_pattern_T + Number_pattern + Treatment_stimulus,
                random = list(~1 | Study_ID,
                              ~1 | Obs_ID),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dt_all)

summary(mr_cons3)

r2_ml(mr_cons3)


## Moderators associated with prey

mr_pre <- rma.mv(yi = lnRR,
                V = VCV, 
                mods = ~ Treatment_stimulus + Type_prey,
                random = list(~1 | Study_ID,
                              ~1 | Obs_ID),
                test = "t",
                method = "REML", 
                sparse = TRUE,
                data = dt_all)

summary(mr_pre)

r2_ml(mr_pre)


# More information
## Bird tree 

trees <- read.nexus(here("data/bird_phy.nex"))
plot(trees[1])


## Check phylogenetic relatedness
## We tested whether phylogeny should be included in our analyses.

# data 
dat_predator <- filter(dat_all, Dataset == "predator")

# check if the tree is ultrametric
# is.ultrametric(trees[[1]])
# [1] TRUE

# example of variance-covariance matrix
vcv(trees[[1]], corr = T)

# calculate lnRR and lnRR_var
dat_pred <- effect_lnRR(dat_predator)
dat_pred$Obs_ID <- 1:nrow(dat_pred)
V <- vcalc(lnRR_var, cluster = Cohort_ID, subgroup = Obs_ID, data = dat_pred)

# create function to run meta-analysis for 50 trees
phy_model <- function(cor_tree = vcv_tree){
            model <- rma.mv(yi = lnRR,
                            V = V,
                            random = list(~1 | Study_ID,
                                          ~1 | Shared_control_ID,
                                          ~1 | Cohort_ID,
                                          ~1 | Obs_ID,
                                          ~1 | Bird_species,
                                          ~1 | Phylogeny),
                            R = list(Phylogeny = cor_tree),
                            test = "t",
                            method = "REML",
                            sparse = TRUE,
                            data = dat_pred)
  model
}

# running 50 meta-analyses with 50 different trees
# We got 1000 phylogenetic trees from BirdTree.org (https://birdtree.org)
# Only 50 trees are used as in Nakagawa & Villemereuil (2019)(https://doi.org/10.1093/sysbio/syy089)

tree_50 <- trees[1:50]
vcv_tree_50 <- map(tree_50, ~vcv(.x, corr = TRUE))
ma_50 <- mclapply(vcv_tree_50, phy_model) # detectCores() = 12

# always best to save RDS files for analysis 
# which takes more than a min or so 
# saveRDS(ma_50, here("Rdata", "ma_50.RDS")) 

ma_50 <- readRDS(here("Rdata", "ma_50.RDS"))

# combining the results
est_50 <- map_dbl(ma_50, ~ .x$b[[1]])
se_50 <-  map_dbl(ma_50, ~ .x$se)
df_50 <- c(rbind(est_50, se_50))

# creating an array required by pool.mi
my_array <- array(df_50, dim = c(1, 2, 50))

pool.mi(my_array)

# extract sigma^2 for averaging variance components
sigma2_mod <- do.call(rbind, lapply(ma_50, function(x) x$sigma2))
sigma2_mod <- data.frame(sigma2_mod)

colnames(sigma2_mod) <- c("sigma^2.1_Study_ID", "sigma^2.2_SharedControl_ID", 
                          "sigma^2.3_Cohort_ID", "sigma^2.4_Obs_ID", 
                          "sigma^2.5_BirdSpecies", "sigma^2.6_Phylogeny")

# easier to undersatnd if you round
round(colMeans(sigma2_mod), 2)


# Figure gallery
# These were the basis for the final figures presented in the paper.

## Discrete moderators

# extract the results of each model and combine them into one data frame for plotting
main  <- mod_results(ma_all, group = "Study_ID")
pattern <- mod_results(mr_eyespot, mod = "Treatment_stimulus", group = "Study_ID")
shape  <- mod_results(mr_prey_type, mod = "Type_prey", group = "Study_ID")

comb <- submerge(main, pattern, shape)

p1 <- orchard_plot(comb,
                    group = "Study_ID",
                    xlab = "log response ratio (lnRR)",
                    tree.order = c( "Artificial", "Real","Non_eyespot", "Eyespot", "Intrcpt"),
                    angle = 45,
                    branch.size = 2) +
                    theme(legend.direction = "horizontal",
                    legend.text = element_text(size = 6)) +
                    scale_colour_brewer(palette = "Pastel1") +
                    scale_fill_brewer(palette = "Pastel1") 
p1

p2 <- orchard_plot(mr_prey_shape,
                  mod = "Shape_prey",
                  group = "Study_ID",
                  xlab = "log response ratio (lnRR)",
                  angle = 45,
                  legend.pos = "bottom.right") +
                  theme(legend.direction = "horizontal",
                  legend.text = element_text(size = 6)) +
                  scale_colour_brewer(palette = "YlGnBu") +
                  scale_fill_brewer(palette = "YlGnBu")
p2

## Continuous moderators

# these figs were created by multi meta-regression model results using the log-transformed area

p3 <- bubble_plot(mr_area,
                  mod = "Log_area",
                  group = "Study_ID",
                  xlab = "Log-transformed area",
                  est.lwd = 1.2,
                  est.col = "deeppink3") +
                  scale_x_continuous(breaks = seq(0,7,1.5)) +
                  scale_y_continuous(breaks = seq(-3.0, 4.5, 1.0))

# number of patterns
p4 <- bubble_plot(mr_num,
            mod = "Number_pattern",
            group = "Study_ID",
            xlab = "Number of patterns",
            est.lwd = 1.2,
            est.col = "turquoise4") +
      scale_x_continuous(breaks = seq(0,11,2)) +
      scale_y_continuous(breaks = seq(-3.0, 4.5, 1.0))

# combine
p3 / p4 + plot_annotation(tag_levels = "a")

# total pattern area
p5 <- bubble_plot(mr_area_T,
            mod = "Log_area_pattern_T",
            group = "Study_ID",
            xlab = "Log-transformed total pattern area")

# diameter
p6 <- bubble_plot(mr_diameter,
            mod = "Log_diameter",
            group = "Study_ID",
            xlab = "Log-transformed diameter")

# background
p7 <- bubble_plot(mr_background,
            mod = "Log_background",
            group = "Study_ID",
            xlab = "Log-transformed total prey surface area")

p5 / p6 /p7 + plot_annotation(tag_levels = "a")


## Publication bias

# funnel plot

#pdf("funnel.pdf")

funnel(ma_all, yaxis = "seinv", 
      xlab = "Standarised residuals",
      ylab = "Precision (inverse of SE)",
      # xlim = c(-4.0, 4.5), 
      # ylim = c(0.01, 60.0), 
      col = c(alpha("black", 0.5)),
      cex = 1.5) 

# pdf("funnel2.pdf", width = 15, height = 10)
# dev.off()

# Egger's test and decline effect
p8 <- bubble_plot(bias_model,
            mod = "sqrt_inv_e_n",
            group = "Study_ID",
            xlab = "Square root of inverse of effective sample size") +
            scale_y_continuous(breaks = seq(-2.5, 4.5, 1.5)) 

p9 <- bubble_plot(year_model,
            mod = "Year",
            group = "Study_ID",
            xlab = "Year of publication") +  scale_y_continuous(breaks = seq(-2.5, 4.0, 1.5)) 

# combine plots created by bubble_plot()
pub <- p8 / p9
pub +  plot_annotation(tag_levels = 'a') 


## Publication bias (multi-moderator)

p10 <- bubble_plot(multi_bias,
            mod = "sqrt_inv_e_n",
            group = "Study_ID",
            xlab = "Square root of inverse of effective sample size")

p11 <- bubble_plot(multi_bias,
            mod = "Year",
            group = "Study_ID",
            xlab = "Year of publication")

p10 / p11 +  plot_annotation(tag_levels = 'a')