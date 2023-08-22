# test whether phylogeny should be included in our analyses

# read libraries
pacman::p_load(ape, here, tidyverse, miWQS, parallel)
source("R/function_2.R")

# data and phylogeny
dat_predator <- read_csv(here("data/predator_22072023.csv"))
trees <- read.nexus(here("data/bird_phy.nex"))
is.ultrametric(trees[[1]])
# [1] TRUE

# example of variance-covariance matrix
vcv(trees[[1]], cor = T)
#                     Gallus_gallus Cyanocitta_cristata Sturnus_vulgaris Ficedula_hypoleuca Emberiza_sulphurata Parus_major Parus_caeruleus
# Gallus_gallus                   1           0.0000000        0.0000000          0.0000000           0.0000000   0.0000000       0.0000000
# Cyanocitta_cristata             0           1.0000000        0.4437522          0.4437522           0.4437522   0.4437522       0.4437522
# Sturnus_vulgaris                0           0.4437522        1.0000000          0.6390457           0.5328950   0.5288567       0.5288567
# Ficedula_hypoleuca              0           0.4437522        0.6390457          1.0000000           0.5328950   0.5288567       0.5288567
# Emberiza_sulphurata             0           0.4437522        0.5328950          0.5328950           1.0000000   0.5288567       0.5288567
# Parus_major                     0           0.4437522        0.5288567          0.5288567           0.5288567   1.0000000       0.8176725
# Parus_caeruleus                 0           0.4437522        0.5288567          0.5288567           0.5288567   0.8176725       1.0000000


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

tree_50 <- trees[1:50]
vcv_tree_50 <- map(tree_50, ~vcv(.x, corr = TRUE))

# running 50 meta-analyses with 50 different trees
ma_50 <- mclapply(vcv_tree_50, phy_model, mc.cores = 8) # detectCores() = 8

# It is not necessary for me - save and load the results
# saveRDS(ma_50, here("data", "ma_50.RDS")) 
# ma_50 <- readRDS(here("data", "ma_50.RDS"))

# combining the results
est_50 <- map_dbl(ma_50, ~ .x$b[[1]])
se_50 <-  map_dbl(ma_50, ~ .x$se)
df_50 <- c(rbind(est_50, se_50))

# creating an array required by pool.mi
my_array <- array(df_50, dim = c(1, 2, 50))

pool.mi(my_array)
# - Pooling estimates from 50 imputed analyses for 1 parameters -
#   pooled.mean pooled.total.se frac.miss.info   CI.1  CI.2 P.value
# 1       0.139           0.119              0 -0.094 0.373   0.242
#   pooled.mean pooled.total.se se.within   se.between relative.inc.var
# 1   0.1393776       0.1192319 0.1192319 2.612486e-09     4.896921e-16
#   frac.miss.info        CI.1      CI.2   p.value
# 1    2.00061e-06 -0.09431294 0.3730681 0.2424192


# extract sigma^2 for averaging variance components?　- I am not sure whether this is correct *間違っている? いらない？
sigma2_mod <- do.call(rbind, lapply(ma_50, function(x) x$sigma2))
sigma2_mod <- data.frame(sigma2_mod)

colnames(sigma2_mod) <- c("sigma^2.1_Study_ID", "sigma^2.2_SharedControl_ID", 
                          "sigma^2.3_Cohort_ID", "sigma^2.4_Obs_ID", 
                          "sigma^2.5_BirdSpecies", "sigma^2.6_Phylogeny")

colMeans(sigma2_mod)
#        sigma^2.1_Study_ID sigma^2.2_SharedControl_ID 
#              3.496353e-09               9.225743e-02 
#       sigma^2.3_Cohort_ID           sigma^2.4_Obs_ID 
#              1.009451e-01               5.323456e-01 
#      sigma^2.5_BirdSpecies       sigma^2.6_Phylogeny 
#              6.604261e-10               2.888791e-10