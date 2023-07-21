# reading datasets

library(tidyverse)
library(here)
library(orchaRd)
library(metafor)

# get data

# prey
dat_prey <- read_csv(here("data/prey_19072023.csv"))
dat_pred <- read_csv(here("data/predator_19072023.csv"))
dim(dat_prey)
dim(dat_pred)

names(dat_prey)

# check whether "T_proportion" and "C_proportion" have 0 or not
dat_prey %>%
  filter(T_proportion == 0 | C_proportion == 0) %>%
  select(T_proportion, C_proportion)
# # A tibble: 3 × 2
#   T_proportion C_proportion
#          <dbl>        <dbl>
# 1        0.184            0
# 2        0.318            0
# 3        0.660            0

# check the samllest values in "C_proportion" apart from 0
dat_prey %>%
  filter(C_proportion != 0) %>%
  arrange(C_proportion) %>%
  head(3) %>%
  select(C_proportion)

#   # A tibble: 3 × 1
#   C_proportion
#          <dbl>
# 1       0.0107
# 2       0.0107
# 3       0.0183

# given this we should repliace 0 in "C_proportion" with 0.01
dat_prey <- dat_prey %>%
  mutate(C_proportion = ifelse(C_proportion == 0, 0.01, C_proportion))

# check whether "T_meam" and "C_mean" have 0 or not 
dat_prey %>%
  filter(T_mean == 0 | C_mean == 0) %>%
  select(T_mean, C_mean)

# check whether "T_sd" and "C_sd" have 0 or not
dat_prey %>%
  filter(T_sd == 0 | C_sd == 0) %>%
  select(T_sd, C_sd)
