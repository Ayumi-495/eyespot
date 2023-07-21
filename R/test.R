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