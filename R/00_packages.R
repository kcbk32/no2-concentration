# R/00_packages.R
suppressPackageStartupMessages({
  library(rstan)
  library(loo)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(splines)
})

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
