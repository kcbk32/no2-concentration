# R/04_fit_model_seasonal_ar.R
source("R/00_packages.R")

dir.create("output/fits", recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("output/prepared_data.rds")

data_m3 <- prep$stan_data_m2
if (is.null(data_m3$S)) data_m3$S <- 7

stan_file <- "stan/m3_harmonic_sar1.stan"

cat("Starting stan at:", format(Sys.time()), "\n")

fit_m3 <- rstan::stan(
  file   = stan_file,
  data   = data_m3,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 1234,
  refresh = 200,
  
  pars    = c("alpha","beta","b_seas","phiS_raw","sigma_sar","sigma_obs","log_lik"),
  include = TRUE,
  
  save_warmup = FALSE,
  control = list(adapt_delta = 0.99, max_treedepth = 18)
)

cat("Stan returned at:", format(Sys.time()), "\n")
cat("Object size (MB):", as.numeric(object.size(fit_m3))/1024^2, "\n")

cat("Starting saveRDS at:", format(Sys.time()), "\n")
saveRDS(fit_m3, "output/fits/fit_m3.rds", compress = TRUE)
cat("Saved: output/fits/fit_m3.rds\n")
