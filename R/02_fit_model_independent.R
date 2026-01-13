# R/02_fit_model_independent.R
# Model 1: Trend + deterministic seasonality, independent Gaussian errors
# Saves a stanfit object for downstream PSIS-LOO and inference.

source("R/00_packages.R")

dir.create("output/fits", recursive = TRUE, showWarnings = FALSE)

# ---- Load prepared data ----
prep <- readRDS("output/prepared_data.rds")
data_m1 <- prep$stan_data_m1

# ---- Stan model ----
stan_file <- "stan/m1_harmonic_only.stan"
if (!file.exists(stan_file)) stop("Stan file not found: ", stan_file)

# ---- Fit args ----
fit_args <- list(
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 1234,
  refresh = 200,
  control = list(adapt_delta = 0.99, max_treedepth = 18)
)

cat("Fitting Model 1 at:", format(Sys.time()), "\n")

# ---- Fit model (stanfit) ----
fit_m1 <- do.call(
  rstan::stan,
  c(list(file = stan_file, data = data_m1), fit_args)
)

cat("Stan returned at:", format(Sys.time()), "\n")

# ---- Diagnostics ----
rstan::check_hmc_diagnostics(fit_m1)

# ---- Save stanfit ----
out_fit <- "output/fits/fit_m1.rds"
saveRDS(fit_m1, out_fit, compress = TRUE)
cat("Saved:", out_fit, "\n")
