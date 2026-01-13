# R/02_fit_model1.R
source("R/00_packages.R")

dir.create("output/fits", recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("output/prepared_data.rds")

data_m1 <- prep$stan_data_m1

stan_file <- "stan/m1_harmonic_only.stan"

fit_args <- list(
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 1234,
  refresh = 200,
  control = list(
    adapt_delta   = 0.99,
    max_treedepth = 18
  )
)

fit_m1 <- do.call(
  rstan::stan,
  c(list(file = stan_file, data = data_m1), fit_args)
)

# ---- keep only needed parameters ----
keep_pars <- c("alpha", "beta", "b_seas", "log_lik", "y_rep")

draws <- rstan::extract(fit_m1, pars = keep_pars, permuted = FALSE)

# Save a compact object (small + portable)
out_path <- "output/fits/fit_m1.rds"
saveRDS(
  list(
    stanfit   = fit_m1,
    draws     = draws,
    summary   = rstan::summary(fit_m1, pars = keep_pars)$summary,
    stan_data = data_m1
  ),
  out_path
)

cat("Saved:", out_path, "\n")
