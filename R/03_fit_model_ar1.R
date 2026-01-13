# R/03_fit_model2.R
source("R/00_packages.R")

dir.create("output/fits", recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("output/prepared_data.rds")

# Model 2 data (must include: N, y, t_std, P, X_seas)
data_m2 <- prep$stan_data_m2

stan_file <- "stan/m2_harmonic_ar1_errors.stan"

fit_args <- list(
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 1234,
  refresh = 200,
  control = list(
    adapt_delta   = 0.99,
    max_treedepth = 20
  )
)

cat("Starting stan at:", format(Sys.time()), "\n")

# ---- Fit model ----
fit_m2 <- do.call(
  rstan::stan,
  c(list(file = stan_file, data = data_m2), fit_args)
)

cat("Stan returned at:", format(Sys.time()), "\n")

# ---- Quick diagnostics in console ----
rstan::check_hmc_diagnostics(fit_m2)
print(fit_m2, pars = c("alpha","beta","phi_raw","sigma_eps"), probs = c(.1,.5,.9))

# ---- keep only needed parameters ----
keep_pars <- c(
  "alpha","beta","b_seas",
  "phi_raw", "sigma_eps",
  "log_lik"
)

# Extract these draws (iterations x chains x parameters)
draws <- rstan::extract(fit_m2, pars = keep_pars, permuted = FALSE)

# Save compact object (+ scaling constants so one can back-transform to original y scale)
saveRDS(list(
  draws     = draws,
  summary   = rstan::summary(fit_m2, pars = keep_pars)$summary,
  stan_data = data_m2,
  prep_meta = list(
    y_mean = prep$y_mean,
    y_sd   = prep$y_sd,
    date   = if (!is.null(prep$df$date)) prep$df$date else NULL
  )
), "output/fits/fit_m2.rds", compress = TRUE)

cat("Saved: output/fits/fit_m2.rds\n")
