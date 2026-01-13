# R/05_analysis_q2_q3_q4.R
source("R/00_packages.R")

dir.create("output/loo", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("output/prepared_data.rds")

# ---- Load fits (ALL are stanfit now)
fit_m1 <- readRDS("output/fits/fit_m1.rds")
fit_m2 <- readRDS("output/fits/fit_m2.rds")
fit_m3 <- readRDS("output/fits/fit_m3.rds")

# ---------- Helpers ----------

akaike_weights_from_loo <- function(loo_list_named) {
  looic <- sapply(loo_list_named, function(x) x$estimates["looic", "Estimate"])
  d <- looic - min(looic)
  w <- exp(-0.5 * d)
  w / sum(w)
}

get_log_lik_matrix <- function(fit) {
  loo::extract_log_lik(
    fit,
    parameter_name = "log_lik",
    merge_chains = TRUE
  )
}

loo_one <- function(fit, name = "Model") {
  ll <- get_log_lik_matrix(fit)
  res <- loo::loo(ll)
  
  k <- res$diagnostics$pareto_k
  if (!is.null(k) && any(k > 0.7, na.rm = TRUE)) {
    message(name, ": Pareto-k > 0.7 detected; retrying with moment_match=TRUE")
    res <- loo::loo(ll, moment_match = TRUE)
  }
  res
}

get_beta_draws <- function(fit) {
  rstan::extract(fit, pars = "beta")$beta
}


# Convert slope to per-year scale
slope_per_year_from_beta <- function(beta_draws, N) {
  sd_t <- stats::sd(seq_len(N))
  (beta_draws / sd_t) * 365.25
}

# ---------- Q2: LOOIC ----------
message("Computing LOO for Model 1 ...")
loo_m1 <- loo_one(fit_m1, "Model1")
saveRDS(loo_m1, "output/loo/loo_m1.rds"); gc()

message("Computing LOO for Model 2 ...")
loo_m2 <- loo_one(fit_m2, "Model2")
saveRDS(loo_m2, "output/loo/loo_m2.rds"); gc()

message("Computing LOO for Model 3 ...")
loo_m3 <- loo_one(fit_m3, "Model3")
saveRDS(loo_m3, "output/loo/loo_m3.rds"); gc()

loo_list <- list(
  Model1 = loo_m1,
  Model2 = loo_m2,
  Model3 = loo_m3
)

loo_table <- data.frame(
  model = names(loo_list),
  looic = sapply(loo_list, function(x) x$estimates["looic", "Estimate"]),
  se    = sapply(loo_list, function(x) x$estimates["looic", "SE"])
) |>
  dplyr::arrange(.data$looic)

print(loo_table)
saveRDS(loo_table, "output/results/loo_table.rds")

# ---------- Q3: Akaike weights ----------
w <- akaike_weights_from_loo(loo_list)

weights_df <- data.frame(
  model = names(w),
  akaike_weight = as.numeric(w)
) |>
  dplyr::arrange(dplyr::desc(.data$akaike_weight))

print(weights_df)
saveRDS(weights_df, "output/results/akaike_weights.rds")

# ---------- Q4: Trend from best model ----------
best_model_name <- loo_table$model[1]
cat("\nBest model by LOOIC:", best_model_name, "\n")

best_fit <- switch(
  best_model_name,
  Model1 = fit_m1,
  Model2 = fit_m2,
  Model3 = fit_m3
)

beta_draws <- get_beta_draws(best_fit)

slope_year <- slope_per_year_from_beta(beta_draws, N = prep$N)
trend_summary <- stats::quantile(slope_year, probs = c(0.05, 0.5, 0.95))

pct_change <- (exp(slope_year) - 1) * 100
pct_summary <- stats::quantile(pct_change, probs = c(0.05, 0.5, 0.95))

cat("\nTrend (dy/year in log1p scale) [5%, 50%, 95%]:\n")
print(trend_summary)

cat("\nApprox annual % change in (1+NO2) [5%, 50%, 95%]:\n")
print(pct_summary)

saveRDS(
  list(
    loo_table            = loo_table,
    akaike_weights       = weights_df,
    best_model           = best_model_name,
    slope_year_summary   = trend_summary,
    pct_change_summary   = pct_summary
  ),
  "output/results/analysis_results.rds"
)

cat("\nSaved: output/results/analysis_results.rds\n")
