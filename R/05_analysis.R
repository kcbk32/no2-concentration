# R/05_analysis_q2_q3_q4.R
source("R/00_packages.R")

dir.create("output/loo", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("output/prepared_data.rds")

# ---- Load fits
# Model 1: stanfit saved directly
# Model 2: compact list with $draws (permuted=FALSE) and includes log_lik + beta
# Model 3: stanfit saved directly
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

# Build (S x N) log-likelihood matrix for loo::loo()
# Supports:
# - rstan stanfit objects
# - compact list objects with $draws from rstan::extract(..., permuted = FALSE)
get_log_lik_matrix <- function(fit_obj) {
  # Case A: stanfit
  if (inherits(fit_obj, "stanfit")) {
    return(loo::extract_log_lik(
      fit_obj,
      parameter_name = "log_lik",
      merge_chains = TRUE
    ))
  }
  
  # Case B: compact list with draws array
  if (is.list(fit_obj) && !is.null(fit_obj$draws)) {
    draws <- fit_obj$draws
    dn <- dimnames(draws)
    if (is.null(dn) || length(dn) < 3 || is.null(dn[[3]])) {
      stop("Compact fit $draws has no parameter dimnames; cannot locate log_lik.")
    }
    par_names <- dn[[3]]
    
    ll_idx <- grep("^log_lik(\\[|$)", par_names)
    if (length(ll_idx) == 0) stop("No log_lik parameters found in compact draws.")
    
    ll_arr <- draws[, , ll_idx, drop = FALSE]  # iter x chain x N
    it <- dim(ll_arr)[1]
    ch <- dim(ll_arr)[2]
    N  <- dim(ll_arr)[3]
    
    ll_mat <- matrix(ll_arr, nrow = it * ch, ncol = N)
    colnames(ll_mat) <- par_names[ll_idx]
    return(ll_mat)
  }
  
  stop("Unsupported fit object: expected stanfit or list(draws=...).")
}

loo_one <- function(fit_obj, name = "Model") {
  ll <- get_log_lik_matrix(fit_obj)
  res <- loo::loo(ll)
  
  k <- res$diagnostics$pareto_k
  if (!is.null(k) && any(k > 0.7, na.rm = TRUE)) {
    message(name, ": Pareto-k > 0.7 detected; retrying with moment_match=TRUE")
    res <- loo::loo(ll, moment_match = TRUE)
  }
  res
}

# Extract beta draws as a numeric vector
# Supports:
# - stanfit (beta parameter)
# - compact list with $draws array containing beta
get_beta_draws <- function(fit_obj) {
  if (inherits(fit_obj, "stanfit")) {
    b <- rstan::extract(fit_obj, pars = "beta")$beta
    return(as.numeric(b))
  }
  
  if (is.list(fit_obj) && !is.null(fit_obj$draws)) {
    draws <- fit_obj$draws
    dn <- dimnames(draws)
    if (is.null(dn) || length(dn) < 3 || is.null(dn[[3]])) {
      stop("Compact fit $draws has no parameter dimnames; cannot locate beta.")
    }
    par_names <- dn[[3]]
    
    if (!("beta" %in% par_names)) stop("No 'beta' found in compact draws.")
    beta_arr <- draws[, , "beta"]  # iter x chain
    return(as.numeric(beta_arr))
  }
  
  stop("Unsupported fit object for beta extraction.")
}

# Convert slope to per-year scale
slope_per_year_from_beta <- function(beta_draws, N) {
  sd_t <- stats::sd(seq_len(N))
  (beta_draws / sd_t) * 365.25
}

# ---------- Q2: LOOIC ----------
message("Computing LOO for Model 1 ...")
loo_m1 <- loo_one(fit_m1, "Model1")
saveRDS(loo_m1, "output/loo/loo_m1.rds", compress = TRUE); gc()

message("Computing LOO for Model 2 ...")
loo_m2 <- loo_one(fit_m2, "Model2")
saveRDS(loo_m2, "output/loo/loo_m2.rds", compress = TRUE); gc()

message("Computing LOO for Model 3 ...")
loo_m3 <- loo_one(fit_m3, "Model3")
saveRDS(loo_m3, "output/loo/loo_m3.rds", compress = TRUE); gc()

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
saveRDS(loo_table, "output/results/loo_table.rds", compress = TRUE)

# ---------- Q3: Akaike weights ----------
w <- akaike_weights_from_loo(loo_list)

weights_df <- data.frame(
  model = names(w),
  akaike_weight = as.numeric(w)
) |>
  dplyr::arrange(dplyr::desc(.data$akaike_weight))

print(weights_df)
saveRDS(weights_df, "output/results/akaike_weights.rds", compress = TRUE)

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
  "output/results/analysis_results.rds",
  compress = TRUE
)

cat("\nSaved: output/results/analysis_results.rds\n")
