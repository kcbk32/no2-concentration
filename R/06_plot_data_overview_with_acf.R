# R/06_plot_data_overview_with_acf.R
# One informative "data overview" figure:
# (A) Raw daily NO2 + 365-day rolling mean
# (B) ACF of log1p(NO2) (model scale)

source("R/00_packages.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(gridExtra)  # for arranging plots (usually available; if not, install.packages("gridExtra"))
})

prep <- readRDS("output/prepared_data.rds")
df   <- prep$df
N    <- prep$N
stopifnot(nrow(df) == N)

# ---- Dates
dates <- NULL
if (!is.null(df$date)) {
  dates <- as.Date(df$date)
} else {
  dates <- as.Date("2000-01-01") + seq_len(N) - 1
}

# ---- Get NO2 on original scale (preferred for the top panel)
# Option A: df$no2 exists (ideal)
# Option B: reconstruct from standardized y (log1p scale) using prep$y_mean, prep$y_sd
inv_y_to_no2 <- function(y_std, y_mean, y_sd) exp(y_std * y_sd + y_mean) - 1

no2 <- NULL
if (!is.null(df$no2)) {
  no2 <- as.numeric(df$no2)
} else {
  # try reconstruct from stan_data y = standardized log1p(NO2)
  y_std <- NULL
  if (!is.null(df$y)) {
    y_std <- as.numeric(df$y)
  } else if (!is.null(prep$stan_data_m2$y)) {
    y_std <- as.numeric(prep$stan_data_m2$y)
  } else if (!is.null(prep$stan_data_m1$y)) {
    y_std <- as.numeric(prep$stan_data_m1$y)
  } else {
    stop("Cannot find df$no2 and cannot find y (standardized log1p) in df$y or stan_data.")
  }
  
  if (is.null(prep$y_mean) || is.null(prep$y_sd)) {
    stop("Need prep$y_mean and prep$y_sd to reconstruct NO2 scale.")
  }
  no2 <- inv_y_to_no2(y_std, prep$y_mean, prep$y_sd)
}

# ---- Series on model scale for ACF (log1p(NO2), standardized if that’s what you modeled)
# If you have y_std already, use it; else compute log1p(no2) and standardize for comparability.
y_std_for_acf <- NULL
if (!is.null(df$y)) {
  y_std_for_acf <- as.numeric(df$y)
} else if (!is.null(prep$stan_data_m2$y)) {
  y_std_for_acf <- as.numeric(prep$stan_data_m2$y)
} else {
  # fallback: compute from no2
  y_raw <- log1p(pmax(no2, 0))
  y_std_for_acf <- as.numeric(scale(y_raw))
}

# ---- Rolling mean (365-day centered moving average)
roll_mean_centered <- function(x, k = 365) {
  # centered moving average with NA padding, using stats::filter
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
}
no2_roll <- roll_mean_centered(no2, k = 365)

plot_df <- tibble(
  date = dates,
  no2  = no2,
  no2_roll = no2_roll
)

# ---- ACF computation (exclude NAs, keep it simple)
max_lag <- 90
acf_obj <- stats::acf(y_std_for_acf, lag.max = max_lag, plot = FALSE, na.action = na.pass)
acf_df <- tibble(
  lag = as.numeric(acf_obj$lag),
  acf = as.numeric(acf_obj$acf)
)

# Approx 95% bounds for white-noise ACF (visual guide)
ci <- 1.96 / sqrt(sum(is.finite(y_std_for_acf)))

# ------------------ Panel A: raw + rolling mean ------------------
p1 <- ggplot(plot_df, aes(date, no2)) +
  geom_point(size = 0.35, alpha = 0.12) +                 # raw points, very light
  geom_line(aes(y = no2_roll), linewidth = 1.0, na.rm = TRUE) +  # rolling mean, prominent
  labs(
    title = "Daily NO2 measurements with rolling annual mean",
    # subtitle = "Light points: daily observations. Thick line: 365-day rolling mean (trend + slow changes)."
    x = NULL, y = "NO2"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# ------------------ Panel B: ACF ------------------
p2 <- ggplot(acf_df, aes(lag, acf)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = c(-ci, ci), linetype = 2, alpha = 0.6) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.7, alpha = 0.9) +
  labs(
    title = "Autocorrelation (ACF) of transformed series",
    # subtitle = "Computed on the model scale (standardized log(1+NO2)). Persistent positive autocorrelation motivates AR errors.",
    x = "Lag (days)", y = "ACF"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# ------------------ Combine and save ------------------
dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)

out_file <- "output/plots/data_overview_no2_rollingmean_acf.png"
png(out_file, width = 1600, height = 900, res = 180)
gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(2.1, 1))
dev.off()

cat("Saved:", out_file, "\n")

# Also print to viewer
gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(2.1, 1))
