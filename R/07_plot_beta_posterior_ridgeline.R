# R/06_plot_pctchange_posterior_ridgeline_all_models.R
# Ridgeline posterior densities of ANNUAL % change in NO2

source("R/00_packages.R")

suppressPackageStartupMessages({
  library(rstan)
  library(ggplot2)
  library(ggridges)
})

# ---- Load fits
fits <- list(
  Model1 = readRDS("output/fits/fit_m1.rds"),  # stanfit
  Model2 = readRDS("output/fits/fit_m2.rds"),  # list(draws=...)
  Model3 = readRDS("output/fits/fit_m3.rds")   # stanfit
)

# ---- Load prep (to get N for sd(t))
prep <- readRDS("output/prepared_data.rds")
N <- prep$N

steps_per_year <- 365.25
sd_t <- stats::sd(seq_len(N))  # since t_std = (t - mean(t))/sd(t)

# ---- Extract beta draws as numeric vector
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

# ---- Convert beta -> annual % change in (1+NO2)
get_pct_change_year <- function(fit_obj) {
  beta <- get_beta_draws(fit_obj)
  beta_year <- (beta / sd_t) * steps_per_year
  100 * (exp(beta_year) - 1)
}

draws_df <- do.call(rbind, lapply(names(fits), function(nm) {
  d <- get_pct_change_year(fits[[nm]])
  data.frame(model = nm, pct_change_year = as.numeric(d))
}))

# ---- Medians for vertical dashed lines
medians <- aggregate(pct_change_year ~ model, draws_df, median)
names(medians)[2] <- "median_pct"

# ---- Axis limits
xlims <- range(draws_df$pct_change_year, na.rm = TRUE)
xpad  <- 0.15 * diff(xlims)
xlims <- c(xlims[1] - xpad, xlims[2] + xpad)

# ---- Plot
p <- ggplot(draws_df, aes(x = pct_change_year, y = model)) +
  ggridges::geom_density_ridges(
    scale = 0.9,
    rel_min_height = 0.01,
    alpha = 0.65,
    linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.6) +
  geom_vline(
    data = medians,
    aes(xintercept = median_pct),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  coord_cartesian(xlim = xlims) +
  labs(
    title = "Posterior distributions of annual % change in NO2",
    x = "Annual % change in (1+NO2)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)
ggsave(
  "output/plots/all_models_pct_change_ridgeline.pdf",
  p, width = 7.5, height = 5.2
)

print(p)

# ---- Numeric summaries
summ <- aggregate(pct_change_year ~ model, draws_df, function(x) {
  c(
    q05 = unname(stats::quantile(x, 0.05)),
    median = unname(stats::median(x)),
    q95 = unname(stats::quantile(x, 0.95))
  )
})
print(summ, row.names = FALSE)
