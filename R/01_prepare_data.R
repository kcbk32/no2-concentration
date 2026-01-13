source("R/00_packages.R")

# ---------- Paths ----------
csv_path <- "data/no2.csv"
out_dir  <- "output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Load + clean ----------
df <- read.csv(csv_path, stringsAsFactors = FALSE) %>%
  mutate(
    date = as.Date(date),
    no2  = as.numeric(no2)
  ) %>%
  filter(!is.na(date), !is.na(no2)) %>%
  arrange(date) %>%
  group_by(date) %>%                              # handle duplicates
  summarise(no2 = mean(no2, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(date = seq(min(date), max(date), by = "day")) %>%  # daily grid
  arrange(date) %>%
  filter(!is.na(no2)) %>%
  mutate(
    no2   = pmax(no2, 0),
    y_raw = log1p(no2)                             # transformed response
  )

# ---------- Standardize response ----------
y_mean <- mean(df$y_raw)
y_sd   <- sd(df$y_raw)

df <- df %>%
  mutate(
    y = as.vector((y_raw - y_mean) / y_sd)
  )

# ---------- Time index ----------
N <- nrow(df)
t_idx <- seq_len(N)
t_std <- as.vector((t_idx - mean(t_idx)) / sd(t_idx))

# ---------- Seasonality (Fourier design) ----------
period_main <- 365.25
K_fourier <- 3

X_sin <- sapply(1:K_fourier, function(k) sin(2 * pi * k * t_idx / period_main))
X_cos <- sapply(1:K_fourier, function(k) cos(2 * pi * k * t_idx / period_main))
X_seas <- cbind(X_sin, X_cos)
P <- ncol(X_seas)

# ---------- Seasonal lag for SAR model ----------
# Model 3: harmonic + seasonal AR(1)
S <- 7   # weekly seasonality for daily pollution data

# ---------- Stan data ----------
# Model 1: harmonic regression
stan_data_m1 <- list(
  N = N,
  y = df$y,
  t_std = t_std,
  P = P,
  X_seas = X_seas
)

# Model 2: harmonic + AR(1)
stan_data_m2 <- stan_data_m1

# Model 3: harmonic + seasonal AR(1)
stan_data_m3 <- c(stan_data_m1, list(S = S))

# ---------- Save ----------
saveRDS(
  list(
    df = df,
    N = N,
    t_idx = t_idx,
    t_std = t_std,
    y_mean = y_mean,
    y_sd = y_sd,
    period_main = period_main,
    K_fourier = K_fourier,
    X_seas = X_seas,
    S = S,
    stan_data_m1 = stan_data_m1,
    stan_data_m2 = stan_data_m2,
    stan_data_m3 = stan_data_m3
  ),
  file = file.path(out_dir, "prepared_data.rds")
)

cat("Saved:", file.path(out_dir, "prepared_data.rds"), "\n")
