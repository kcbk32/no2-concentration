data {
  int<lower=1> N;
  vector[N] y;

  vector[N] t_std;

  int<lower=1> P;
  matrix[N, P] X_seas;

  int<lower=2> S;               // seasonal period in time steps (e.g. 365 for daily)
}
parameters {
  real alpha;
  real beta;
  vector[P] b_seas;

  real phiS_raw;                // unconstrained seasonal AR(1) coefficient
  real<lower=0> sigma_sar;       // seasonal innovation sd

  vector[N] z;                   // N(0,1) innovations (non-centered)
  real<lower=0> sigma_obs;       // observation noise
}
transformed parameters {
  real phiS = tanh(phiS_raw);    // maps to (-1, 1)
  vector[N] r;

  // Seasonal AR(1): r[t] = phiS * r[t-S] + eta[t]
  // Initialize first S states as stationary-ish:
  {
    real denom = sqrt(fmax(1 - phiS * phiS, 1e-6));
    for (t in 1:S) {
      r[t] = z[t] * sigma_sar / denom;
    }
  }
  for (t in (S+1):N) {
    r[t] = phiS * r[t - S] + z[t] * sigma_sar;
  }
}
model {
  // Regression priors
  alpha  ~ normal(0, 2);
  beta   ~ normal(0, 1);
  b_seas ~ normal(0, 0.5);

  // Seasonal AR coefficient: keeps away from |phiS| ~ 1
  phiS_raw ~ normal(0, 0.5);

  // Scale priors
  sigma_sar ~ exponential(2);    // mean 0.5
  sigma_obs ~ exponential(2);    // mean 0.5

  z ~ normal(0, 1);

  y ~ normal(alpha + beta * t_std + X_seas * b_seas + r, sigma_obs);
}
generated quantities {
  vector[N] log_lik;

  for (t in 1:N) {
    real mu = alpha + beta * t_std[t] + dot_product(X_seas[t], b_seas) + r[t];
    log_lik[t] = normal_lpdf(y[t] | mu, sigma_obs);
  }
}
