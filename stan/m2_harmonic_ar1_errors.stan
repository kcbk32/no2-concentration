data {
  int<lower=1> N;
  vector[N] y;
  vector[N] t_std;

  int<lower=1> P;
  matrix[N, P] X_seas;
}

parameters {
  real alpha;
  real beta;
  vector[P] b_seas;

  real phi_raw;              // unconstrained
  real<lower=0> sigma_eps;
}

transformed parameters {
  real phi = tanh(phi_raw);  // (-1, 1)
}

model {
  vector[N] mu;

  mu = alpha + beta * t_std + X_seas * b_seas;

  // priors
  alpha ~ normal(0, 2);
  beta  ~ normal(0, 1);
  b_seas ~ normal(0, 0.5);

  phi_raw ~ normal(0, 0.5);
  sigma_eps ~ normal(0, 0.5);

  // AR(1) errors — CONDITIONAL likelihood (stable!)
  for (t in 2:N) {
    y[t] - mu[t] ~ normal(
      phi * (y[t-1] - mu[t-1]),
      sigma_eps
    );
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] mu;

  mu = alpha + beta * t_std + X_seas * b_seas;

  // First point: weak conditional contribution
  log_lik[1] = 0;

  for (t in 2:N) {
    log_lik[t] = normal_lpdf(
      y[t] - mu[t] |
      phi * (y[t-1] - mu[t-1]),
      sigma_eps
    );
  }
}
