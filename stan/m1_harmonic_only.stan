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

  real<lower=0> sigma;
}
model {
  // weakly informative priors
  alpha ~ normal(0, 2);
  beta  ~ normal(0, 1);
  b_seas ~ normal(0, 0.5);
  sigma ~ normal(0, 0.5);

  y ~ normal(alpha + beta * t_std + X_seas * b_seas, sigma);
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;

  for (t in 1:N) {
    real mu = alpha + beta * t_std[t] + dot_product(X_seas[t], b_seas);
    log_lik[t] = normal_lpdf(y[t] | mu, sigma);
    y_rep[t]   = normal_rng(mu, sigma);
  }
}
