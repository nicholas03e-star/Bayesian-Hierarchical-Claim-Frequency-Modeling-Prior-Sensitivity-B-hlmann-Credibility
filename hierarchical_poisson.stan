data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> R;
  array[N] int<lower=0> y;
  vector[N] exposure;
  matrix[N, K] X;
  array[N] int<lower=1, upper=R> region;
  real<lower=0> sigma_scale;  // prior scale for sigma_u
}

parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma_u;
  vector[R] u_tilde;
}

transformed parameters {
  vector[R] u;
  u = sigma_u * u_tilde;
}

model {
  // priors
  alpha   ~ normal(-2, 0.5);
  beta    ~ normal(0, 1);
  sigma_u ~ normal(0, sigma_scale);
  u_tilde ~ normal(0, 1);

  // likelihood
  vector[N] log_lambda;
  for (n in 1:N)
    log_lambda[n] = alpha + X[n] * beta + u[region[n]] + log(exposure[n]);
  y ~ poisson_log(log_lambda);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = poisson_log_lpmf(y[n] | alpha + X[n] * beta + u[region[n]] + log(exposure[n]));
}