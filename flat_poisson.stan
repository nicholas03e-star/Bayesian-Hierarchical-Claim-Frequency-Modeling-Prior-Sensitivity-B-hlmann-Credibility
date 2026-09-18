data {
  int<lower=0> N;
  int<lower=0> K;
  array[N] int<lower=0> y;
  vector[N] exposure;
  matrix[N, K] X;
}

parameters {
  real alpha;
  vector[K] beta;
}

model {
  // priors
  alpha ~ normal(-2, .5);
  beta  ~ normal(0, 1);

  // likelihood
  y ~ poisson_log(alpha + X * beta + log(exposure));
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = poisson_log_lpmf(y[n] | alpha + X[n] * beta + log(exposure[n]));
}