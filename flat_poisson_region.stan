data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=1> R;
  array[N] int<lower=0> y;
  vector[N] exposure;
  matrix[N, K] X;
  array[N] int<lower=1, upper=R> region;
}
parameters {
  real alpha;
  vector[K] beta;
  vector[R] gamma;
}
model {
  alpha ~ normal(-2, .5);
  beta  ~ normal(0, 1);
  gamma ~ normal(0, 1);

  vector[N] log_lambda;
  for (n in 1:N)
    log_lambda[n] = alpha + X[n] * beta + gamma[region[n]] + log(exposure[n]);

  y ~ poisson_log(log_lambda);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = poisson_log_lpmf(y[n] | alpha + X[n] * beta + gamma[region[n]] + log(exposure[n]));
}