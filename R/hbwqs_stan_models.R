#'
#' STAN definitions for Hierarchical BWQS (HBWQS) with continuous outcome Y.
#'
#'
hbwqs_covmod_contY <- "
data {
  int<lower=0> N;     \\ total number of observations
  int<lower=0> J;     \\ number of cohorts
  int<lower=0> K;     \\ number of elements in exposure mixture
  int<lower=0> M;     \\ number of other covariates
  matrix[N,K] X_matr; \\ exposure matrix
  matrix[N,M] Z_matr; \\ covariates matrix
  int cohort[N];      \\ cohort indicator vector
  vector[N] y;        \\ continuous outcome
}
parameters {
  vector[J] a;            \\ cohort-specific intercepts
  vector[J] b;            \\ cohort-specific slopes
  real mu_a;              \\ overall intercept mean
  real<lower=0> sigma_a;  \\ overall intercept sd
  real mu_b;              \\ overall slope mean
  real<lower=0> sigma_b;  \\ overall slope sd
  vector[M] delta;        \\ covariates coefs
  simplex[K] w;           \\ exposure mixture weights
  real<lower=0> sigma;    \\ sd param for outcome model
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}
transformed parameters {
  vector[N] mu;
  mu = a[cohort] + b[cohort].*(X_matr*w) + Z_matr*delta;
}
model {
  mu_a ~ normal(0, 10);
  mu_b ~ normal(0, 10);
  sigma_a ~ inv_gamma(0.1,0.1);
  sigma_b ~ inv_gamma(0.1,0.1);
  w ~ dirichlet(Dir_alpha);
  sigma ~ inv_gamma(0.01,0.01);
  delta ~ normal(0,50);

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ normal(mu, sigma);
}
generated quantities {
  vector[N] log_lik;
  for (nn in 1:N)
    log_lik[nn] = normal_lpdf(y[nn]| mu[nn], sigma);
}
"
hbwqs_nocovmod_contY <- "
data {
  int<lower=0> N;     \\ total number of observations
  int<lower=0> J;     \\ number of cohorts
  int<lower=0> K;     \\ number of elements in exposure mixture
  matrix[N,K] X_matr; \\ exposure matrix
  int cohort[N];      \\ cohort indicator vector
  vector[N] y;        \\ continuous outcome
}
parameters {
  vector[J] a;            \\ cohort-specific intercepts
  vector[J] b;            \\ cohort-specific slopes
  real mu_a;              \\ overall intercept mean
  real<lower=0> sigma_a;  \\ overall intercept sd
  real mu_b;              \\ overall slope mean
  real<lower=0> sigma_b;  \\ overall slope sd
  simplex[K] w;           \\ exposure mixture weights
  real<lower=0> sigma;    \\ sd param for outcome model
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}
transformed parameters {
  vector[N] mu;
  mu = a[cohort] + b[cohort].*(X_matr*w);
}
model {
  mu_a ~ normal(0, 10);
  mu_b ~ normal(0, 10);
  sigma_a ~ inv_gamma(0.1,0.1);
  sigma_b ~ inv_gamma(0.1,0.1);
  w ~ dirichlet(Dir_alpha);
  sigma ~ inv_gamma(0.01,0.01);

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ normal(mu, sigma);
}
generated quantities {
  vector[N] log_lik;
  for (nn in 1:N)
    log_lik[nn] = normal_lpdf(y[nn]| mu[nn], sigma);
}
"

#'
#' STAN definitions for Hierarchical BWQS (HBWQS) with binary outcome Y.
#'
#'
hbwqs_covmod_binY <- "
data {
  int<lower=0> N;     \\ total number of observations
  int<lower=0> J;     \\ number of cohorts
  int<lower=0> K;     \\ number of elements in exposure mixture
  int<lower=0> M;     \\ number of other covariates
  matrix[N,K] X_matr; \\ exposure matrix
  matrix[N,M] Z_matr; \\ covariates matrix
  int cohort[N];      \\ cohort indicator vector
  int<lower=0, upper=1> y[N];        \\ binary outcome
}
parameters {
  vector[J] a;            \\ cohort-specific intercepts
  vector[J] b;            \\ cohort-specific slopes
  real mu_a;              \\ overall intercept mean
  real<lower=0> sigma_a;  \\ overall intercept sd
  real mu_b;              \\ overall slope mean
  real<lower=0> sigma_b;  \\ overall slope sd
  vector[M] delta;        \\ covariates coefs
  simplex[K] w;           \\ exposure mixture weights
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}
transformed parameters {
  vector[N] mu;
  mu = a[cohort] + b[cohort].*(X_matr*w) + Z_matr*delta;
}
model {
  mu_a ~ normal(0, 10);
  mu_b ~ normal(0, 10);
  sigma_a ~ inv_gamma(0.1,0.1);
  sigma_b ~ inv_gamma(0.1,0.1);
  w ~ dirichlet(Dir_alpha);
  delta ~ normal(0,50);

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ bernoulli_logit(mu);
}
generated quantities {
  vector[N] log_lik;
  for (nn in 1:N)
    log_lik[nn] = bernoulli_logit_lpmf(y[nn]| mu[nn]);
}
"

hbwqs_nocovmod_binY <- "
data {
  int<lower=0> N;     \\ total number of observations
  int<lower=0> J;     \\ number of cohorts
  int<lower=0> K;     \\ number of elements in exposure mixture
  matrix[N,K] X_matr; \\ exposure matrix
  int cohort[N];      \\ cohort indicator vector
  int<lower=0, upper=1> y[N];        \\ binary outcome
}
parameters {
  vector[J] a;            \\ cohort-specific intercepts
  vector[J] b;            \\ cohort-specific slopes
  real mu_a;              \\ overall intercept mean
  real<lower=0> sigma_a;  \\ overall intercept sd
  real mu_b;              \\ overall slope mean
  real<lower=0> sigma_b;  \\ overall slope sd
  simplex[K] w;           \\ exposure mixture weights
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}
transformed parameters {
  vector[N] mu;
  mu = a[cohort] + b[cohort].*(X_matr*w);
}
model {
  mu_a ~ normal(0, 10);
  mu_b ~ normal(0, 10);
  sigma_a ~ inv_gamma(0.1,0.1);
  sigma_b ~ inv_gamma(0.1,0.1);
  w ~ dirichlet(Dir_alpha);

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ bernoulli_logit(mu);
}
generated quantities {
  vector[N] log_lik;
  for (nn in 1:N)
    log_lik[nn] = bernoulli_logit_lpmf(y[nn]| mu[nn]);
}
"

#'
#' STAN definitions of Hierarchical BWQS (HBWQS) models for count outcome Y.
#'

hbwqs_covmod_countY <- "
data {
  int<lower=0> N;     \\ total number of observations
  int<lower=0> J;     \\ number of cohorts
  int<lower=0> K;     \\ number of elements in exposure mixture
  int<lower=0> M;     \\ number of other covariates
  matrix[N,K] X_matr; \\ exposure matrix
  matrix[N,M] Z_matr; \\ covariates matrix
  int cohort[N];      \\ cohort indicator vector
  int<lower=0> y[N];        \\ count outcome
}
parameters {
  vector[J] a;            \\ cohort-specific intercepts
  vector[J] b;            \\ cohort-specific slopes
  real mu_a;              \\ overall intercept mean
  real<lower=0> sigma_a;  \\ overall intercept sd
  real mu_b;              \\ overall slope mean
  real<lower=0> sigma_b;  \\ overall slope sd
  vector[M] delta;        \\ covariates coefs
  simplex[K] w;           \\ exposure mixture weights
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}
transformed parameters {
  vector[N] mu;
  mu = exp(a[cohort] + b[cohort].*(X_matr*w) + Z_matr*delta);
}
model {
  mu_a ~ normal(0, 10);
  mu_b ~ normal(0, 10);
  sigma_a ~ inv_gamma(0.1,0.1);
  sigma_b ~ inv_gamma(0.1,0.1);
  w ~ dirichlet(Dir_alpha);
  delta ~ normal(0,50);

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ poisson(mu);
}
generated quantities {
  vector[N] log_lik;
  for (nn in 1:N)
    log_lik[nn] = poisson_lpmf(y[nn]| mu[nn]);
}
"
hbwqs_nocovmod_countY <- "
data {
  int<lower=0> N;     \\ total number of observations
  int<lower=0> J;     \\ number of cohorts
  int<lower=0> K;     \\ number of elements in exposure mixture
  matrix[N,K] X_matr; \\ exposure matrix
  int cohort[N];      \\ cohort indicator vector
  int<lower=0> y[N];        \\ count outcome
}
parameters {
  vector[J] a;            \\ cohort-specific intercepts
  vector[J] b;            \\ cohort-specific slopes
  real mu_a;              \\ overall intercept mean
  real<lower=0> sigma_a;  \\ overall intercept sd
  real mu_b;              \\ overall slope mean
  real<lower=0> sigma_b;  \\ overall slope sd
  simplex[K] w;           \\ exposure mixture weights
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}
transformed parameters {
  vector[N] mu;
  mu = exp(a[cohort] + b[cohort].*(X_matr*w));
}
model {
  mu_a ~ normal(0, 10);
  mu_b ~ normal(0, 10);
  sigma_a ~ inv_gamma(0.1,0.1);
  sigma_b ~ inv_gamma(0.1,0.1);
  w ~ dirichlet(Dir_alpha);

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ poisson(mu);
}
generated quantities {
  vector[N] log_lik;
  for (nn in 1:N)
    log_lik[nn] = poisson_lpmf(y[nn]| mu[nn]);
}
"
