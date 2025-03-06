#'
#' STAN definitions of BWQS models for continuous outcome Y.
#'
bwqs_covmod_contY <- "
data {
  int<lower=1> N;       // total number of observations
  int<lower=2> K;       // number of elements in the exposure mixture
  int<lower=1> M;       // number of other covariates
  matrix[N, K] X_matr;  // exposure matrix
  matrix[N, M] Z_matr;  // Other covariates matrix
  real y[N];            // continuous outcome

}

parameters {
  real intercept;           // intercept
  real coef;                // mixture index coefficient
  vector[M] delta;          // other covariates coefficients
  simplex[K] w;             // weights for the exposure mixture
  real<lower=0.0001> sig;     // standard deviation parameter
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}

transformed parameters {
  vector[N] Xb;
  Xb = intercept + coef * (X*w) + Z_matr * delta;
}

model {
  // Priors
  sig ~ inv_gamma(0.01,0.01);
  intercept ~ normal(0,10);
  coef ~ normal(0,10);
  delta ~ normal(0,10);
  Dir_alpha ~ gamma(2, 2);    // each alpha inside Dir_alpha is ~ gamma(2, 2)
  w ~ dirichlet(Dir_alpha);
  y ~ normal(Xb, sig);
}
generated quantities{
  # calculate likelihood
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = normal_lpdf(y[n] | Xb[n], sig);
  }
}
"

bwqs_nocovmod_contY <- "
data {
  int<lower=1> N;       // total number of observations
  int<lower=2> K;       // number of elements in the exposure mixture
  matrix[N, K] X_matr;  // exposure matrix
  real y[N];            // continuous outcome
}

parameters {
  real intercept;           // intercept
  real coef;                // mixture index coefficient
  simplex[K] w;             // weights for the exposure mixture
  real<lower=0.0001> sig;     // standard deviation parameter
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}

transformed parameters {
  vector[N] Xb;
  Xb = intercept + coef * (X*w);
}

model {
  // Priors
  sig ~ inv_gamma(0.01,0.01);
  intercept ~ normal(0,10);
  coef ~ normal(0,10);
  Dir_alpha ~ gamma(2, 2);    // each alpha inside Dir_alpha is ~ gamma(2, 2)
  w ~ dirichlet(Dir_alpha);
  y ~ normal(Xb, sig);
}
generated quantities{
  # calculate likelihood
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = normal_lpdf(y[n] | Xb[n], sig);
  }
}
"

#'
#' STAN definitions of BWQS models for binary outcome Y.
#'
bwqs_covmod_binY <- "
data {
  int<lower=1> N;       // total number of observations
  int<lower=2> K;       // number of elements in the exposure mixture
  int<lower=1> M;       // number of other covariates
  matrix[N, K] X_matr;  // exposure matrix
  matrix[N, M] Z_matr;  // Other covariates matrix
  int<lower=0, upper=1> y[N];            // binary outcome
}

parameters {
  real intercept;           // intercept
  real coef;                // mixture index coefficient
  vector[M] delta;          // other covariates coefficients
  simplex[K] w;             // weights for the exposure mixture
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}

transformed parameters {
  vector[N] Xb;
  Xb = intercept + coef * (X*w) + Z_matr * delta;
}

model {
  // Priors
  intercept ~ normal(0,10);
  coef ~ normal(0,10);
  delta ~ normal(0,10);
  Dir_alpha ~ gamma(2, 2);    // each alpha inside Dir_alpha is ~ gamma(2, 2)
  w ~ dirichlet(Dir_alpha);
  y ~ bernoulli_logit(Xb);
}
generated quantities{
  # calculate likelihood
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = bernoulli_logit_lpmf(y[n]| Xb[n]);
  }
}
"
bwqs_nocovmod_binY <- "
data {
  int<lower=1> N;       // total number of observations
  int<lower=2> K;       // number of elements in the exposure mixture
  matrix[N, K] X_matr;  // exposure matrix
  int<lower=0, upper=1> y[N];            // binary outcome
}

parameters {
  real intercept;           // intercept
  real coef;                // mixture index coefficient
  simplex[K] w;             // weights for the exposure mixture
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}

transformed parameters {
  vector[N] Xb;
  Xb = intercept + coef * (X*w);
}

model {
  // Priors
  intercept ~ normal(0,10);
  coef ~ normal(0,10);
  Dir_alpha ~ gamma(2, 2);    // each alpha inside Dir_alpha is ~ gamma(2, 2)
  w ~ dirichlet(Dir_alpha);
  y ~ bernoulli_logit(Xb);
}
generated quantities{
  # calculate likelihood
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = bernoulli_logit_lpmf(y[n]| Xb[n]);
  }
}
"

#'
#' STAN definitions of BWQS models for count outcome Y.
#'
bwqs_covmod_countY <- "
data {
  int<lower=1> N;       // total number of observations
  int<lower=2> K;       // number of elements in the exposure mixture
  int<lower=1> M;       // number of other covariates
  matrix[N, K] X_matr;  // exposure matrix
  matrix[N, M] Z_matr;  // Other covariates matrix
  int<lower=0> y[N];            // count outcome
}

parameters {
  real intercept;           // intercept
  real coef;                // mixture index coefficient
  vector[M] delta;          // other covariates coefficients
  simplex[K] w;             // weights for the exposure mixture
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}

transformed parameters {
  vector[N] mu;
  mu = exp(intercept + coef * (X*w) + Z_matr * delta);
}

model {
  // Priors
  intercept ~ normal(0,10);
  coef ~ normal(0,10);
  delta ~ normal(0,10);
  Dir_alpha ~ gamma(2, 2);    // each alpha inside Dir_alpha is ~ gamma(2, 2)
  w ~ dirichlet(Dir_alpha);
  y ~ poisson(mu);

}
generated quantities{
  # calculate likelihood
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = poisson_lpmf(y[n]| mu[n]);
  }
}
"
bwqs_nocovmod_countY <- "
data {
  int<lower=1> N;       // total number of observations
  int<lower=2> K;       // number of elements in the exposure mixture
  matrix[N, K] X_matr;  // exposure matrix
  int<lower=0> y[N];            // count outcome
}

parameters {
  real intercept;           // intercept
  real coef;                // mixture index coefficient
  simplex[K] w;             // weights for the exposure mixture
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}

transformed parameters {
  vector[N] mu;
  mu = exp(intercept + coef * (X*w));
}

model {
  // Priors
  intercept ~ normal(0,10);
  coef ~ normal(0,10);
  Dir_alpha ~ gamma(2, 2);    // each alpha inside Dir_alpha is ~ gamma(2, 2)
  w ~ dirichlet(Dir_alpha);
  y ~ poisson(mu);

}
generated quantities{
  # calculate likelihood
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = poisson_lpmf(y[n]| mu[n]);
  }
}
"
