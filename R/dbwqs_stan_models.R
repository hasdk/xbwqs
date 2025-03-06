#'
#' STAN definitions for Dirichlet-Bayesian Weighted Quantile Sum regression models.
#'
#'
#'
dbwqs_covmod <- "
functions {
  /* dirichlet-logit log-PDF
   * Args:
   *   y: vector of real response values
   *   mu: vector of category logit probabilities
   *   phi: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real dirichlet_logit_lpdf(vector y, vector mu, real phi) {
     return dirichlet_lpdf(y | softmax(mu) * phi);
   }

}


data {
  int<lower=1> N;       // total number of observations
  int<lower=2> C;       // number of elements in the compositional outcome
  int<lower=2> K;       // number of elements in the exposure mixture
  int<lower=1> M;       // number of other covariates
  vector[C] Y[N];       // compositional outcome matrix
  matrix[N, K] X_matr;  // exposure matrix
  matrix[N, M] Z_matr;  // Other covariates matrix
}

parameters {
  vector[C-1] rel_coefs;    // mixture index coefficients
  vector[C-1] intercepts; // mixture index intercepts
  real<lower=0.1> phi;  // precision parameter
  matrix[M,C-1] delta;  // covariates coefficients
  simplex[K] w;             // weights for the X mixture
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}

transformed parameters {
  vector[N] S_bwqs = X_matr * w;
  matrix[N,C-1] Zprod = Z_matr * delta;
}

model {
  // Priors
  phi ~ gamma(0.01,0.01);
  intercepts ~ normal(0,10);
  rel_coefs ~ normal(0,10);
  Dir_alpha ~ gamma(2, 2);    // each alpha inside Dir_alpha is ~ gamma(2, 2)
  w ~ dirichlet(Dir_alpha);
  to_vector(delta) ~ normal(0, 10);

  // likelihood
 {
    vector[C] mu[N];
    for (n in 1:N) {
      mu[n][1] = 0;  // Reference category
      mu[n][2:C] = intercepts + rel_coefs .* S_bwqs[n] + Zprod[n]';
      target += dirichlet_logit_lpdf(Y[n] | mu[n], phi);
    }
  }
}

// add block to calculate effect estimate of reference proportion for unit increase in exposure
generated quantities {
    vector<lower=0.0001>[N] num_i;
    vector<lower=0.0001>[N] denom_i;
    vector<lower=0.0001>[N] ref_prop_i;
    vector<lower = 0.0001>[C] abs_coefs;

    for(n in 1:N){
      num_i[n] = 1.0;
      denom_i[n] = 1.0;
      for(c in 1:(C-1)){
        num_i[n] += exp(intercepts[c] + rel_coefs[c]*S_bwqs[n] + Zprod[n,c]);
        denom_i[n] += exp(intercepts[c] + rel_coefs[c]* (S_bwqs[n]+1)  + Zprod[n,c]);
      }
      ref_prop_i[n] = num_i[n] / denom_i[n];
    }
    abs_coefs[1] = mean(ref_prop_i);
    abs_coefs[2:C] = exp(rel_coefs) * abs_coefs[1];
}
"

dbwqs_nocovmod <- "
functions {
  /* dirichlet-logit log-PDF
   * Args:
   *   y: vector of real response values
   *   mu: vector of category logit probabilities
   *   phi: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real dirichlet_logit_lpdf(vector y, vector mu, real phi) {
     return dirichlet_lpdf(y | softmax(mu) * phi);
   }

}


data {
  int<lower=1> N;       // total number of observations
  int<lower=2> C;       // number of elements in the compositional outcome
  int<lower=2> K;       // number of elements in the exposure mixture
  vector[C] Y[N];       // compositional outcome matrix
  matrix[N, K] X_matr;  // exposure matrix
}

parameters {
  vector[C-1] rel_coefs;    // mixture index coefficients
  vector[C-1] intercepts; // mixture index intercepts
  real<lower=0.1> phi;  // precision parameter
  simplex[K] w;             // weights for the X mixture
  vector<lower = 0.0001>[K] Dir_alpha; // alphas parameters of the weights
}

transformed parameters {
  vector[N] S_bwqs = X_matr * w;
}

model {
  // Priors
  phi ~ gamma(0.01,0.01);
  intercepts ~ normal(0,10);
  rel_coefs ~ normal(0,10);
  Dir_alpha ~ gamma(2, 2);    // each alpha inside Dir_alpha is ~ gamma(2, 2)
  w ~ dirichlet(Dir_alpha);

  // likelihood
 {
    vector[C] mu[N];
    for (n in 1:N) {
      mu[n][1] = 0;  // Reference category
      mu[n][2:C] = intercepts + rel_coefs .* S_bwqs[n] ;
      target += dirichlet_logit_lpdf(Y[n] | mu[n], phi);
    }
  }
}

// add block to calculate effect estimate of reference proportion for unit increase in exposure
generated quantities {
    vector<lower=0.0001>[N] num_i;
    vector<lower=0.0001>[N] denom_i;
    vector<lower=0.0001>[N] ref_prop_i;
    vector<lower = 0.0001>[C] abs_coefs;

    for(n in 1:N){
      num_i[n] = 1.0;
      denom_i[n] = 1.0;
      for(c in 1:(C-1)){
        num_i[n] += exp(intercepts[c] + rel_coefs[c]*S_bwqs[n]);
        denom_i[n] += exp(intercepts[c] + rel_coefs[c]* (S_bwqs[n]+1) );
      }
      ref_prop_i[n] = num_i[n] / denom_i[n];
    }
    abs_coefs[1] = mean(ref_prop_i);
    abs_coefs[2:C] = exp(rel_coefs) * abs_coefs[1];
}
"
