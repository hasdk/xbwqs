#' Hierarchical Bayesian Weighted Quantile Sum regression models for cross-cohort mixtures.
#'
#' Fits the Hierarchical Bayesian Weighted Quantile Sum (HBWQS) regression for analyzing the effect of exposure mixtures across multiple cohort on continuous, binary or count outcomes. The function \code{bwqs} uses the package \code{rstan} which allows the connection with the STAN software written in C++ for bayesian inference; for further information, please consult https://mc-stan.org/.
#'
#'
#' @param data Data frame containing outcome, exposures and covariates measured in multiple cohorts. The data frame should also include a column with cohort indices for each row in `data`.
#' @param out Character vector of outcome name corresponding to the one in `data`.
#' @param exp Character vector of exposure names included in the mixture, corresponding to those in `data`.
#' @param cov Character vector covariate names, corresponding to those in `data`. Note that the categorical covariates need to be included as factors in `data`; otherwise, they will be treated as numeric in the model.
#' @param cohort Character representing the name of the cohort indices column in `data`.
#' @param family Character representing the outcome distribution used in the model. Currently, the supported outcome distributions include: `gaussian` for continuous, `binomial` for binary, and `poisson` for count outcomes. For compositional or proportion outcomes, please consult the `dbwqs()` function.
#' @param q Integer specifying how exposure variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). Default is \code{q=4}.
#' @param chains Integer specifying the number of chains in Hamiltonian Monte Carlo algorithm.
#' Default value \code{chains = 1}.
#' @param iter Integer specifying the number of iterations in HMC algorithm. Default value is \code{iter = 2000}. This value is chosen to make testing faster; it is recommended to use 10000 iterations in a full analysis.
#' @param warmup Integer spcifying the number of warmup iterations in HMC algorithm. Default value is \code{warmup = 1000}.
#' @param thin Integer specifying the thinning parameter in Hamiltonian Monte Carlo algorithm. Default value is \code{thin=1}.
#' @param seed Integer value to fix the seed for reproducibility of results.
#' @param verbose Boolean value indicating if user wants to display additional information during fitting. Default value is \code{verbose=FALSE}.
#' @param max_treedepth Integer specifying maximum tree depth for HMC algorithm. Default value is \code{max_treedepth=20}.
#' @param adapt_delta Numeric value specifying the target average acceptance probability for HMC algorithm. Default value is \code{adapt_delta = 0.999999}.
#'
#' @return
#' \code{hbwqs} returns a list containing the following elements:
#' - `hbwqs_fit`: An \code{S4} object with all details of the Hamiltonian Monte Carlo, all the extractions
#' from the posterior distribution and all values of the parameters.
#' - `cohort_slopes`: Table with the statistics of the cohort-specific regression coefficients: mean, standard deviation, lower and upper values for the 80% and 95% credible interval, effective sample size, and Rhat.
#' - `cohort_intercepts`: Table with the statistics of the cohort-specific intercepts: mean, standard error of the mean,
#' standard deviation, lower and upper values for the 80% and 95% credible interval, effective sample size, and Rhat.
#' - `mixweights`: Table with the statistics of the mixture weights: mean, standard error of the mean,
#' standard deviation, lower and upper values for the 80% and 95% credible interval, effective sample size, and Rhat.
#' - `error_flag`: Flag indicating whether the fitting process was successful (0:success / -1: failed).
#'
#' @author
#' Hachem Saddiki
#'
#' @import rstan
#' @import Rcpp
#' @import fastDummies
#'
#' @export
hbwqs <- function(data, out, exp, cov=NULL, cohort, family, q=4, chains=1,
                 seed=99221, warmup=1000, iter=2000, thin=1, verbose=F,
                 max_treedepth=20,  adapt_delta = 0.999999){
  # assert dataframe
  data <- as.data.frame(data)

  # quantize exposures
  exp_q = data.frame(quantile_split(data, mix_name = exp, q=q))

  # initialize cohort
  if(is.null(cohort)){cohort = rep(1, nrow(data))}

  # set DBWQS appropriate model
  if(length(cov) == 0){

    # prepare input data list
    data_stan = list(N = nrow(data), # sample size
                     K = length(exp), # number of exposures in the mixture
                     J = length(unique(cohort)), # number of cohorts
                     y = data[,out],  # outcome vector
                     cohort = cohort, # vector of cohort indices
                     X_matr = exp_q) # exposures quartile matrix

    # import STAN model based on outcome family
    if(family == 'gaussian'){
      hbwqs_stanmod <- hbwqs_nocovmod_contY
    } else if(family == 'binomial'){
      hbwqs_stanmod <- hbwqs_nocovmod_binY
    } else if(family == 'poisson'){
      hbwqs_stanmod <- hbwqs_nocovmod_countY
    } else{
      stop("Family type not recognized. Please input one of the following: 'gaussian', 'binomial', or 'poisson'. ")
    }
  } else{

    # check for factor variables
    fact_idx <- which(sapply(data[, cov], is.factor) == TRUE)

    # construct new covariate data frame with binary dummy variables representing the factor variables
    if(length(fact_idx) >0){
      cov_df <- fastDummies::dummy_cols(data[, cov], cov[fact_idx],
                                        remove_most_frequent_dummy = TRUE,
                                        remove_selected_columns = TRUE)
    } else{cov_df <- data[,cov]}

    # prepare input data list
    data_stan = list(N = nrow(data), # sample size
                     K = length(exp), # number of exposures in the mixture
                     M = ncol(cov_df), # number of covariates
                     J = length(unique(cohort)), # number of cohorts
                     cohort = cohort, # vector of cohort indices
                     y = data[, out],  # outcome vector
                     X_matr = exp_q, # exposures quartile matrix
                     Z_matr = cov_df) # adj. covariates matrix

    # import STAN model based on outcome family
    if(family == 'gaussian'){
      hbwqs_stanmod <- hbwqs_covmod_contY
    } else if(family == 'binomial'){
      hbwqs_stanmod <- hbwqs_covmod_binY
    } else if(family == 'poisson'){
      hbwqs_stanmod <- hbwqs_covmod_countY
    } else{
      stop("Family type not recognized. Please input one of the following: 'gaussian', 'binomial', or 'poisson'. ")
    }
  }

  # run STAN with NUTS
  if(verbose){print('Starting HBWQS model fitting ...')}
  rstan_options(auto_write = TRUE)

  hbwqs_fit = suppressWarnings(rstan::stan(model_name="hbwqs_model",
                                          model_code = hbwqs_stanmod,
                                          data = data_stan,
                                          warmup = warmup,
                                          iter = iter,
                                          thin = thin,
                                          chains = chains,
                                          seed = seed,
                                          algorithm = 'NUTS',
                                          control = list(max_treedepth=max_treedepth,
                                                         adapt_delta = adapt_delta)))

  # return error code -1 if dbwqs fit failed
  if(is.null(summary(hbwqs_fit))){
    return(list(error_flag=-1))
  } else{
    if(verbose){
      print('HBWQS model fitting finished.')
    }
    # extract effect estimates cohort-specific intercepts and slopes
    cohort_slopes <- data.frame(summary(hbwqs_fit, pars=c('b'), probs=c(0.025,0.975,0.1,0.9))$summary)
    cohort_intercepts <- data.frame(summary(hbwqs_fit, pars=c('a'), probs=c(0.025,0.975,0.1,0.9))$summary)
    colnames(cohort_slopes)[4:7] <- colnames(cohort_intercepts)[4:7] <- c('CrI95_low', 'CrI95_high','CrI80_low', 'CrI80_high')

    # extract estimated exposure mixture weights
    mixweights <- data.frame(summary(hbwqs_fit, pars=c('w'), probs=c(0.025,0.975,0.1,0.9))$summary)
    colnames(mixweights)[4:7] <- c('CrI95_low', 'CrI95_high','CrI80_low', 'CrI80_high')
    rownames(mixweights) <- exp

    return(list(hbwqs_fit = hbwqs_fit, cohort_slopes=cohort_slopes, cohort_intercepts=cohort_intercepts,
                mixweights=mixweights, error_flag=0))
  }
}


