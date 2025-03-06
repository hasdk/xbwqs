#' Bayesian Weighted Quantile Sum regression models.
#'
#' Fits the Bayesian Weighted Quantile Sum (BWQS) regression for continuous, binary or count outcomes. The function \code{bwqs} uses the package \code{rstan} which allows the connection with the STAN software written in C++ for bayesian inference; for further information, please consult https://mc-stan.org/.
#'
#'
#' @param data Data frame containing outcomes, exposures and covariates.
#' @param out Character vector of outcome name corresponding to the one in `data`.
#' @param exp Character vector of exposure names included in the mixture, corresponding to those in `data`.
#' @param cov Character vector covariate names, corresponding to those in `data`. Note that the categorical covariates need to be included as factors in `data`; otherwise, they will be treated as numeric in the model.
#' @param family Character representing the outcome distribution used in the model. Currently, the supported outcome distributions include: `gaussian` for continuous, `binomial` for binary, and `poisson` for count outcomes. For compositional or proportion outcomes, please consult the `dbwqs()` function.
#' @param q Integer specifying how exposure variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). Default is \code{q=4}.
#' @param chains Integer specifying the number of chains in Hamiltonian Monte Carlo algorithm.
#' Default value \code{chains = 1}.
#' @param iter Integer specifying the number of iterations in HMC algorithm. Default value is \code{iter = 2000}. This value is chosen to make testing faster; it is recommended to use 10000 iterations in a full analysis.
#' @param warmup Integer spcifying the number of warmup iterations in HMC algorithm. Default value is \code{warmup = 2000}.
#' @param thin Integer specifying the thinning parameter in Hamiltonian Monte Carlo algorithm. Default value is \code{thin=1}.
#' @param seed Integer value to fix the seed for reproducibility of results.
#' @param verbose Boolean value indicating if user wants to display additional information during fitting. Default value is \code{verbose=FALSE}.
#' @param max_treedepth Integer specifying maximum tree depth for HMC algorithm. Default value is \code{max_treedepth=20}.
#' @param adapt_delta Numeric value specifying the target average acceptance probability for HMC algorithm. Default value is \code{adapt_delta = 0.999999}.
#'
#' @return
#' \code{bwqs} returns a list containing the following elements:
#' - `bwqs_fit`: An \code{S4} object with all details of the Hamiltonian Monte Carlo, all the extractions from the posterior distribution and all values of the parameters.
#' - `coefs`: Table with the statistics of the regression coefficients: mean, standard deviation, lower and upper values for the 80% and 95% credible interval, effective sample size, and Rhat.
#' - `mixweights`: Table with the statistics of the mixture weights: mean, standard error of the mean,
#' standard deviation, lower and upper values for the 80% and 95% credible interval, effective sample size, and Rhat.
#' - `error_flag`: Flag indicating whether the fitting process was successful (0:success / -1: failed).
#'
#' @author
#' Hachem Saddiki
#'
#' @import rstan
#' @import Rcpp
#'
#' @export
bwqs <- function(data, out, exp, cov=NULL, family, q=4, chains=1,
                  seed=99221, warmup=1000, iter=2000, thin=1, verbose=F,
                  max_treedepth=20,  adapt_delta = 0.999999){
  # assert dataframe
  data <- as.data.frame(data)

  # quantize exposures
  exp_q = data.frame(quantile_split(data, mix_name = exp, q=q))

  # set DBWQS appropriate model
  if(length(cov) == 0){

    # prepare input data list
    data_stan = list(N = nrow(data), # sample size
                     K = length(exp), # number of exposures in the mixture
                     y = data[,out],  # outcome vector
                     X_matr = exp_q) # exposures quartile matrix

    # import STAN model based on outcome family
    if(family == 'gaussian'){
      bwqs_stanmod <- bwqs_nocovmod_contY
    } else if(family == 'binomial'){
      bwqs_stanmod <- bwqs_nocovmod_binY
    } else if(family == 'poisson'){
      bwqs_stanmod <- bwqs_nocovmod_countY
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
                     y = data[, out],  # outcome vector
                     X_matr = exp_q, # exposures quartile matrix
                     Z_matr = cov_df) # adj. covariates matrix

    # import STAN model based on outcome family
    if(family == 'gaussian'){
      bwqs_stanmod <- bwqs_covmod_contY
    } else if(family == 'binomial'){
      bwqs_stanmod <- bwqs_covmod_binY
    } else if(family == 'poisson'){
      bwqs_stanmod <- bwqs_covmod_countY
    } else{
      stop("Family type not recognized. Please input one of the following: 'gaussian', 'binomial', or 'poisson'. ")
    }
  }

  # run STAN with NUTS
  if(verbose){print('Starting BWQS model fitting ...')}
  rstan_options(auto_write = TRUE)

  bwqs_fit = suppressWarnings(rstan::stan(model_name="bwqs_model",
                                           model_code = bwqs_stanmod,
                                           data = data_stan,
                                           warmup = warmup,
                                           iter = iter,
                                           thin = thin,
                                           chains = chains,
                                           seed = seed,
                                           algorithm = 'NUTS',
                                           control = list(max_treedepth=max_treedepth,
                                                          adapt_delta = adapt_delta)))

  # return error code -1 if bwqs fit failed
  if(is.null(summary(bwqs_fit))){
    return(list(error_flag=-1))
  } else{
    if(verbose){
      print('BWQS model fitting finished.')
    }
    # extract effect estimates of intercept and slope
    coefs <- data.frame(summary(bwqs_fit, pars=c('intercept','coef'), probs=c(0.025,0.975,0.1,0.9))$summary)
    colnames(coefs)[4:7] <- c('CrI95_low', 'CrI95_high','CrI80_low', 'CrI80_high')

    # extract estimated exposure mixture weights
    mixweights <- data.frame(summary(bwqs_fit, pars=c('w'), probs=c(0.025,0.975,0.1,0.9))$summary)
    colnames(mixweights)[4:7] <- c('CrI95_low', 'CrI95_high','CrI80_low', 'CrI80_high')
    rownames(mixweights) <- exp

    return(list(bwqs_fit = bwqs_fit, coefs = coefs,
                mixweights = mixweights, error_flag=0))
  }
}


