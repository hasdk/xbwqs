#' Dirichlet-Bayesian Weighted Quantile Sum regression models.
#'
#' Fits the Dirichlet Bayesian Weighted Quantile Sum (DBWQS) regression for compositional outcomes. The function \code{dbwqs} uses the package \code{rstan} which allows the connection with the STAN software written in C++ for bayesian inference; for further information, please consult https://mc-stan.org/.
#'
#'
#' @param data Data frame containing outcomes, exposures and covariates.
#' @param out Character vector of compositional outcome names corresponding to those in `data`.
#' @param exp Character vector of exposure names included in the mixture, corresponding to those in `data`.
#' @param cov Character vector covariate names, corresponding to those in `data`. Note that the categorical covariates need to be included as factors in `data`; otherwise, they will be treated as numeric in the model.
#' @param ref Character specifying the name of the outcome variable to be used as reference for the multinomial logit link. By default (`ref=NULL`), the first outcome variable in `out` is taken to be the reference.
#' @param q Integer specifying how exposure variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). Default is \code{q=4}.
#' @param chains Integer specifying the number of chains in Hamiltonian Monte Carlo algorithm.
#' Default value \code{chains = 1}.
#' @param iter Integer specifying the number of iterations in HMC algorithm. Default value is \code{iter = 2000}. This value is chosen to make testing faster; it is recommended to use 10000 iterations in a full analysis.
#' @param warmup Integer spcifying the number of warmup iterations in HMC algorithm. Default value is \code{warmup = 1000}.
#' @param thin Integer specifying the thinning parameter in Hamiltonian Monte Carlo algorithm. Default value is \code{thin=1}.
#' @param impute_zeros Boolean value indicating if zero values in the compositional outcomes should be imputed. Default is \code{impute_zeros=FALSE}.
#' @param eps Numeric value defining a small threshold below which the outcome proportion will be rounded to zero. Default value is \code{eps=0.0001}. This parameter is only used when \code{impute_zeros=TRUE}.
#' @param seed Integer value to fix the seed for reproducibility of results.
#' @param verbose Boolean value indicating if user wants to display additional information during fitting. Default value is \code{verbose=FALSE}.
#' @param max_treedepth Integer specifying maximum tree depth for HMC algorithm. Default value is \code{max_treedepth=20}.
#' @param adapt_delta Numeric value specifying the target average acceptance probability for HMC algorithm. Default value is \code{adapt_delta = 0.999999}.
#'
#' @return
#' \code{dbwqs} returns a list containing the following elements:
#' \item{dbwqs_fit}{An \code{S4} object with all details of the Hamiltonian Monte Carlo, all the extractions
#' from the posterior distribution and all values of the parameters}
#' \item{absprop}{Table with the statistics of the mixture index coefficients in absolute proportion scale: mean, standard error of the mean,
#' standard deviation, lower and upper values for the 80% and 95% credible interval, n_eff and Rhat.}
#' \item{relprop}{Table with the statistics of the mixture index coefficients in relative proportion scale: mean, standard error of the mean,
#' standard deviation, lower and upper values for the 80% and 95% credible interval, n_eff and Rhat.}
#' \item{mixprop}{Table with the statistics of the mixture weights: mean, standard error of the mean,
#' standard deviation, lower and upper values for the 80% and 95% credible interval, n_eff and Rhat.}
#' \item{phi}{Table with the statistics of the precision parameter: mean, standard error of the mean,
#' standard deviation, lower and upper values for the 80% and 95% credible interval, n_eff and Rhat.}
#' \item{error_flag}{Flag indicating whether the fitting process was successful (0:success / -1: failed).}
#'
#' @author
#' Hachem Saddiki
#'
#' @import rstan
#' @import Rcpp
#'
#' @export
dbwqs <- function(data, out, exp, cov=NULL, ref=NULL, q=4, chains=1, impute_zeros=F, eps=1e-4,
                  seed=99221, warmup=1000, iter=2000, thin=1, verbose=F,
                  max_treedepth=20,  adapt_delta = 0.999999){
  # assert dataframe
  data <- as.data.frame(data)

  # quantize exposures
  exp_q = data.frame(quantile_split(data, mix_name = exp, q=q))

  # set outcome reference for multinomial logit link
  if(!is.null(ref)){
    ref_idx <- which(out == ref)
    out <- c(out[ref_idx], out[-ref_idx])
  }

  if(impute_zeros){
    # impute rounded zeros in outcome proportions
    out_df <- impute_zeros(data[,out], eps=eps)
  } else{
    out_df <- data[,out]
  }

  # ensure cell proportions sum to 1 exactly
  row_totals <- rowSums(out_df)
  for(o in 1:ncol(out_df)){out_df[,o] <- out_df[,o] / row_totals}

  # set DBWQS appropriate model
  if(length(cov) == 0){
    # import stan model with no covariates
    dbwqs_stanmod <- dbwqs_nocovmod

    # prepare input data list
    data_stan = list(N = nrow(data), # sample size
                     C = length(out), # Number of outcome categories
                     K = length(exp), # number of exposures in the mixture
                     Y = out_df,  # outcome matrix
                     X_matr = exp_q) # exposures quartile matrix
  } else{
    # import stan model with covariates
    dbwqs_stanmod <- dbwqs_covmod

    # check for factor variables
    fact_idx <- which(sapply(data[, cov], is.factor) == TRUE)

    # construct new covariate data frame with binary dummy variables representing the factor variables
    if(length(fact_idx) >0){
      cov_df <- fastDummies::dummy_cols(data[, cov], cov[fact_idx],
                           remove_most_frequent_dummy = TRUE,
                           remove_selected_columns = TRUE)
    } else{
      cov_df <- data[,cov]
    }

    # prepare input data list
    data_stan = list(N = nrow(data), # sample size
                     C = length(out), # Number of outcome categories
                     K = length(exp), # number of exposures in the mixture
                     M = ncol(cov_df), # number of covariates
                     Y = out_df,  # outcome matrix
                     X_matr = exp_q, # exposures quartile matrix
                     Z_matr = cov_df) # adj. covariates matrix
  }


  # run STAN with NUTS
  if(verbose){print('Starting DBWQS model fitting ...')}
  rstan_options(auto_write = TRUE)

  dbwqs_fit = suppressWarnings(rstan::stan(model_name="dbwqs_model",
                   model_code = dbwqs_stanmod,
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
  if(is.null(summary(dbwqs_fit))){
    return(list(error_flag=-1))
  } else{
    if(verbose){
      print('DBWQS model fitting finished.')
    }
    # extract effect estimates in relative proportion (reference = first outcome category in the list)
    relprop <- data.frame(summary(dbwqs_fit, pars=c('rel_coefs'), probs=c(0.025,0.975,0.1,0.9))$summary)
    colnames(relprop)[4:7] <- c('CrI95_low', 'CrI95_high','CrI80_low', 'CrI80_high')
    rownames(relprop) <- out[-1]

    # extract effect estimates in absolute proportion (%)
    absprop <- data.frame(summary(dbwqs_fit, pars=c('abs_coefs'),
                                      probs=c(0.025,0.975,0.1,0.9))$summary)
    colnames(absprop)[4:7] <- c('CrI95_low', 'CrI95_high','CrI80_low', 'CrI80_high')
    rownames(absprop) <- out
    absprop[,c(1,4:7)] = apply(absprop[,c(1,4:7)], 2 , function(x) (x-1)*100)

    # extract estimated exposure mixture weights
    mixweights <- data.frame(summary(dbwqs_fit, pars=c('w'), probs=c(0.025,0.975,0.1,0.9))$summary)
    colnames(mixweights)[4:7] <- c('CrI95_low', 'CrI95_high','CrI80_low', 'CrI80_high')
    rownames(mixweights) <- exp

    # extract estimated precision parameter phi
    phi <- data.frame(summary(dbwqs_fit, pars=c('phi'), probs=c(0.025,0.975,0.1,0.9))$summary)
    colnames(phi)[4:7] <- c('CrI95_low', 'CrI95_high','CrI80_low', 'CrI80_high')

    return(list(dbwqs_fit = dbwqs_fit, absprop=absprop, relprop=relprop,
                  mixweights=mixweights, phi_res=phi, error_flag=0))
  }
}


