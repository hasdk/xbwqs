#' Generate simulated compositional regression data set with exposure mixtures, covariates, and compositional outcomes.
#'
#' @param N Integer specifying the number of simulated observations (default=200).
#' @param C Integer specifying the number of components in the simulated compositional outcomes (default=3).
#' @param K Integer specifying the number of exposures included in the simulated exposure mixture (default=4).
#' @param J Integer specifying the number of simulated covariates (default = 2).
#' @param q Integer defining the exposure quantile split (default = 4).
#' @param rho Numeric value in (0,1) specifying the correlation between simulated exposures (default=0.5).
#' @param phi Numeric value greater than 0 specifying the true precision parameter for the DBWQS model (default=5).
#' @param thetas Numeric vector of size `C` representing the true regression coefficients of the exposure mixture index. The first element of `thetas` is equal to 0 by definition. If `NULL` (default), coefficients are randomly drawn from a uniform distribution.
#' @param betas Numeric matrix of size `J` by `C` representing the true regression coefficients of the simulated covariates. The first column of `betas` is equal to 0 by definition. If `NULL` (default), coefficients are randomly drawn from a uniform distribution.
#' @param wt Numeric vector of size `K` representing the exposure mixture weights. These weights should all be positive and sum to 1. If `NULL` (default), weights are randomly drawn from a Dirichlet distribution.
#' @param seed Integer specifying seed value for reproducibility.
#'
#' @return This function returns a list containing the following:
#' \item{obsdata}{Data frame of simulated observations of exposures, covariates, and compositional outcomes.}
#' \item{thetas}{Vector of true effect estimates associated with the exposure mixture index in relative proportion scale.}
#' \item{absprops_per}{Vector of true effect estimates associated with exposure mixture index in absolute proportion percentage scale}
#' \item{wt}{Vector of true exposure mixture weights used to simulate the exposue mixture index.}
#' \item{betas}{Matrix of true regression coefficients associated with the covariates.}
#'
#' @author
#' Hachem Saddiki
#'
#' @import LaplacesDemon
#'
#' @export
gen_dbwqs_simdata <- function(thetas=NULL, betas=NULL,wt=NULL,N=200,C=3,K=4,J=2,rho=0.5,
                        phi=5, q=4, seed=99221){

  ## 1) simulate J covariates and add column of ones for intercept
  if(J > 0){
    X = cbind(rep(1,N), LaplacesDemon::rmvn(n=N, mu=rep(0, J), Sigma=diag(J)))
  } else{
    X = as.matrix(rep(1,N))
  }

  ## 2) simulate K exposures

  # var-covar matrix
  Sig = diag(K) + matrix(rho, nrow=K, ncol=K)

  # draw K exposures from mvnorm
  Ex = LaplacesDemon::rmvn(n=N, mu=rep(0, K), Sigma=Sig)

  # calculate K exposure quantiles
  Q <- quantile_split(data=Ex,q=q)

  ## 3) calculate weighted quantile sum

  # draw weights from Dirichlet(pi)
  if(is.null(wt)){
    wt_pi = c(4,rep(2,K-1))
    wt <- as.vector(LaplacesDemon::rdirichlet(1, wt_pi))
  }
  S <- Q %*% wt

  ## 4) simulate regression coefficients

  # simulate true covariate coefficients from Uniform(-1,1), including intercept
  if(is.null(betas)){
    if(J > 0){
      betas <- t(cbind(c(0, stats::runif(C-1, -0.98, 0.99)),
                       rbind(rep(0,J), matrix(stats::runif((C-1)*J, -0.98, 0.98),
                                              nrow=C-1, ncol=J))))
    } else{
      betas <- t(as.matrix(c(0, stats::runif(C-1, -0.98, 0.99))))
    }
  }

  # simulate true exposure mixture coefficient
  if(is.null(thetas)){
    thetas <- stats::runif(C-2, -0.98,0.98)
    thetas <- ifelse(abs(thetas)<0.1, 0.1,thetas)
    thetas <- c(0,0,thetas)
  }

  ## 5) calculate Dirichlet regression mean vectors
  num <- sapply(1:C, function(i){exp(S*thetas[i] + X %*% betas[,i])})
  denom <- rowSums(num)
  mu <- t(sapply(1:N, function(i){num[i,]/denom[i]}))

  ## 6) Calculate Dirichlet regression alpha parameter
  alpha_mat <- mu * phi

  ## 7) Simulate Dirichlet regression outcome
  Y <- t(apply(alpha_mat, 1, function(i){LaplacesDemon::rdirichlet(n=1, alpha=i)}))

  ## 8) calculate true abs prop effects for 1 unit increase in S
  # calculate b_ref (eq.44)
  bref_num = rowSums(sapply(1:C, function(i){exp(S*thetas[i] + X %*% betas[,i])}))
  bref_denom = rowSums(sapply(1:C, function(i){exp((S+1)*thetas[i] + X %*% betas[,i])}))
  bref = mean(bref_num / bref_denom)
  absprops = sapply(1:C, function(i) exp(thetas[i]) * bref)
  absprops_per = (absprops-1)*100

  ## pack and return simulated data and params
  if(J==0){
    obsdata = data.frame(cbind(Y, Ex))
    colnames(obsdata) <- c(paste0('Y',1:C), paste0('Ex',1:K))
  } else{
    obsdata = data.frame(cbind(Y, Ex, X[,-1]))
    colnames(obsdata) <- c(paste0('Y',1:C), paste0('Ex',1:K), paste0('Cov',1:J))
  }

  return(list(obsdata=obsdata, betas=betas, theta=thetas, wt=wt,
              absprops_per=absprops_per))

}
