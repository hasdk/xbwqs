#' Impute zero proportions in compositional outcomes with an estimated lower limit of detection
#'
#' @param df Data frame containing compositional outcomes, one for each column.
#' @param eps Numeric value defining a small threshold below which the outcome proportion will be rounded to zero (`default=1e-4`).
#'
#' @return Data frame containing compositional outcomes where all zeros are imputed with an estimated lower limit of detection.
#'
#' @import robCompositions
#' @export
impute_zeros <- function(df, eps=0.0001){
  dl = rep(0, ncol(df))
  for(i in 1:ncol(df)){
    idx = which(df[,i] < eps)
    df[idx,i] = 0
    if(length(idx) > 0){
      dl[i] = min(df[which(df[,i] >= eps),i])
    }
  }
  df_imp <- robCompositions::imputeBDLs(as.matrix(df),
                    dl= matrix(dl,ncol=ncol(df),nrow=1),
                    maxit=50,eps=0.1,R=50,method="lmrob", variation=F)
  return(df_imp$x)
}

#' Quantize one or multiple continuous variables using ecdf.
#'
#' @param data Data frame with at least one continuous variable.
#' @param mix_name Character vector containing names of variable to quantize. If null, all variables in `data` are quantized.
#' @param q Integer defining the quantile split (default = 4).
#'
#' @return Data frame with original variables replaced by their quantized version.
#' @export
quantile_split <- function(data, mix_name = NULL, q=4){
  if( (q>0) & (q<=100)) {
    if(is.null(mix_name)){
      idx <- 1:ncol(data)
    } else{
      data=data.frame(data)
      idx <- which(colnames(data)%in% mix_name)
    }
    out <- apply(data[,idx], 2, function(i){stats::ecdf(i)(i)*q})
    return(out)
  } else{return(data)}

}

#' Check that the installed rstan can compile models with the installed StanHeaders.
#'
#' rstan versions before 2.39 force the C++14 standard, while StanHeaders 2.33
#' and later require C++17. With that combination every Stan model fails to
#' compile with an opaque C++ error, so stop early with an actionable message.
#'
#' @noRd
check_rstan_compat <- function(rstan_ver = utils::packageVersion("rstan"),
                               headers_ver = utils::packageVersion("StanHeaders")){
  rstan_ver <- package_version(rstan_ver)
  headers_ver <- package_version(headers_ver)
  if(rstan_ver < "2.39.0" && headers_ver >= "2.33.0"){
    stop(paste0(
      "Incompatible Stan installation: rstan ", rstan_ver, " cannot compile models with StanHeaders ",
      headers_ver, ".\n",
      "Install a newer rstan from the Stan R-universe, restart R, and try again:\n",
      "  install.packages(\"rstan\", repos = c(\"https://stan-dev.r-universe.dev\", getOption(\"repos\")))"),
      call. = FALSE)
  }
  invisible(TRUE)
}



