#' DBWQS change plot.
#' Composition proportions distribution plot.
#'
#' Construct stacked bar plot representing the distribution of a given compositional variable.
#'
#' @param ccomp data frame containing the compositional variables.
#' @param sort_by Character variable representing the name of the variable to sort by in `ccomp`. Default value is `NULL` (no sorting).
#' @param decreasing `TRUE/FALSE`: should the sorting be decreasing or increasing? Default is `TRUE`. Only considered when `sort_by` is not `NULL`.
#'
#' @return ggplot object containing the composition proportions distribution plot.
#'
#' @importFrom ggplot2 ggplot geom_col aes theme_bw xlab ylab theme ggtitle
#' @import reshape2
#'
#' @export
ccomp_plot <- function(ccomp, sort_by=NULL, decreasing=TRUE){

  # assert data frame
  ccomp <- as.data.frame(ccomp)

  # sort if needed
  if(!is.null(sort_by)){
    ccomp <- ccomp[order(ccomp[,sort_by], decreasing=decreasing),]
  }

  # reshape data frame to long and add sample id
  ccomp$SampleID = 1:nrow(ccomp)
  melt.df <- reshape2::melt(ccomp,id.vars = 'SampleID',variable.name = 'Comp')

  # construct bar plot
  plt_ccomp <- ggplot2::ggplot(data=melt.df, ggplot2::aes(x = SampleID, y = value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Comp), width =2)+
    ggplot2::ylab("Composition Proportions")+
    ggplot2::xlab("Sample ID")+ggplot2::theme_bw()+
    ggplot2::theme(axis.title=element_text(size=13),
          axis.text.y=element_text(size=13),
          axis.text.x = element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(size=13),
          panel.background = element_blank()) +
    ggplot2::ggtitle(sprintf('Composition Distribution (N=%d)',nrow(ccomp)))

  return(plt_ccomp)
}

#' Auto-correlation plots from fitted DBWQS model.
#'
#' Construct auto-correlation plots using R-STAN plotting functions
#'
#' @param dbwqs_obj Object containing the fitted DBWQS model (from `dbwqs` output).
#'
#' @return Returns a list containins 4 plot objects:
#' - `relcoef_acplot`: Auto-correlation plots for coefficients associated with the exposure mixture index in relative proportion scale.
#' - `abscoef_acplot`: Auto-correlation plots for coefficients associated with the exposure mixture index in absolute proportion scale.
#' - `mixweights_acplot`: Auto-correlation plots for the estimated exposure mixture weights.
#' - `phi_acplot`: Auto-correlation plot for the precision parameter `phi`.
#'
#' @import ggplot2
#' @import rstan
#'
#' @export
dbwqs_acplot <- function(dbwqs_obj){

  # stan auto correlation plots for relative scale coefficients
  relcoef_acplot <- rstan::stan_ac(dbwqs_obj$dbwqs_fit, pars=c('rel_coefs'))
  levels(relcoef_acplot$data$parameters) <- rownames(dbwqs_obj$relprop)

  # stan auto correlation plots for absolute scale coefficients
  abscoef_acplot <- rstan::stan_ac(dbwqs_obj$dbwqs_fit, pars=c('abs_coefs'))
  levels(abscoef_acplot$data$parameters) <- rownames(dbwqs_obj$absprop)

  # stan auto correlation plots for exposure mixture weights
  mixweights_acplot <- rstan::stan_ac(dbwqs_obj$dbwqs_fit, pars=c('w'))
  levels(mixweights_acplot$data$parameters) <- rownames(dbwqs_obj$mixweights)

  # stan auto correlation plot for precision parameter phi
  phi_acplot <- rstan::stan_ac(dbwqs_obj$dbwqs_fit, pars=c('phi'))

  # pack and return list of plots
  return(list(relcoef_acplot=relcoef_acplot, abscoef_acplot=abscoef_acplot,
              mixweights_acplot=mixweights_acplot, phi_acplot=phi_acplot))
}

#' Trace plots from fitted DBWQS model.
#'
#' Construct trace plots using R-STAN plotting functions.
#'
#' @param dbwqs_obj Object containing the fitted DBWQS model (from `dbwqs` output).
#'
#' @return Returns a list containins 4 plot objects:
#'  - `relcoef_traceplot`: Trace plots for coefficients associated with the exposure mixture index in relative proportion scale.
#'  - `abscoef_traceplot`: Trace plots for coefficients associated with the exposure mixture index in absolute proportion scale.
#'  - `mixweights_traceplot`: Trace plots for the estimated exposure mixture weights.
#'  - `phi_traceplot`: Trace plot for the precision parameter `phi`.
#'
#' @import ggplot2
#' @import rstan
#'
#' @export
dbwqs_traceplot <- function(dbwqs_obj){

  # cleanup clutter in x-axis by only showing min and max ticks
  max_tick <- dbwqs_obj$dbwqs_fit@stan_args[[1]]$iter
  min_tick <- max_tick - dbwqs_obj$dbwqs_fit@stan_args[[1]]$warmup

  # stan trace plots for relative scale coefficients
  relcoef_traceplot <- rstan::stan_trace(dbwqs_obj$dbwqs_fit, pars=c('rel_coefs'))+
    ggplot2::scale_x_continuous(breaks=c(min_tick, max_tick))
  levels(relcoef_traceplot$data$parameter) <- rownames(dbwqs_obj$relprop)

  # stan trace plots for absolute scale coefficients
  abscoef_traceplot <- rstan::stan_trace(dbwqs_obj$dbwqs_fit, pars=c('abs_coefs'))+
    ggplot2::scale_x_continuous(breaks=c(min_tick, max_tick))
  levels(abscoef_traceplot$data$parameter) <- rownames(dbwqs_obj$absprop)

  # stan auto correlation plots for exposure mixture weights
  mixweights_traceplot <- rstan::stan_trace(dbwqs_obj$dbwqs_fit, pars=c('w'))+
    ggplot2::scale_x_continuous(breaks=c(min_tick, max_tick))
  levels(mixweights_traceplot$data$parameter) <- rownames(dbwqs_obj$mixweights)

  # stan auto correlation plot for precision parameter phi
  phi_traceplot <- rstan::stan_trace(dbwqs_obj$dbwqs_fit, pars=c('phi'))+
    ggplot2::scale_x_continuous(breaks=c(min_tick, max_tick))

  # pack and return list of plots
  return(list(relcoef_traceplot=relcoef_traceplot,
              abscoef_traceplot=abscoef_traceplot,
              mixweights_traceplot=mixweights_traceplot,
              phi_traceplot=phi_traceplot))
}
