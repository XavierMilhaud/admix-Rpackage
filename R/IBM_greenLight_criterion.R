#' Green-light criterion to decide whether to perform full equality test between unknown components between two admixture models
#'
#' Indicate whether there is need to perform the statistical test of equality between unknown components when comparing the unknown
#' components of two samples following admixture models. Based on the IBM approach, see more in 'Details' below.
#'
#' @param estim.obj Object obtained from the estimation of the component weights related to the proportions of the unknown component
#'                  in each of the two admixture models studied.
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param min_size (optional, NULL by default) In the k-sample case, useful to provide the minimal size among all samples
#'                 (needed to take into account the correction factor for variance-covariance assessment). Otherwise, useless.
#' @param alpha Confidence level at which the criterion is assessed (used to compute the confidence bands of the estimators
#'              of the unknown component weights).
#'
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.science/hal-03201760
#'
#' @return A boolean indicating whether it is useful or useless to tabulate the contrast distribution in order to answer
#'         the testing problem (f1 = f2).
#'
#' @examples
#' \donttest{
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=550, unknownComp_weight=0.7, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=450, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the unknown component weights in the two admixture models in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' estim <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], known.prop = NULL,
#'                        comp.dist = list.comp, comp.param = list.param,
#'                        with.correction = FALSE, n.integ = 1000)
#' IBM_greenLight_criterion(estim.obj = estim, sample1 = sample1[['mixt.data']],
#'                         sample2 = sample2[['mixt.data']], comp.dist = list.comp,
#'                         comp.param = list.param, min_size = NULL, alpha = 0.05)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export


IBM_greenLight_criterion <- function(estim.obj, sample1, sample2, comp.dist = NULL, comp.param = NULL, min_size = NULL, alpha = 0.05)
{
  if (is.null(min_size)) {
    min_sample_size <- min(length(sample1), length(sample2))
  } else {
    min_sample_size <- min_size
  }

  length.support <- length(estim.obj[["integ.supp"]])
  z <- estim.obj[["integ.supp"]][round(floor(length.support/2))]
  varCov_estim <- IBM_estimVarCov_gaussVect(x = z, y = z, estim.obj = estim.obj, fixed.p1 = estim.obj[["p.X.fixed"]], known.p = NULL,
                                            sample1 = sample1, sample2 = sample2, min_size = min_size,
                                            comp.dist = comp.dist, comp.param = comp.param)

  if (length(estim.obj[["prop.estim"]]) == 2) {

    inf_bound.p1 <- estim.obj[["prop.estim"]][1] - sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p1 <- estim.obj[["prop.estim"]][1] + sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p1 <- c(inf_bound.p1, sup_bound.p1)
    inf_bound.p2 <- estim.obj[["prop.estim"]][2] - sqrt(varCov_estim[2,2]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p2 <- estim.obj[["prop.estim"]][2] + sqrt(varCov_estim[2,2]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p2 <- c(inf_bound.p2, sup_bound.p2)
    green_light_crit <- max(conf_interval.p1[1],conf_interval.p2[1]) <= 1

  } else {

    inf_bound.p <- estim.obj[["prop.estim"]][1] - sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p <- estim.obj[["prop.estim"]][1] + sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p2 <- c(inf_bound.p, sup_bound.p)
    conf_interval.p1 <- NULL
    green_light_crit <- conf_interval.p2[1] <= 1

  }

  return( list(green_light = green_light_crit, conf_interval_p1 = conf_interval.p1, conf_interval_p2 = conf_interval.p2) )
}
