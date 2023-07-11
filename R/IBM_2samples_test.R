#' Equality test of unknown component distributions in two admixture models with IBM approach
#'
#' Two-sample test of the unknown component distribution in admixture models using Inversion - Best Matching
#' (IBM) method. Recall that we have two admixture models with respective probability density functions (pdf)
#' l1 = p1 f1 + (1-p1) g1 and l2 = p2 f2 + (1-p2) g2, where g1 and g2 are known pdf and l1 and l2 are observed.
#' Perform the following hypothesis test: H0 : f1 = f2 versus H1 : f1 differs from f2.
#'
#' @param samples A list of the two observed samples, where each sample follows the mixture distribution given by l = p*f + (1-p)*g,
#'                with f and p unknown and g known.
#' @param known.p (default to NULL) Numeric vector with two elements, the known (true) mixture weights.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param sim_U Random draws of the inner convergence part of the contrast as defined in the IBM approach (see 'Details' below).
#' @param n_sim_tab Number of simulated gaussian processes used in the tabulation of the inner convergence distribution in the IBM approach.
#' @param min_size (default to NULL, only used with 'ICV' testing method in the k-sample case, otherwise useless) Minimal size among all samples (needed
#'                  to take into account the correction factor for the variance-covariance assessment).
#' @param conf.level The confidence level of the 2-samples test, i.e. the quantile level to which the test statistic is compared.
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when parallelizing.
#'
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.science/hal-03201760
#'
#' @return A list of five elements, containing : 1) the test statistic value; 2) the rejection decision; 3) the p-value of the
#'         test, 4) the estimated weights of the unknown component for each of the two admixture models, 5) the simulated distribution
#'         of the inner convergence regime (useful to perform the test when comparing to the extreme quantile of this distribution).
#'
#' @examples
#' \donttest{
#' ####### Under the null hypothesis H0 :
#' ## Simulate data:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm")
#' list.param <- list(f1 = list(mean = 1, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 1, sd = 1), g2 = list(mean = 3, sd = 1.2))
#' X.sim <- rsimmix(n= 1100, unknownComp_weight=0.85, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' Y.sim <- rsimmix(n= 1200, unknownComp_weight=0.75, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 3, sd = 1.2))
#' IBM_2samples_test(samples = list(X.sim, Y.sim), known.p = NULL, comp.dist = list.comp,
#'                   comp.param = list.param, sim_U = NULL, n_sim_tab = 6, min_size = NULL,
#'                   conf.level = 0.95, parallel = FALSE, n_cpu = 2)
#'
#' ####### Under the alternative H1 :
#' ## Simulate data:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm")
#' list.param <- list(f1 = list(mean = 1, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 2, sd = 1), g2 = list(mean = 3, sd = 1.2))
#' X.sim <- rsimmix(n= 1100, unknownComp_weight=0.85, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' Y.sim <- rsimmix(n= 1200, unknownComp_weight=0.75, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 3, sd = 1.2))
#' IBM_2samples_test(samples = list(X.sim, Y.sim), known.p = NULL, comp.dist = list.comp,
#'                   comp.param = list.param, sim_U = NULL, n_sim_tab = 6, min_size = NULL,
#'                   conf.level = 0.95, parallel = FALSE, n_cpu = 2)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_2samples_test <- function(samples, known.p = NULL, comp.dist = NULL, comp.param = NULL, sim_U = NULL,
                              n_sim_tab = 50, min_size = NULL, conf.level = 0.95, parallel = FALSE, n_cpu = 2)
{
  ## Estimate the proportions of the mixtures:
  estim <- try(suppressWarnings(IBM_estimProp(sample1 = samples[[1]], sample2 = samples[[2]], known.prop = known.p, comp.dist = comp.dist,
                                              comp.param = comp.param, with.correction = FALSE, n.integ = 1000)),
               silent = TRUE)
  count_error <- 0
  while ((inherits(x = estim, what = "try-error", which = FALSE)) & (count_error < 3)) {
    estim <- NULL
    estim <- try(suppressWarnings(IBM_estimProp(sample1 = samples[[1]], sample2 = samples[[2]], known.prop = known.p, comp.dist = comp.dist,
                                                comp.param = comp.param, with.correction = FALSE, n.integ = 1000)),
                 silent = TRUE)
    count_error <- count_error + 1
  }
  if (inherits(x = estim, what = "try-error", which = FALSE)) {
    estim.weights <- 100
    contrast_val <- NA
  } else {
    estim.weights <- estim[["prop.estim"]]
    if (is.null(min_size)) {
      sample.size <- min(length(samples[[1]]), length(samples[[2]]))
    } else {
      sample.size <- min_size
    }
    contrast_val <- sample.size * IBM_empirical_contrast(par = estim.weights, fixed.p.X = estim[["p.X.fixed"]], sample1 = samples[[1]],
                                                         sample2 = samples[[2]], G = estim[["integ.supp"]], comp.dist = comp.dist, comp.param = comp.param)
    ## Identify boolean 'green light' criterion to known whether we need to perform the test with stochastic integral tabulation:
    # green_light <- IBM_greenLight_criterion(estim.obj = estim, sample1 = samples[[1]], sample2 = samples[[2]], comp.dist = comp.dist,
    #                                        comp.param = comp.param, min_size = min_size, alpha = 0.05)
    # if (green_light) {
    #   if (is.null(sim_U)) {
    #     tab_dist <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 50, sample1 = samples[[1]], sample2 = samples[[2]], min_size = min_size,
    #                                          comp.dist = comp.dist, comp.param = comp.param, parallel = parallel, n_cpu = n_cpu)
    #     sim_U <- tab_dist$U_sim
    #     contrast_val <- tab_dist$contrast_value
    #   }
    #   reject.decision <- contrast_val > quantile(sim_U, 0.95)
    #   CDF_U <- ecdf(sim_U)
    #   p_value <- 1 - CDF_U(contrast_val)
    # } else {
    #   reject.decision <- TRUE
    #   p_value <- 1e-16
    # }
  }

  ## Earn computation time using this soft version of the green light criterion:
  if (any(abs(estim.weights) > 1)) {
    reject.decision <- TRUE
    p_value <- 1e-16
    sim_U <- NA
  } else {
    if (is.null(sim_U)) {
      tab_dist <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 60, sample1 = samples[[1]], sample2 = samples[[2]], min_size = min_size,
                                            comp.dist = comp.dist, comp.param = comp.param, parallel = parallel, n_cpu = n_cpu)
      sim_U <- tab_dist$U_sim
      contrast_val <- tab_dist$contrast_value
    }
    reject.decision <- contrast_val > stats::quantile(sim_U, conf.level)
    CDF_U <- stats::ecdf(sim_U)
    p_value <- 1 - CDF_U(contrast_val)
  }

  return( list(confidence_level = conf.level, rejection_rule = reject.decision, p_value = p_value,
               test.stat = contrast_val, weights = estim.weights, sim_U = sim_U) )
}

