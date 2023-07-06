#' Estimate the unknown parameters of the admixture model(s) under study
#'
#' Estimate the component weights, the location shift parameter (in case of a symmetric unknown component density),
#' and the unknown component distribution using different estimation techniques. We remind that the i-th admixture
#' model has probability density function (pdf) l_i such that:
#'    l_i = p_i * f_i + (1-p_i) * g_i, where g_i is the known component density.
#' The unknown quantities p_i and f_i then have to be estimated.
#'
#' @param samples A list of the K samples to be studied, all following admixture distributions.
#' @param sym.f A boolean indicating whether the unknown component densities are assumed to be symmetric or not.
#' @param est.method The estimation method to be applied. Can be one of 'BVdk' (Bordes and Vandekerkhove estimator), 'PS' (Patra and Sen
#'         estimator), or 'IBM' (Inversion Best-Matching approach). The same estimation method is performed on each sample.
#'         Important note: estimation by 'IBM' is unbiased only under H0, meaning that choosing this method requires to perform
#'         previously the test hypothesis between the pairs of samples. For further details, see section 'Details' below.
#' @param comp.dist A list with 2*K elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the K admixture models. Elements, grouped by 2, refer to the unknown and known components of each admixture model,
#'                  If there are unknown elements, they must be specified as 'NULL' objects. For instance, 'comp.dist' could be specified
#'                  as follows with K = 3: list(f1 = NULL, g1 = 'norm', f2 = NULL, g2 = 'norm', f3 = NULL, g3 = 'norm').
#' @param comp.param A list with 2*K elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Elements, grouped by 2, refer to the parameters of unknown and known components of each admixture model.
#'                   If there are unknown elements, they must be specified as 'NULL' objects. For instance, 'comp.param' could
#'                   be specified as follows (with K = 3):
#'                   list(f1 = NULL, g1 = list(mean=0,sd=1), f2 = NULL, g2 = list(mean=3,sd=1.1), f3 = NULL, g3 = list(mean=-2,sd=0.6)).
#'
#' @details For further details on the different estimation techniques, see i) IBM approach at https://hal.science/hal-03201760 ;
#'          ii) Patra and Sen estimator: Patra, R.K. and Sen, B. (2016); Estimation of a Two-component Mixture Model with Applications to
#'          Multiple Testing; JRSS Series B, 78, pp. 869--893. ; iii) BVdk estimator: Bordes, L. and Vandekerkhove, P. (2010);
#'          Semiparametric two-component mixture model when a component is known: an asymptotically normal estimator; Math. Meth. Stat.; 19, pp. 22--41.
#'
#' @return A list containing the estimated weight of every unknown component distribution among admixture samples.
#'
#' @examples
#' ##### On a simulated example to see whether the true parameters are well estimated.
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm")
#' list.param <- list(f1 = list(mean = 0, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 0, sd = 1), g2 = list(mean = -3, sd = 1.1))
#' ## Simulate data:
#' sim1 <- rsimmix(n = 2100, unknownComp_weight = 0.8, comp.dist = list(list.comp$f1,list.comp$g1),
#'                 comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' sim2 <- rsimmix(n= 2000, unknownComp_weight = 0.85, comp.dist = list(list.comp$f2,list.comp$g2),
#'                 comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' ## Estimate the mixture weights of the admixture models:
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = -3, sd = 1.1))
#' estim <- admix_estim(samples = list(sim1,sim2), sym.f = TRUE, est.method = 'IBM',
#'                      comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_estim <- function(samples = NULL, sym.f = FALSE, est.method = c("PS","BVdk","IBM"), comp.dist = NULL, comp.param = NULL)
{
  stopifnot( length(comp.dist) == (2*length(samples)) )
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    if ( (!all(sapply(comp.param, is.null)[seq.int(from = 2, to = length(comp.dist), by = 2)] == FALSE)) |
         (!all(sapply(comp.param, is.null)[seq.int(from = 1, to = length(comp.dist), by = 2)] == TRUE)) ) {
      stop("Component distributions and/or parameters must have been badly specified in the admixture models.")
    }
  }

  n_samples <- length(samples)
  meth <- match.arg(est.method)
  ## Check right specification of arguments:
  if ((n_samples == 1) & (meth == "IBM")) stop("Estimation by IBM technique requires at least two samples.")
  if ((sym.f == FALSE) & (meth == "BVdk")) stop("Estimation by BVdk estimator requires the unknown component density to be assumed symmetric.")

  n_obs <- sapply(X = samples, FUN = length)
  estimate <- vector(mode = "list", length = n_samples)
  if (meth == "BVdk") {
    for (k in 1:n_samples) {
      estimate[[k]] <- BVdk_estimParam(data = samples[[k]], method = 'L-BFGS-B', list(comp.dist[[2*k-1]],comp.dist[[2*k]]), list(comp.param[[2*k-1]],comp.param[[2*k]]))
    }
    estim_weight <- sapply(X = estimate, "[[", 1)
    estim_loc <- sapply(X = estimate, "[[", 2)
  } else if (meth == "PS") {
    data_transfo <- vector(mode = "list", length = n_samples)
    for (k in 1:n_samples) {
      data_transfo[[k]] <- knownComp_to_uniform(data = samples[[k]], list(comp.dist[[2*k-1]],comp.dist[[2*k]]), list(comp.param[[2*k-1]],comp.param[[2*k]]))
      estimate[[k]] <- PatraSen_est_mix_model(data = data_transfo[[k]], method = 'fixed', c.n = 0.1*log(log(n_obs[k])), gridsize = 2000)
    }
    estim_weight <- sapply(X = estimate, "[[", "alp.hat")
    estim_loc <- NA
  } else if (meth == "IBM") {
    warning("Do not forget that estimators of proportions are reliable only if unknown component distributions are tested equal!")
    for (k in 2:n_samples) {
      estimate[[k]] <- IBM_estimProp(sample1 = samples[[1]], sample2 = samples[[k]], known.prop = NULL,
                                     comp.dist = list(comp.dist[[1]],comp.dist[[2]],comp.dist[[2*k-1]],comp.dist[[2*k]]),
                                     comp.param = list(comp.param[[1]],comp.param[[2]],comp.param[[2*k-1]],comp.param[[2*k]]),
                                     with.correction = F, n.integ = 1000)
    }
    estim_weight <- c(estimate[[k]]$prop.estim[1], unlist(sapply(X = lapply(X = estimate, "[[", "prop.estim"), FUN = "[[", 2)))
    estim_loc <- NA
  } else stop("Please choose appropriately the arguments of the function.")

  estimators <- list(Estimated_weight = estim_weight,
                     Estimated_location = estim_loc,
                     Estimation_method = meth,
                     Unknown_density_symmetry_assumption = sym.f,
                     n_populations = n_samples)
  class(estimators) <- "admix_estim"
  estimators$call <- match.call()

  return(estimators)
}
