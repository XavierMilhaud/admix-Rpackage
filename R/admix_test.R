#' Hypothesis test between unknown components of the admixture models under study
#'
#' Perform hypothesis test between unknown components of a list of admixture models, where we remind that the i-th admixture
#' model has probability density function (pdf) l_i such that:
#'    l_i = p_i * f_i + (1-p_i) * g_i, with g_i the known component density.
#' The unknown quantities p_i and f_i are thus estimated, leading to the test given by the following null and alternative hypothesis:
#' H0: f_i = f_j for all i != j   against H1 : there exists at least i != j such that f_i differs from f_j.
#' The test can be performed using two methods, either the comparison of coefficients obtained through polynomial basis expansions
#' of the component densities, or by the inner-convergence property obtained using the IBM approach. See 'Details' below for further information.
#'
#' @param samples A list of the K samples to be studied, all following admixture distributions.
#' @param sym.f A boolean indicating whether the unknown component densities are assumed to be symmetric or not.
#' @param test.method The testing method to be applied. Can be either 'Poly' (polynomial basis expansion) or 'ICV' (inner
#'                    convergence from IBM). The same testing method is performed between all samples. In the one-sample case,
#'                    only 'Poly' is available and the test is a gaussianity test. For further details, see section 'Details' below.
#' @param sim_U (Used only with 'ICV' testing method, otherwise useless) Random draws of the inner convergence part of the contrast
#'               as defined in the IBM approach (see 'Details' below).
#' @param n_sim_tab (Used only with 'ICV' testing method, otherwise useless) Number of simulated gaussian processes used in the
#'                   tabulation of the inner convergence distribution in the IBM approach.
#' @param min_size (Potentially used with 'ICV' testing method, otherwise useless) Minimal size among all samples (needed to take
#'                  into account the correction factor for the variance-covariance assessment).
#' @param comp.dist A list with 2*K elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the K admixture models. Elements, grouped by 2, refer to the unknown and known components of each admixture model,
#'                  If there are unknown elements, they must be specified as 'NULL' objects. For instance, 'comp.dist' could be specified
#'                  as follows with K = 3: list(f1 = NULL, g1 = 'norm', f2 = NULL, g2 = 'norm', f3 = NULL, g3 = 'rnorm').
#' @param comp.param A list with 2*K elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Elements, grouped by 2, refer to the parameters of unknown and known components of each admixture model.
#'                   If there are unknown elements, they must be specified as 'NULL' objects. For instance, 'comp.param' could
#'                   be specified as follows (with K = 3):
#'                   list(f1 = NULL, g1 = list(mean=0,sd=1), f2 = NULL, g2 = list(mean=3,sd=1.1), f3 = NULL, g3 = list(mean=-2,sd=0.6)).
#' @param support (Potentially used with 'Poly' testing method, otherwise useless) The support of the observations; one of "Real",
#'                 "Integer", "Positive", or "Bounded.continuous".
#' @param conf.level The confidence level of the K-sample test.
#' @param parallel (default to FALSE) Boolean indicating whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when parallelizing.
#'
#' @details For further details on hypothesis techniques, see i) Inner convergence through IBM approach at
#'          https://hal.science/hal-03201760 ; ii) Polynomial expansions at 'False Discovery Rate model
#'          Gaussianity test' (EJS, Pommeret & Vanderkerkhove, 2017), or 'Semiparametric two-sample admixture components comparison test:
#'          the symmetric case' (JSPI, Milhaud & al., 2021).
#'
#' @return A list containing the decision of the test (reject or not), the confidence level at which the test is performed,
#'         the p-value of the test, and the value of the test statistic (following a chi2 distribution with one degree of freedom
#'         under the null).
#'
#' @examples
#' ##### On a simulated example, with 1 sample (gaussianity test):
#' list.comp <- list(f1 = "norm", g1 = "norm")
#' list.param <- list(f1 = list(mean = 0, sd = 1), g1 = list(mean = 2, sd = 0.7))
#' ## Simulate data:
#' sim1 <- rsimmix(n = 300, unknownComp_weight = 0.85, comp.dist = list(list.comp$f1,list.comp$g1),
#'                 comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' ## Perform the test hypothesis:
#' list.comp <- list(f1 = NULL, g1 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7))
#' gaussTest <- admix_test(samples = list(sim1), sym.f = TRUE, test.method = 'Poly', sim_U = NULL,
#'                         n_sim_tab = 50, min_size = NULL, comp.dist = list.comp,
#'                         comp.param = list.param, support = "Real", conf.level = 0.95,
#'                         parallel = FALSE, n_cpu = 2)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_test <- function(samples = NULL, sym.f = FALSE, test.method = c("Poly","ICV"), sim_U = NULL, n_sim_tab = 50,
                       min_size = NULL, comp.dist = NULL, comp.param = NULL, support = c("Real","Integer","Positive","Bounded.continuous"),
                       conf.level = 0.95, parallel = FALSE, n_cpu = 2)
{
  stopifnot( length(comp.dist) == (2*length(samples)) )
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    if ( (!all(sapply(comp.param, is.null)[seq.int(from = 2, to = length(comp.dist), by = 2)] == FALSE)) |
         (!all(sapply(comp.param, is.null)[seq.int(from = 1, to = length(comp.dist), by = 2)] == TRUE)) ) {
      stop("Component distributions and/or parameters must have been badly specified in the admixture models.")
    }
  }

  n_samples <- length(samples)
  meth <- match.arg(test.method)
  supp <- match.arg(support)
  ## Check right specification of arguments:
  if ((n_samples == 1) & (meth == "ICV")) stop("Testing through the inner convergence property obtained using IBM approach requires at least two samples.")
  if ((sym.f == FALSE) & (meth == "Poly")) stop("Testing using polynomial basis expansions requires a square-root n consistent estimation
                                                of the unknown components weights, and thus symmetric unknown components.")

  if (meth == "ICV") {

    if (n_samples == 2) {
      U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 100, sample1 = samples[[1]], sample2 = samples[[2]], min_size = NULL,
                                     comp.dist = comp.dist, comp.param = comp.param, parallel = parallel, n_cpu = n_cpu)
      test_res <- IBM_2samples_test(samples = samples, known.p = NULL, comp.dist = comp.dist, comp.param = comp.param, sim_U = U[["U_sim"]],
                                    min_size=NULL, conf.level = conf.level, parallel = parallel, n_cpu = n_cpu)
    } else if (n_samples > 2) {
      test_res <- IBM_k_samples_test(samples = samples, sim_U = NULL, n_sim_tab = n_sim_tab, min_size = NULL,
                                     comp.dist = comp.dist, comp.param = comp.param, conf.level = conf.level,
                                     tune.penalty = TRUE, parallel = parallel, n_cpu = n_cpu)
    } else stop("Incorrect number of samples under study!")

  } else if (meth == "Poly") {

    if (n_samples == 1) {
      warnings("In the one-sample case, with polynomial basis expansions, this is a gaussianity test!")
      test_res <- gaussianity_test(sample1 = samples[[1]], comp.dist = comp.dist, comp.param = comp.param, K = 3,
                                   lambda = 0.1, support = supp)
    } else {
      warnings("This is a pairwise testing among all the samples provided.")
      ## Look for all possible couples on which the test will be performed :
      model.list <- lapply(X = seq.int(from = 1, to = length(comp.dist), by = 2), FUN = seq.int, length.out = 2)
      couples.list <- NULL
      for (i in 1:(length(samples)-1)) { for (j in (i+1):length(samples)) { couples.list <- rbind(couples.list,c(i,j)) } }
      test_res <- couples.expr <- couples.param <- vector(mode = "list", length = nrow(couples.list))
      for (k in 1:nrow(couples.list)) {
        couples.expr[[k]] <- comp.dist[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
        couples.param[[k]] <- comp.param[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
        names(couples.expr[[k]]) <- names(couples.param[[k]]) <- c("f1","g1","f2","g2")
        ## Hypothesis test :
        test_res[[k]] <- orthoBasis_test_H0(samples = list(samples[[couples.list[k, ][1]]], samples[[couples.list[k, ][2]]]),
                                            known.p = NULL, comp.dist = couples.expr[[k]], comp.param = couples.param[[k]],
                                            known.coef = NULL, K = 3, nb.ssEch = 2, s = 0.49, var.explicit = TRUE,
                                            nb.echBoot = NULL, support = supp, bounds.supp = NULL, est.method = "BVdk")
      }
    }

  } else stop("Please choose appropriately the arguments of the function.")

  obj_res <- list(Reject_decision = test_res$rejection_rule,
                  Confidence_level = conf.level,
                  P_value = test_res$p_value,
                  Statistic_value = test_res$test.stat)
  class(obj_res) <- "admix_test"
  obj_res$call <- match.call()

  return(obj_res)
}
