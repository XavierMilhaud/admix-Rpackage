#' Equality test for the unknown components of admixture models
#'
#' Perform hypothesis test between unknown components of a list of admixture models, where we remind that the i-th admixture
#' model has probability density function (pdf) l_i such that:
#'    l_i = p_i * f_i + (1-p_i) * g_i, with g_i the known component density.
#' The unknown quantities p_i and f_i are thus estimated, leading to the test given by the following null and alternative hypothesis:
#' H0: f_i = f_j for all i != j   against H1 : there exists at least i != j such that f_i differs from f_j.
#' The test can be performed using two methods, either the comparison of coefficients obtained through polynomial basis expansions
#' of the component densities, or by the inner-convergence property obtained using the IBM approach. See 'Details' below for further information.
#'
#' @param samples A list of the K (K > 0) samples to be studied, each one assumed to follow a mixture distribution.
#' @param admixMod A list of objects of class \link[admix]{admix_model}, containing useful information about distributions and parameters
#'                 of the contamination / admixture models under study.
#' @param test_method The testing method to be applied. Can be either 'poly' (polynomial basis expansion) or 'icv' (inner
#'                    convergence from IBM). The same testing method is performed between all samples. In the one-sample case,
#'                    only 'poly' is available and the test is a gaussianity test. For further details, see section 'Details' below.
#' @param conf_level The confidence level of the K-sample test.
#' @param ... Depending on the choice made by the user for the test method ('poly' or 'icv'), optional arguments to
#'            \link[admix]{gaussianity_test} or \link[admix]{orthobasis_test} (in case of 'poly'), and to \link[admix]{IBM_k_samples_test}
#'            in case of 'icv'.
#'
#' @details For further details on implemented hypothesis tests, see the references hereafter.
#'          .
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024b}{admix}
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2022}{admix}
#' \insertRef{PommeretVandekerkhove2019}{admix}
#'
#' @return An object of class 'htest' containing the classical attributes of the latter class, as well as
#'         other attributes specific to the inherited object class. Usually, the test decision (reject the null hypothesis or not);
#'         the confidence level of the test (1-alpha, where alpha denotes the level of the test or equivalently the type-I error);
#'         the number of samples under study; the respective size of each sample; the information about known mixture components.
#'
#' @examples
#' ####### Example with 2 samples
#' mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 350, weight = 0.85,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = -1, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' data2 <- get_mixture_data(mixt2)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'            conf_level = 0.95, test_method = "poly", ask_poly_param = FALSE, support = "Real")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_test <- function(samples, admixMod, test_method = c("poly","icv"), conf_level = 0.95, ...)
{
  if (!is.list(samples))
    stop("Please provide sample(s) and admixture model(s) in lists, also with only one sample!")
  if (!all(sapply(X = admixMod, FUN = inherits, what = "admix_model")))
    stop("Argument 'admixMod' is not correctly specified. See ?admix_model.")

  meth <- match.arg(test_method)
  n_samples <- length(samples)

  ## Check right specification of arguments:
  if ((n_samples > 2) & (meth == "poly")) stop("Testing using polynomial basis expansions involves at most TWO samples.\n")
  if ((n_samples == 1) & (meth == "icv")) stop("Testing using the inner convergence property (obtained from IBM estimation) requires at least TWO samples.\n")
  if (meth == "poly") message("  Testing using polynomial basis expansions requires in theory a square-root n consistent estimation
  of the proportions of the unknown component distributions (thus using 'BVdk' estimation by default,
  associated to unknown component distributions with symmetric densities). However, it is allowed to
  use 'PS' estimation in practice (argument 'est_method' has therefore to be set to 'PS'. In this case,
  the variance of estimators is obtained by boostrapping.\n")

  old_options_warn <- base::options()$warn
  on.exit(base::options(warn = old_options_warn))
  base::options(warn = -1)

  if (meth == "icv") {
    if (n_samples >= 2) {
      test_res <- IBM_k_samples_test(samples = samples, admixMod = admixMod, conf_level = conf_level, ...)
      specific_class <- "IBM_test"
    } else stop("Incorrect number of samples under study (should be > 1).")

  } else if (meth == "poly") {
    if (n_samples == 1) {
      test_res <- gaussianity_test(sample = samples[[1]], admixMod = admixMod[[1]], conf_level = conf_level, ...)
      specific_class <- "gaussianity_test"
    } else {  # case when n_samples == 2
      test_res <- orthobasis_test(samples = samples, admixMod = admixMod, conf_level = conf_level, ...)
      specific_class <- "orthobasis_test"
    }

  } else stop("Please choose appropriately the arguments of the function.")

  class(test_res) <- c(specific_class, "htest")
  #class(test_res) <- "htest"
  test_res$call <- match.call()
  return(test_res)
}

