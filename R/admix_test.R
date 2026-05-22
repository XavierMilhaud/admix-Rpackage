#' Hypothesis test for the unknown component(s) in admixture model(s)
#'
#' Perform hypothesis test on the unknown component(s) of a list of admixture model(s), where we remind that the \eqn{i}-th admixture
#' model has probability density function (pdf) \eqn{\ell_i} such that:
#' \deqn{
#'   \ell_i = p_i f_i + (1 - p_i) g_i,
#' }
#' with \eqn{g_i} the known component density, and where \eqn{\ell_i} can be estimated consistently thanks to the observations.
#' The test is made on the \eqn{f_i}'s, can be performed using two methods: either the comparison of coefficients obtained through
#' polynomial basis expansions of the component densities, or by the inner-convergence property obtained using the IBM approach.
#' See 'Details' below for further information.
#'
#' @param samples A list of the K (K > 0) samples to be studied, each one assumed to follow a mixture distribution.
#' @param admixMod A list of objects of class \link[admix]{admix_model}, with information about known distributions and known parameters.
#' @param test_method The testing method to be applied. Can be either 'poly' (polynomial basis expansion) or 'icv' (inner
#'                    convergence from IBM). The same testing method is performed between all samples. In the one-sample case,
#'                    only 'poly' is available and the test is a gaussianity test. For further details, see section 'Details' below.
#' @param conf_level The confidence level of the K-sample test.
#' @param ... Depending on the choice made by the user for the test method ('poly' or 'icv'), optional arguments to
#'            \link[admix]{gaussianity_test}, \link[admix]{orthobasis_test} (in case of 'poly'), or \link[admix]{IBM_k_samples_test}
#'            in case of 'icv'.
#'
#' @details For further details on implemented hypothesis tests, see the references hereafter. When choosing the 'icv'
#'          testing method, it is recommended to use parallel computing.
#'
#' @seealso [gaussianity_test()], [orthobasis_test()], [IBM_k_samples_test()], [get_known_component()], [get_mixing_weights()],
#'          [reject_nullHyp()], [which_rank()]
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024b}{admix}
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2022}{admix}
#' \insertRef{PommeretVandekerkhove2019}{admix}
#'
#' @return An object of class \code{gaussianity_test}, \code{orthobasis_test}, or \code{IBM_test} (that inherits from
#'         class \code{htest}), containing attributes specific to the object class (in addition to classical attributes
#'         from \code{htest}). Usually, the test decision (reject the null hypothesis or not); the confidence level of
#'         the test (1-alpha, where alpha denotes the level of the test or equivalently the type-I error); the number
#'         of samples under study; the respective size of each sample; the information about known mixture components.
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
  if (!is.list(samples) || !is.list(admixMod))
    stop("Please provide sample(s) AND admixture model(s) in a list, also with only one sample!")
  if (!all(sapply(X = admixMod, FUN = inherits, what = "admix_model")))
    stop("Argument 'admixMod' is not correctly specified. See ?admix_model.")

  meth <- match.arg(test_method)
  n_samples <- length(samples)

  ## Check right specification of arguments:
  if ((n_samples > 2) & (meth == "poly")) stop("Testing using polynomial basis expansions ('poly') involves at most TWO samples.\n")
  if ((n_samples == 1) & (meth == "icv")) stop("Testing using the Inner ConVergence property ('icv') requires at least TWO samples.\n")
  if (meth == "poly") message("  Default estimation method is 'BVdk' when testing with polynomial basis expansions (ensuring
  theoretical guarantees, but relying on symmetric unknown component densities). To consider
  other frameworks in the 2-sample case, use 'PS' estimator (setting argument 'est_method' to 'PS').")

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
  test_res$call <- match.call()
  return(test_res)
}
