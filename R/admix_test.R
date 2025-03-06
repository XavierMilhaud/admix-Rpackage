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
#'                    only 'Poly' is available and the test is a gaussianity test. For further details, see section 'Details' below.
#' @param conf_level The confidence level of the K-sample test.
#' @param ... Depending on the choice made by the user for the test method ('poly' or 'icv'), optional arguments to
#'            \link[admix]{gaussianity_test} or \link[admix]{orthobasis_test} (in case of 'poly'), and to \link[admix]{IBM_k_samples_test}
#'            in case of 'icv'.
#'            .
#'
#' @details For further details on implemented hypothesis tests, see the references hereafter.
#'          .
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024b}{admix}
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2022}{admix}
#' \insertRef{PommeretVandekerkhove2019}{admix}
#'
#' @return An object of class \link[admix]{admix_test}, containing 8 attributes: 1) the test decision (reject the null hypothesis or not);
#'         2) the p-value of the test; 3) the confidence level of the test (1-alpha, where alpha denotes the level of the test
#'         or equivalently the type-I error); 4) the value of the test statistic; 5) the number of samples under study; 6) the
#'         respective size of each sample; 7) the information about mixture components (distributions and parameters); 8) the
#'         chosen testing method (either based on polynomial basis expansions, or on the inner convergence property; see given
#'         references).
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
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'            conf_level = 0.95, test_method = "poly", ask_poly_param = FALSE, support = "Real")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_test <- function(samples, admixMod, test_method = c("poly","icv"), conf_level = 0.95, ...) {

  meth <- match.arg(test_method)

  n_samples <- length(samples)
  ## Check right specification of arguments:
  if ((n_samples == 1) & (meth == "icv")) stop("Testing using the inner convergence property (obtained from IBM estimation) requires at least TWO samples.\n")
  if (meth == "poly") message("  Testing using polynomial basis expansions requires in theory a square-root n consistent estimation
  of the proportions of the unknown component distributions (thus using 'BVdk' estimation by default, associated to unknown
  component distributions with symmetric densities). However, it is allowed to use 'PS' estimation in practice (argument 'est_method'
  has therefore to be set to 'PS'. In this case, the variance of estimators is obtained by boostrapping.\n")

  if (meth == "icv") {
    options(warn = -1)
    if (n_samples >= 2) {
      test_res <- IBM_k_samples_test(samples = samples, admixMod = admixMod, conf_level = conf_level, ...)
    } else stop("Incorrect number of samples under study (should be > 1).")

  } else if (meth == "poly") {
    if (n_samples == 1) {
      message("In the one-sample case, testing using polynomial basis expansions corresponds to a gaussianity test.\n")
      test_res <- gaussianity_test(samples = samples[[1]], admixMod = admixMod[[1]], conf_level = conf_level, ...)
    } else if (n_samples == 2) {
      test_res <- orthobasis_test(samples = samples, admixMod = admixMod, conf_level = conf_level, ...)
    } else {
      stop("Testing using polynomial basis expansions is limited to ONE-sample or TWO-samples tests!")
    }

  } else stop("Please choose appropriately the arguments of the function.")

  obj_res <- list(
    reject_decision = test_res$reject_decision,
    p_value = test_res$p_value,
    confidence_level = test_res$confidence_level,
    test_statistic_value = test_res$test_statistic_value,
    n_populations = n_samples,
    population_sizes = sapply(X = samples, FUN = length),
    admixture_models = admixMod,
    testing_meth = switch(meth, "poly" = "Polynomial expansion of the density",
                          "icv" = "Inner convergence property (following IBM)")
  )

  class(obj_res) <- "admix_test"
  obj_res$call <- match.call()

  options(warn = 0)
  return(obj_res)
}


#' Print method for objects 'admix_test'
#'
#' @param x An object of class 'admix_test'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_test <- function(x, ...)
{
  cat("Call:")
  print(x$call)
  cat("\n")
  cat("Is the null hypothesis H0 rejected? ", ifelse(x$reject_decision, "Yes", "No"), "\n", sep="")
  if (round(x$p_value, 3) == 0) {
    cat("p-value of the test: 1e-12 \n", sep="")
  } else {
    cat("p-value of the test: ", round(x$p_value, 3), "\n", sep="")
  }
}


#' Summary method for 'admix_test' objects
#'
#' Print the decision (as well as other useful information) of the statistical test with null hypothesis corresponding to
#' the equality of unknown component distributions in admixture models. More precisely, given two (or more) admixture models
#' with cumulative distribution functions (CDF) L1 and L2, where Li = pi*Fi + (1-pi)*Gi i=1,2 and Gi are the known CDFs, the
#' function performs the test: H0: F1 = F2 versus H1: F1 != F2.
#'
#' @param object An object of class 'admix_test' (see ?admix_test).
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.admix_test <- function(object, ...)
{
  cat("Call:\n")
  print(object$call)
  cat("\n------- About samples -------\n")
  cat(paste("Size of sample ", 1:object$n_populations, ": ", object$population_sizes, sep = ""), sep = "\n")
  cat("\n------ About contamination (admixture) models -----")
  cat("\n")
  if (object$n_populations == 1) {
    cat("-> Distribution and parameters of the known component \n for the admixture model: ", sep="")
    cat(object$admixture_models[[1]]$comp.dist$known, "\n")
    print(unlist(object$admixture_models[[1]]$comp.param$known, use.names = TRUE))
  } else {
    for (k in 1:object$n_populations) {
      cat("-> Distribution and parameters of the known component \n for admixture model #", k, ": ", sep="")
      cat(paste(sapply(object$admixture_models[[k]], "[[", "known")[1:2], collapse = " - "))
      cat("\n")
    }
  }
  cat("\n--------- About testing results ---------\n")
  cat("Method: ", object$testing_meth, "\n", sep = "")
  cat("Is the null hypothesis rejected? ", ifelse(object$reject_decision, "Yes", "No"), "\n", sep = "")
  cat("Type-I error is fixed to ", (1-object$confidence_level)*100, "%\n", sep = "")
  if (round(object$p_value, 3) == 0) {
    cat("p-value of the test: 1e-12 \n", sep="")
  } else {
    cat("p-value of the test: ", round(object$p_value, 3), "\n", sep="")
  }
  cat("Value of the test statistic: ", round(object$test_statistic_value,3), "\n\n", sep = "")
}
