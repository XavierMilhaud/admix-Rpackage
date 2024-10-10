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
#' @param samples (list) A list of the K (K>0) samples to be studied, each one assumed to follow a mixture distributions.
#' @param admixMod (list) A list of objects of class 'admix_model', containing useful information about distributions and parameters.
#' @param test_method The testing method to be applied. Can be either 'poly' (polynomial basis expansion) or 'icv' (inner
#'                    convergence from IBM). The same testing method is performed between all samples. In the one-sample case,
#'                    only 'Poly' is available and the test is a gaussianity test. For further details, see section 'Details' below.
#' @param sim_U (Only with 'icv' testing method, otherwise useless) Random draws of the inner convergence part of the contrast
#'               as defined in the IBM approach (see 'Details' below).
#' @param n_sim_tab (Only with 'icv' testing method, otherwise useless) Number of simulated gaussian processes used in the
#'                   tabulation of the inner convergence distribution in the IBM approach.
#' @param ICV_tunePenalty (Only with 'icv' testing method, otherwise useless. Default to TRUE) Boolean used to tune the penalty term in the
#'                        k-sample test (k=2,3,...,K) when using Inversion Best Matching (IBM) approach coupled to Inner ConVergence (icv)
#'                        property. Particularly useful when studying unbalanced samples (in terms of sample size) or small-sized samples.
#' @param ask_poly_param (Only with 'poly' testing method, otherwise useless. Boolean, default to FALSE) If TRUE, ask the user to choose
#'                        both the order 'K' of expansion coefficients in the orthonormal polynomial basis, and the penalization rate 's'
#'                        involved on the penalization rule for the test. Default values for these two parameters are 'K=3' and 's=0.25'.
#' @param support (Used with 'poly' testing method, otherwise useless) The support of the observations; one of "Real",
#'                 "Integer", "Positive", or "Bounded.continuous".
#' @param conf_level The confidence level of the K-sample test.
#' @param parallel (default to FALSE) Boolean indicating whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when parallelizing.
#'
#' @details For further details on hypothesis tests, see i) Inner convergence through IBM approach ; ii) Polynomial expansions.
#'          .
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2022}{admix}
#' \insertRef{PommeretVandekerkhove2019}{admix}
#'
#' @return A list containing the decision of the test (reject or not), the confidence level at which the test is performed,
#'         the p-value of the test, and the value of the test statistic (following a chi2 distribution with one degree of freedom
#'         under the null).
#'
#' @examples
#' ####### Example with 1 sample (gaussianity test):
#' mixt1 <- twoComp_mixt(n = 450, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(c("mean" = -2, "sd" = 0.5),
#'                                         c("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' admixMod <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admix_test(samples = list(data1), admixMod = list(admixMod),
#'            test_method = "poly", ask_poly_param = FALSE, support = "Real",
#'            conf_level = 0.95, parallel = FALSE, n_cpu = 2)
#'
#' ####### Example with 2 samples
#' mixt1 <- twoComp_mixt(n = 450, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(c("mean" = -2, "sd" = 0.5),
#'                                         c("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 450, weight = 0.24,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(c("mean" = -2, "sd" = 0.5),
#'                                         c("mean" = -1, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'            test_method = "poly", ask_poly_param = FALSE, support = "Real",
#'            conf_level = 0.95, parallel = FALSE, n_cpu = 2)
#'
#' ####### Example with 3 samples
#' mixt1 <- twoComp_mixt(n = 450, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(c("mean" = -2, "sd" = 0.5),
#'                                         c("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 600, weight = 0.24,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(c("mean" = -2, "sd" = 0.5),
#'                                         c("mean" = -1, "sd" = 1)))
#' mixt3 <- twoComp_mixt(n = 400, weight = 0.53,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(c("mean" = -2, "sd" = 0.5),
#'                                         c("mean" = 2, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' data3 <- getmixtData(mixt3)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
#'                          knownComp_param = mixt3$comp.param[[2]])
#'
#' admix_test(samples = list(data1, data2, data3),
#'            admixMod = list(admixMod1, admixMod2, admixMod3),
#'            test_method = "icv", n_sim_tab = 10, ICV_tunePenalty = FALSE,
#'            conf_level = 0.95, parallel = FALSE, n_cpu = 2)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_test <- function(samples, admixMod, test_method = c("poly","icv"), sim_U = NULL, n_sim_tab = 50,
                       ICV_tunePenalty = TRUE, ask_poly_param = FALSE,
                       support = c("Real","Integer","Positive","Bounded.continuous"),
                       conf_level = 0.95, parallel = FALSE, n_cpu = 2)
{
  meth <- match.arg(test_method)
  supp <- match.arg(support)

  n_samples <- length(samples)
  ## Check right specification of arguments:
  if ((n_samples == 1) & (meth == "icv")) stop("Testing through the inner convergence property obtained using IBM approach requires at least TWO samples.\n")
  if (meth == "poly") message("Testing using polynomial basis expansions should require a square-root n consistent estimation
  of the proportions of the unknown component distributions, and thus symmetric unknown densities.\n")

  if (meth == "icv") {

    if (n_samples == 2) {
      U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 100, samples = samples, admixMod = admixMod,
                                     min_size = NULL, parallel = parallel, n_cpu = n_cpu)
      test_res <- IBM_2samples_test(samples = samples, admixMod = admixMod, sim_U = U[["U_sim"]],
                                    conf_level = conf_level, parallel = parallel, n_cpu = n_cpu)
    } else if (n_samples > 2) {
      test_res <- IBM_k_samples_test(samples = samples, admixMod = admixMod, sim_U = NULL, n_sim_tab = n_sim_tab,
                                     conf_level = conf_level, tune_penalty = ICV_tunePenalty,
                                     parallel = parallel, n_cpu = n_cpu)
    } else stop("Incorrect number of samples under study (should be > 1).")

  } else if (meth == "poly") {
    if (ask_poly_param) {
      K.user <- base::readline("Please enter 'K' (integer), the order for the polynomial expansion in the orthonormal basis: ")
      s.user <- base::readline("Please enter 's' in ]0,0.5[, involved in the penalization rule for model selection where lower values of 's' lead to more powerful tests: ")
      if (n_samples > 1) {
        est_tech <- base::readline("Choose between 'BVdk' or 'PS' estimation method, given that theoretically speaking
        'BVdk' should be prefered in a testing perspective but is valid only if unknown component densities are symetric: ")
      }
    } else {
      K.user <- 3
      s.user <- 0.25
      est_tech <- "BVdk"
    }
    if (n_samples == 1) {
      message("In the one-sample case, with polynomial basis expansions, the implemented test is a gaussianity test.\n")
      test_res <- gaussianity_test(sample1 = samples[[1]], admixMod = admixMod[[1]],
                                   K = as.numeric(K.user), s = as.numeric(s.user), support = supp)
    } else {
      message("Pairwise tests among all the samples provided.\n")
      ## Look for all possible couples on which the test will be performed :
      couples.list <- NULL
      for (i in 1:(length(samples)-1)) { for (j in (i+1):length(samples)) { couples.list <- rbind(couples.list,c(i,j)) } }

      test_res <- vector(mode = "list", length = nrow(couples.list))
      for (k in 1:nrow(couples.list)) {
        test_res[[k]] <- orthobasis_test(samples = samples[as.numeric(couples.list[k, ])],
                                         admixMod = admixMod[as.numeric(couples.list[k, ])],
                                         K = as.numeric(K.user), s = as.numeric(s.user),
                                         est_method = est_tech, conf_level = conf_level, nb_echBoot = 100, support = supp)[1:7]
        test_res[[k]]$pair_tested <- paste("Sample", as.character(couples.list[k, ]), collapse = " / ")
      }
    }

  } else stop("Please choose appropriately the arguments of the function.")

  if ((meth == "poly") & (n_samples > 1)) {
    obj_res <- list(
      pairs_test_order = sapply(test_res, "[[", "pair_tested"),
      reject_decision = sapply(test_res, "[[", "reject_decision"),
      p_value = sapply(test_res, "[[", "p_value"),
      confidence_level = sapply(test_res, "[[", "confidence_level"),
      test_statistic_value = sapply(test_res, "[[", "test_statistic_value")
    )
  } else {
    obj_res <- list(
      reject_decision = test_res$reject_decision,
      p_value = test_res$p_value,
      confidence_level = test_res$confidence_level,
      test_statistic_value = test_res$test_statistic_value
    )
  }

  #class(obj_res) <- "admix_test"
  #obj_res$call <- match.call()

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
  print(x)
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
  cat("\nIs the null hypothesis rejected? ", object$reject_decision, "\n", sep = "")
  cat("The type-I error is fixed to ", (1-object$confidence_level)*100, "%\n", sep = "")
  cat("The p-value of the test equals ", object$p_value, "\n", sep = "")
  cat("The value of the test statistics is ", object$test_statistic_value, "\n\n", sep = "")
}
