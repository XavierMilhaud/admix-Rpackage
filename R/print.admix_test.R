#' Print the results of statistical test for equality of unknown component distributions in admixture models
#'
#' Print the decision (as well as other useful information) of the statistical test with null hypothesis corresponding to
#' the equality of unknown component distributions in admixture models. More precisely, given two (or more) admixture models
#' with cumulative distribution functions (CDF) L1 and L2, where Li = pi*Fi + (1-pi)*Gi i=1,2 and Gi are the known CDFs, the
#' function performs the test: H0: F1 = F2 versus H1: F1 != F2.
#'
#' @param x An object of class 'admix_test' (see ?admix_test).
#' @param ... further arguments passed to or from other methods.
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
#' print(gaussTest)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_test <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nIs the null hypothesis rejected? ", x$Reject_decision, "\n", sep = "")
  cat("The type-I error is fixed to ", (1-x$Confidence_level)*100, "%\n", sep = "")
  cat("The p-value of the test equals ", x$P_value, "\n", sep = "")
  cat("The value of the test statistics is ", x$Statistic_value, "\n\n", sep = "")
}
