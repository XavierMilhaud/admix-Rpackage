#' Print the results of estimated parameters from K admixture models
#'
#' Print the estimated weight p of the unknown component in the admixture model under study Recall that an admixture model
#' follows the cumulative distribution function (CDF) L, where L = p*F + (1-p)*G, with g a known CDF and p and f unknown quantities.
#'
#' @param x An object of class 'admix_estim' (see ?admix_estim).
#' @param ... further arguments passed to or from other methods.
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
#' print(x = estim)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_estim <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nThe number of populations/samples under study is ", x$n_populations, ".\n", sep = "")
  cat("\n\nThe list of estimated weight(s) for the unknown component distribution(s) is :\n",
      paste("  - estimated weight of the unknown component distribution for population ", 1:x$n_populations, ": ", x$Estimated_weight, collapse="\n"))
  cat("\n\nThe list of estimated location(s) for the unknown component distribution(s) is :\n",
      paste("  - estimated location of the unknown component distribution for population ", 1:x$n_populations, ": ", x$Estimated_location, collapse="\n"))
  cat("\n\nThe chosen estimation technique is ", x$Estimation_method, ".\n", sep = "")
  cat("Was the unknown density assumed to be symmetric (not important unless BVdk estimation is performed)? ",
      x$Unknown_density_symmetry_assumption, ".\n\n", sep = "")
}
