#' Create an object of class 'admix_model'
#'
#' Create an admixture model, also known as (aka) a contamination model. Such a model is a two-component
#' mixture model with one known component. Both the second component distribution and the mixing weight
#' are unknown.
#'
#' @param knownComp_dist (Character) The name of the distribution (specified as in R glossary) of the known component
#'                        of the admixture model
#' @param knownComp_param (Character) A vector of the names of the parameters (specified as in R glossary) involved in
#'                        the chosen known distribution, with their values.
#'
#' @return An object of class \link[admix]{admix_model}, containing 2 attributes: 1) a list that gives the information about the distributions
#'         involved in the two-component mixture model (the unknown and the known ones); 2) a list that gives the information about
#'         the corresponding parameters of those distributions.
#'
#' @examples
#' admix_model(knownComp_dist = "norm", knownComp_param = list("mean"=0, "sd"=1))
#' admix_model(knownComp_dist = "exp", knownComp_param = list("rate"=2))
#' admix_model(knownComp_dist = "multinom", knownComp_param = list("size"=1, "prob"=c(0.2,0.8,0.1)))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_model <- function(knownComp_dist, knownComp_param)
{
  stopifnot("Specify only ONE distribution" = length(knownComp_dist) == 1)
  #dist_nm <- base::readline(prompt = paste("Please enter the name (in R glossary) of the known component distribution: "))
  dist_table <- EnvStats::Distribution.df[ ,c("Name", "Type", "Number.parameters", "Parameter.1",
                                              "Parameter.2", "Parameter.3", "Parameter.4", "Parameter.5")]
  stopifnot("Unknown specified distribution" = knownComp_dist %in% rownames(dist_table) | knownComp_dist == "multinom" | knownComp_dist == "gompertz")

  if (knownComp_dist == "multinom" | knownComp_dist == "gompertz") { nparam_theo <- 2
  } else { nparam_theo <- dist_table[rownames(dist_table) == knownComp_dist, "Number.parameters"] }
  stopifnot("Mispecification of parameters" = nparam_theo == length(knownComp_param))

  if (knownComp_dist == "multinom") {
    stopifnot("Name of parameters not appropriate" = all(names(knownComp_param) == c("size","prob")))
  } else if (knownComp_dist == "gompertz") {
    stopifnot("Name of parameters not appropriate" = all(names(knownComp_param) == c("shape","rate")))
  } else {
    if (!all(as.character(dist_table[rownames(dist_table) == knownComp_dist, 4:(4+nparam_theo-1)]) == names(knownComp_param))) {
      cat("Name of parameters not appropriate, please provide the following parameters /",
          as.character(dist_table[rownames(dist_table) == knownComp_dist, 4:(4+nparam_theo-1)]), sep = " / ")
      cat("\n")
      stop()
    }
  }

  ## Create object:
  obj <- list(
    comp.dist = list(unknown = NULL, known = knownComp_dist),
    comp.param = list(unknown = NULL, known = knownComp_param)
  )
  class(obj) <- c("admix_model", "twoComp_mixt")
  obj$call <- match.call()
  return(obj)
}


#' Print method for objects of class 'admix_model'
#'
#' Print an object of class 'admix_mod'. An admixture model has probability density function (pdf) l_i such that:
#'    l_i = p_i * f_i + (1-p_i) * g_i, with g_i the known component density.
#' The unknown quantities are therefore p_i and f_i.
#'
#' @param x An object of class 'admix_model'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_model <- function(x, ...)
{
  cat("Call:")
  print(x$call)
  cat("\n")
  cat("Known component distribution: ", x$comp.dist$known, "\n")
  cat("Value of the known component parameters:\n")
  print(x$comp.param$known)
}
