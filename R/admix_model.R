#' Define the distribution/parameter(s) of the known component
#'
#' Create an object of class \code{admix_model}, containing the information about the known component distribution in the admixture model.
#' An admixture (aka contamination) model is a two-component mixture model with one known component.
#' Both the second component distribution and the mixing weight are unknown.
#'
#' @param knownComp_dist (Character) The name of the distribution (specified as in R glossary) of the known component
#'                        of the admixture model.
#' @param knownComp_param (Character) A list of the names of the parameters (specified as in R glossary) involved in
#'                        the chosen known distribution, with their values.
#'
#' @return An object of class \link[admix]{admix_model}, containing 2 attributes: 1) a list that gives the information about the distributions
#'         involved in the two-component mixture model (the unknown and the known ones); 2) a list that gives the information about
#'         the corresponding parameters of those distributions.
#'
#' @examples
#' admix_model(knownComp_dist = "norm", knownComp_param = list("mean"=0, "sd"=1))
#' admix_model(knownComp_dist = "exp", knownComp_param = list("rate"=2))
#' admix_model(knownComp_dist = "pois", knownComp_param = list("lambda"=5))
#' admix_model(knownComp_dist = "multinom", knownComp_param = list("size"=1, "prob"=c(0.1,0.8,0.1)))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_model <- function(knownComp_dist, knownComp_param)
{
  if (!is.character(knownComp_dist) || length(knownComp_dist) != 1) {
    stop("`knownComp_dist` must be a character string of length 1.")
  }
  if (!is.list(knownComp_param) || is.null(names(knownComp_param))) {
    stop("`knownComp_param` must be a named list.")
  }

  ## --- Distribution / parameter validation -------
  validate_distribution(dist = knownComp_dist, params = knownComp_param)

  ## --- Object creation ------
  obj <- structure(
    list(
      comp.dist = list(unknown = NULL, known = knownComp_dist),
      comp.param = list(unknown = NULL, known = knownComp_param),
      call = match.call()
    ),
    class = "admix_model"
  )
  return(obj)
}


#' Print method for objects of class \code{admix_model}
#'
#' Print the information about the distribution of the known component of an admixture model,
#' as well as the known parameters for this distribution.
#'
#' @param x An object of class \code{admix_model}.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_model <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nAdmixture model\n")
  cat("----------------\n")
  cat("Known component distribution: ", x$comp.dist$known,"\n",sep = "")
  cat("Known parameters:\n")
  for (nm in names(x$comp.param$known)) {
    value <- x$comp.param$known[[nm]]
    if (length(value) > 1) {
      value <- paste(value, collapse = ", ")
    }
    cat("  - ", nm, " = ", value, "\n", sep = "")
  }
  invisible(x)
}


#' Plot method for objects of class \code{admix_model}
#'
#' Plots the probability density function of the known component of the admixture model.
#'
#' @param x An object of class \code{admix_model}.
#' @param n The number of 'x' values to consider for plotting the pdf in the continuous case.
#' @param main The title of the plot.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @examples
#' plot(admix_model(knownComp_dist = "norm", knownComp_param = list("mean"=0, "sd"=1)))
#' plot(admix_model(knownComp_dist = "pois", knownComp_param = list("lambda"=1.5)))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot.admix_model <- function(x, n = 1000, main = "Known component distribution", ...)
{
  dist <- x$comp.dist$known
  params <- x$comp.param$known
  multivariate_dist <- dist %in% c("multinom")
  if (multivariate_dist) {
    x_range <- seq_along(params$prob)
    graphics::barplot(height = params$prob, names = as.character(x_range),
                      space = 0.1, main = main, ...)
  } else {
    qfun <- get(paste0("q", dist), mode = "function")
    dfun <- get(paste0("d", dist), mode = "function")
    x_min <- do.call(qfun, c(list(p = 0.001), params))
    x_max <- do.call(qfun, c(list(p = 0.999), params))
    discrete_dist <- dist %in% c("binom","pois","geom","nbinom","hyper")
    if (discrete_dist) {
      xx <- seq(floor(x_min), ceiling(x_max), by = 1)
      yy <- do.call(dfun, c(list(x = xx), params))
      plot(xx, yy, type = "h", lwd = 2, xlab = "x", ylab = "probability", main = main, ...)
    } else {
      xx <- seq(x_min, x_max, length.out = n)
      yy <- do.call(dfun, c(list(x = xx), params))
      plot(xx, yy, type = "l", lwd = 2, xlab = "x", ylab = "density", main = main, ...)
    }
  }
  invisible(x)
}

#' Summary method for objects of class \code{admix_model}
#'
#' Summarizes the information related to the known component of the two-component mixture.
#'
#' @param object An object of class \code{admix_model}.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.admix_model <- function(object, ...)
{
  print(object)
}
