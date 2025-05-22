#' Define the distribution/parameter(s) of the known component
#'
#' Create an object of class 'admix_model', containing the information about the known component distribution in the admixture model.
#' An admixture (aka contamination) model is a two-component mixture model with one known component.
#' Both the second component distribution and the mixing weight are unknown.
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
#' admix_model(knownComp_dist = "pois", knownComp_param = list("lambda"=5))
#' admix_model(knownComp_dist = "multinom", knownComp_param = list("size"=1, "prob"=c(0.2,0.8,0.1)))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_model <- function(knownComp_dist, knownComp_param)
{
  stopifnot("Specify only ONE distribution" = length(knownComp_dist) == 1)
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
  class(obj) <- "admix_model"
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
  cat("Known component parameters:",
      paste(names(x$comp.param$known), "=", unlist(x$comp.param$known, use.names=FALSE), sep=""))
  cat("\n\n")
}


#' Plot method for objects of class 'admix_model'
#'
#' Plots the probability density function of the known component of the admixture model, where
#' we recall that an admixture model has probability density function (pdf) l_i such that:
#'    l_i = p_i * f_i + (1-p_i) * g_i, with g_i the known component density.
#' The unknown quantities are therefore p_i and f_i.
#'
#' @param x An object of class 'admix_model'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot.admix_model <- function(x, ...)
{
  sim_txt <- paste("r", x$comp.dist$known, "(n=10000, ",
                   paste(names(x$comp.param$known), "=", x$comp.param$known,
                         collapse=", ", sep=""), ")", sep = "", collapse = "")
  sim_values <- eval(parse(text = sim_txt))

  if (length(sim_values) == 10000) {
    x_range <- c(min(sim_values), max(sim_values))
    supp <- detect_support_type(sample1 = sim_values)
    if (supp == "Continuous") {
      txt <- paste("d", x$comp.dist$known, "(x=seq(from=",x_range[1], ", to=", x_range[2], ", length.out=1000),",
                   paste(names(x$comp.param$known), "=", x$comp.param$known, collapse=", ", sep=""), ")",
                   sep = "", collapse = "")
      plot(x = seq(from = x_range[1], to = x_range[2], length.out = 1000),
           y = eval(parse(text = txt)), ...)
    } else {  # the discrete case
      txt <- paste("d", x$comp.dist$known, "(x=seq(from=",x_range[1], ", to=", x_range[2], ", by=1),",
                   paste(names(x$comp.param$known), "=", x$comp.param$known, collapse=", ", sep=""), ")",
                   sep = "", collapse = "")
      plot(x = seq(from = x_range[1], to = x_range[2], by = 1), y = eval(parse(text = txt)), ...)
    }

  } else {  # case of multinomial distribution
    x_range <- 1:length(rowSums(sim_values))
    graphics::barplot(height = rowSums(sim_values)/sum(rowSums(sim_values)), names = as.character(x_range), space = 0.1, ...)
  }
}


#' Summary method for objects of class 'admix_model'
#'
#' Summarizes the information related to the known component of the two-component mixture.
#' An admixture model has probability density function (pdf) l_i such that:
#'    l_i = p_i * f_i + (1-p_i) * g_i, with g_i the known component density.
#' The unknown quantities are therefore p_i and f_i.
#'
#' @param object An object of class 'admix_model'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.admix_model <- function(object, ...)
{
  cat("Call:")
  print(object$call)
  cat("\n")
  cat("Known component distribution: ", object$comp.dist$known, "\n")
  cat("Value of the known component parameters:\n")
  print(object$comp.param$known)
  cat("\n")
}
