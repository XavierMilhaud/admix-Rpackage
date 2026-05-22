#' Simulation of a two-component mixture model
#'
#' Simulate a two-component mixture model following the probability density function (pdf) \eqn{\ell} such that
#' \deqn{
#'   \ell = p f + (1 - p) g,
#' }
#' with \eqn{f} and \eqn{g} the mixture component distributions, and \eqn{p} the mixing weight.
#'
#' @param n Number of observations to be simulated.
#' @param weight Weight of the first component distribution (distribution f) in the mixture.
#' @param comp.dist A list of two elements corresponding to the component distributions (with available names listed in object 'Distribution.df'
#'                  in package EnvStats) involved in the mixture model. These elements respectively refer to the two component distributions f and g.
#'                  By convention, in the framework of admixture models where one of the two components is unknown, the first element of the list
#'                  corresponds to the 'unknown' component distribution, whereas the second one refers to the known one.
#' @param comp.param A list of two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in each list must correspond to the available parameters listed in object 'Distribution.df'
#'                   in package EnvStats. These elements respectively refer to the parameters of f and g distributions of the mixture model.
#'                   By convention, in the framework of admixture models where one of the two components is unknown, the first element of the list
#'                   corresponds to the 'unknown' component parameters, whereas the second one refers to the known ones.
#'
#' @return An object of class \link[admix]{twoComp_mixt}, containing eight attributes: 1) the number of simulated observations, 2) the simulated mixture
#'         data, 3) the support of the distributions, 4) the name of the component distributions, 5) the name of the parameters of the
#'         component distributions and their values, 6) the mixing proportion, 7) the observations coming from the first component,
#'         8) the observations coming from the second component.
#'
#' @seealso [get_mixture_data()] to access the simulated mixture data.
#'
#' @examples
#' ## Mixture of continuous random variables:
#' sim.X <- twoComp_mixt(n = 1200, weight = 0.7,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean"=-3, "sd"=0.5),
#'                                         list("rate"=1)))
#' print(sim.X)
#' data.X <- get_mixture_data(sim.X)
#' plot(density(data.X))
#'
#' ## Mixture of discrete random variables:
#' sim.Y <- twoComp_mixt(n = 1800, weight = 0.7,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size"=1, "prob"=c(0.3,0.4,0.3)),
#'                                         list("size"=1, "prob"=c(0.6,0.2,0.2))))
#' print(sim.Y)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

twoComp_mixt <- function(n = 1000, weight = 0.5, comp.dist = list("norm", "norm"),
                         comp.param = list(list(mean = 0, sd = 1), list(mean = 2, sd = 1)))
{
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n %% 1 != 0) { stop("`n` must be a positive integer.") }
  if (!is.numeric(weight) || weight <= 0 || weight >= 1) { stop("`weight` must belong to (0,1).") }
  if (length(comp.dist) != 2 || length(comp.param) != 2) { stop("Please provide exactly two component distributions.") }

  for (k in 1:2) { validate_distribution(dist = comp.dist[[k]], params = comp.param[[k]] ) }

  ## Distribution type
  dist.type <- vapply(comp.dist, distribution_type, character(1))
  ## Compatibility checks
  if ("multinom" %in% comp.dist && !all(comp.dist == "multinom")) {
    stop("`multinom` can only be mixed with another multinomial distribution.")
  }
  if ("gompertz" %in% comp.dist && !all(comp.dist == "gompertz")) {
    stop("`gompertz` can only be mixed with another Gompertz distribution.")
  }

  ## Random generation
  rfun <- lapply(comp.dist, function(d) get(paste0("r", d), mode = "function"))
  ## Component labels
  z <- sample(x = 1:2, size = n, replace = TRUE, prob = c(weight, 1 - weight))
  ## Generate observations
  res <- vector("list", n)
  for (i in seq_len(n)) {
    k <- z[i]
    res[[i]] <- do.call(rfun[[k]], c(list(n = 1), comp.param[[k]]) )
  }

  ## Output
  if (all(comp.dist == "multinom")) {
    res_tmp <- vapply(res, function(x) which(x == 1), integer(1))
    res <- res_tmp
    #res_tmp <- rowSums(res)
    #res <- unlist( apply( as.data.frame(1:length(res_tmp)), 1, function(k) { rep(k, res_tmp[k]) } ) )
  } else {
    res <- unlist(res)
  }

  obj_res <- list(
    n = n,
    mixt.data = res,
    dist.type = dist.type,
    comp.dist = comp.dist,
    comp.param = comp.param,
    mix.prop = weight,
    comp1.data = res[z == 1],
    comp2.data = res[z == 2],
    call = match.call()
  )
  class(obj_res) <- "twoComp_mixt"
  obj_res
}


#' Print method for objects \code{twoComp_mixt}
#'
#' Print an object of class \code{twoComp_mixt}. A two-component mixture model has probability density function (pdf) l such that:
#'    l = p * f + (1-p) * g,
#' where p is the mixing proportion, and f and g are the component distributions.
#'
#' @param x An object of class \code{twoComp_mixt}.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.twoComp_mixt <- function(x, ...)
{
  cat("\nCall:")
  print(x$call)
  cat("\n")
  cat("Number of observations:", x$n, "\n")
  cat("\n")
  if (any(x$comp.dist == "multinom")) {
    cat("Obtained multinomial mixture distribution: \n", table(x$mixt.data), "\n")
  } else {
    cat("Simulated data (first 5 obs.): \n", utils::head(x$mixt.data, 5))
    cat("\n")
    cat("Simulated observations coming from the 1st component (first 5 obs.): \n", utils::head(x$comp1.data, 5))
    cat("\n")
    cat("Simulated observations coming from the 2nd component (first 5 obs.): \n", utils::head(x$comp2.data, 5))
  }
  cat("\n\n")
}


#' Summary method for objects \code{twoComp_mixt}
#'
#' Provides statistical indicators of an object of class \code{twoComp_mixt}.
#' A two-component mixture model has probability density function (pdf) l such that:
#'    l = p * f + (1-p) * g,
#' where p is the mixing proportion, and f and g are the component distributions.
#'
#' @param object An object of class \code{twoComp_mixt}.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.twoComp_mixt <- function(object, ...)
{
  cat("\nCall:")
  print(object$call)
  cat("\n")
  cat("Component distributions: ", unlist(object$comp.dist), "\n")
  cat("Support of the component distributions: ", unlist(object$dist.type), "\n\n")
  cat("Value of the component parameters: \n")
  print(unlist(object$comp.param))
  cat("\nMixing proportion:", object$mix.prop, "\n")
  cat("\n")
  cat("Number of observations: ", object$n, "\n")
  cat("\n")
  if (any(object$comp.dist == "multinom")) {
    cat("Obtained multinomial mixture distribution: \n", table(object$mixt.data), "\n")
  } else {
    cat("Statistics related to the simulated sample: \n")
    print(summary(object$mixt.data))
    cat("\n")
    cat("Statistics related to the first component of the two-component mixture: \n")
    print(summary(object$comp1.data))
    cat("\n")
    cat("Statistics related to the second component of the two-component mixture: \n")
    print(summary(object$comp2.data))
  }
  cat("\n")
}


#' Plot the empirical mixture pdf
#'
#' Plots the empirical densities of the samples provided, with optional arguments to improve the visualization.
#'
#' @param x Object of class \code{twoComp_mixt} from which the density will be plotted.
#' @param add_plot (default to FALSE) Option to plot another mixture distribution on the same graph.
#' @param offset Numeric. Position of the bars relative to the labels on the x-axis.
#' @param bar_width Width of bars to be plotted.
#' @param main The title of the plot.
#' @param ... further classical arguments and graphical parameters for methods plot and hist.
#'
#' @return A plot with the densities of the samples provided as inputs.
#'
#' @examples
#' ## Mixture of continuous random variables:
#' sim.X <- twoComp_mixt(n = 2000, weight = 0.5,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean"=3, "sd"=0.5),
#'                                         list("mean"=0, "sd"=1)))
#' sim.Y <- twoComp_mixt(n = 1200, weight = 0.7,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean"=-3, "sd"=0.5),
#'                                         list("rate"=1)))
#' plot(sim.X, xlim=c(-5,5), ylim=c(0,0.5))
#' plot(sim.Y, add_plot = TRUE, xlim=c(-5,5), ylim=c(0,0.5), col = "red")
#' legend("topright", legend = c("sim.X","sim.Y"), col = c("black","red"),
#'        lty = rep(1,2), bty = "n")
#'
#' ## Mixture of discrete random variables:
#' sim.X <- twoComp_mixt(n = 2000, weight = 0.5,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size"=1, "prob"=c(0.3,0.4,0.3)),
#'                                         list("size"=1, "prob"=c(0.1,0.2,0.7))))
#' sim.Y <- twoComp_mixt(n = 1800, weight = 0.7,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size"=1, "prob"=c(0.3,0.4,0.3)),
#'                                         list("size"=1, "prob"=c(0.6,0.2,0.2))))
#' sim.Z <- twoComp_mixt(n = 1800, weight = 0.3,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size"=1, "prob"=c(0.2,0.1,0.7)),
#'                                         list("size"=1, "prob"=c(1/3,1/3,1/3))))
#' plot(sim.X, offset = -0.05, bar_width = 0.05, col = "steelblue")
#' plot(sim.Y, add_plot = TRUE, offset = 0, bar_width = 0.05, col = "orange")
#' plot(sim.Z, add_plot = TRUE, offset = +0.05, bar_width = 0.05, col = "red")
#' legend("topleft", legend = c("sim.X","sim.Y","sim.Z"), col = c("steelblue","orange","red"),
#'        lty = rep(1,3), bty = "n")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
plot.twoComp_mixt <- function(x, add_plot = FALSE, offset = 0, bar_width = 0.2,
                              main = "Mixture distribution (density or mass function)", ...)
{
  if (all(x$dist.type == "Discrete") | all(x$dist.type == "Multivariate")) {
    ## Discrete data and density
    freq <- as.numeric(table(x$mixt.data))
    x_val <- as.numeric(names(table(x$mixt.data)))
    heights <- freq / sum(freq)
    if (!add_plot) {
      ## Initialise graphic window
      plot(range(x_val), range(0, heights * 1.1), type="n", xaxt="n",
           xlab="support", ylab="probability mass", main = main, ...)
      graphics::axis(1, at=x_val, labels=as.character(x_val))
    }
    ## Bars
    for (i in seq_along(x_val)) {
      graphics::rect(xleft  = x_val[i] - bar_width/2 + offset,
                     xright = x_val[i] + bar_width/2 + offset,
                     ybottom = 0, ytop = heights[i], ...)
    }

  } else {
    ## Continuous case: densities
    densities <- stats::density(x$mixt.data)
    if (!add_plot) { plot(densities, xlab = "support", main = main, ...)
    } else {
      graphics::lines(densities, main = "", ...)
    }
  }
}
