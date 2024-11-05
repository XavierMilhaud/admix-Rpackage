#' Simulation of a two-component mixture model
#'
#' Simulate a two-component mixture model following the probability density function (pdf) l such that l = p*f + (1-p)*g,
#' with f and g the mixture component distributions, and p the mixing weight.
#'
#' @param n Number of observations to be simulated.
#' @param weight Weight of the first component distribution (distribution f) in the mixture.
#' @param comp.dist A list of two elements corresponding to the component distributions (specified with R native names)
#'                  involved in the mixture model. These elements respectively refer to the two component distributions f and g.
#' @param comp.param A list of two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in each list must correspond to the native R argument names for these distributions.
#'                   These elements respectively refer to the parameters of f and g distributions of the mixture model.
#'
#' @return An object of class 'twoComp_mixt', containing eight attributes: 1) the number of simulated observations, 2) the simulated mixture
#'         data, 3) the support of the distributions, 4) the name of the component distributions, 5) the name of the parameters of the
#'         component distributions and their values, 6) the mixing proportion, 7) the observations coming from the first component,
#'         8) the observations coming from the second component.
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
#' print(sim.X)
#' plot(sim.X, xlim=c(-5,5), ylim=c(0,0.5))
#' plot(sim.Y, add.plot = TRUE, xlim=c(-5,5), ylim=c(0,0.5), col = "red")
#'
#' ## Mixture of discrete random variables:
#' sim.X <- twoComp_mixt(n = 2000, weight = 0.5,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size"=1, "prob"=c(0.3,0.4,0.3)),
#'                                         list("size"=1, "prob"=c(0.1,0.2,0.7))))
#' print(sim.X)
#' plot(sim.X)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

twoComp_mixt <- function(n = 1000, weight = 0.5, comp.dist = list("norm", "norm"),
                         comp.param = list(list("mean"=0,"sd"=1), list("mean"=2,"sd"=1)))
{
  ## Some arguments to check:
  if ( !((weight > 0) & (weight < 1)) ) stop("Mixing proportion in a mixture model must belong to ]0,1[.")
  if ((length(comp.dist) != 2) | (length(comp.param) != 2)) stop("Please provide TWO distributions with corresponding parameters,
                                                                 one for each of the mixture components.")
  dist_table <- EnvStats::Distribution.df[ ,c("Name", "Type", "Number.parameters", "Parameter.1",
                                              "Parameter.2", "Parameter.3", "Parameter.4", "Parameter.5")]
  stopifnot("Unknown specified distribution" = any(comp.dist %in% rownames(dist_table)) | any(comp.dist == "multinom"))
  dist.type <- dist_table[match(unlist(comp.dist), rownames(dist_table)), "Type"]
  if (any(comp.dist == "multinom")) {
    dist.type <- c("Discrete", "Discrete")
    if (!(all(comp.dist == "multinom"))) stop("Multinomial distribution can only be mixed with another multinomial distribution")
  }

  nparam_theo <- numeric(length = 2L)
  for (k in 1:length(nparam_theo)) {
    if (any(comp.dist == "multinom")) { nparam_theo[k] <- 2
    } else { nparam_theo[k] <- dist_table[rownames(dist_table) == comp.dist[[k]], "Number.parameters"] }
  }

  stopifnot("Mispecification of parameters" = all(nparam_theo == sapply(comp.param, length)))
  for (k in 1:length(nparam_theo)) {
    if (any(comp.dist == "multinom")) {
      stopifnot("Name of parameters not appropriate" = all(names(comp.param[[k]]) == c("size","prob")))
    } else {
      if (!all(as.character(dist_table[rownames(dist_table) == comp.dist[[k]], 4:(4+nparam_theo[k]-1)]) == names(comp.param[[k]]))) {
        cat("Name of parameters not appropriate, please provide the following parameters /",
              as.character(dist_table[rownames(dist_table) == comp.dist[[k]], 4:(4+nparam_theo[k]-1)]), sep = " / ")
        cat("\n")
        stop()
      }
    }
  }

  ## Extracts the information on component distributions:
  comp.dist_sim <- paste0("r", comp.dist)
  comp_sim <- sapply(X = comp.dist_sim, FUN = get, mode = "function")
  for (i in 1:length(comp_sim)) assign(x = names(comp_sim)[i], value = comp_sim[[i]])
  ## Check if arguments of R core functions were correctly specified:
  #arg.names <- sapply(X = comp_sim, FUN = methods::formalArgs)

  ## Creates the expression allowing further to generate the right data:
  make.expr_sim <- function(i) {
    paste(names(comp_sim)[i], "(n=1,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  }
  expr_sim <- sapply(1:length(comp.dist), make.expr_sim)

  ## Generates the label for each observation:
  z <- sample(x = 2, size = n, replace = TRUE, prob = c(weight, 1-weight))

  ## Generates the final mixture data:
  data.gen <- parse(text = expr_sim[z])
  res <- sapply(data.gen, eval)
  if (any(comp.dist_sim == "rmultinom")) {
    res_tmp <- rowSums(res)
    res <- unlist( apply( as.data.frame(1:length(res_tmp)), 1, function(k) { rep(k, res_tmp[k]) } ) )
  }

  obj_res <- list(n = n,
                  mixt.data = res,
                  dist.type = dist.type,
                  comp.dist = comp.dist,
                  comp.param = comp.param,
                  mix.prop = weight,
                  comp1.data = res[z == 1],
                  comp2.data = res[z == 2]
                  )
  class(obj_res) <- "twoComp_mixt"
  obj_res$call <- match.call()
  return(obj_res)
}


#' Plots several mixture densities on the same graph
#'
#' Plots the empirical densities of the samples with optional arguments to improve the visualization.
#'
#' @param x Object of class 'twoComp_mixt' from which the density will be plotted.
#' @param add.plot (default to FALSE) Option to plot another mixture distribution on the same graph.
#' @param ... further classical arguments and graphical parameters for methods plot and hist.
#'
#' @return a plot with the densities of the samples provided as inputs.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot.twoComp_mixt <- function(x, add.plot = FALSE, ...)
{
  if (add.plot) {
    graphics::par(new = TRUE)
  }
  n <- x$n
  densities <- stats::density(x$mixt.data)
  x.axis.lim <- c(min(densities$x), max(densities$x))
  if (all(x$dist.type != "Discrete")) {
    base::plot(densities, ...)
  } else {
    ## FIXME: pour breaks' defini comme ci-dessous, cela ne marche que si la loi discrete est a support positif!
    #graphics::hist(x$mixt.data, freq = TRUE, breaks = 0:ceiling(x.axis.lim[2]), ...)
    base::plot(as.numeric(names(table(x$mixt.data))),
               as.numeric(table(x$mixt.data)) / sum(as.numeric(table(x$mixt.data))),
               type = "h", ...)
  }
  on.exit(graphics::par(new = FALSE))
}


#' Print method for objects 'twoComp_mixt'
#'
#' Print an object of class 'twoComp_mixt'. A two-component mixture model has probability density function (pdf) l such that:
#'    l = p * f + (1-p) * g,
#' where p is the mixing proportion, and f and g are the component distributions.
#'
#' @param x An object of class 'twoComp_mixt'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.twoComp_mixt <- function(x, ...)
{
  cat("\nCall:")
  print(x$call)
  cat("\n")
  cat("Component distributions: ", unlist(x$comp.dist), "\n")
  cat("Support of the component distributions: ", unlist(x$dist.type), "\n\n")
  cat("Value of the component parameters: \n")
  print(unlist(x$comp.param))
  cat("\nMixing proportion:", x$mix.prop, "\n")
  cat("\n")
  cat("Number of observations: ", x$n, "\n")
  cat("\n")
  if (any(x$comp.dist == "multinom")) {
    cat("Obtained multinomial mixture distribution: \n", table(x$mixt.data), "\n")
  } else {
    cat("Simulated data (first 10 obs.): \n", utils::head(x$mixt.data, 10), "\n")
    cat("\n")
    cat("Simulated observations coming from the 1st component (first 10 obs.): \n", utils::head(x$comp1.data, 10), "\n")
    cat("\n")
    cat("Simulated observations coming from the 2nd component (first 10 obs.): \n", utils::head(x$comp2.data, 10), "\n")
  }
  cat("\n")
}


#' Extractor for object of class 'twoComp_mixt'
#'
#' Get the mixture data generated from method 'twoComp_mixt'.
#'
#' @param x An object of class 'twoComp_mixt'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

getmixtData <- function(x)
{
  x$mixt.data
}
