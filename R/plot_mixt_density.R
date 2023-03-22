#' Plot the density of some given sample(s) with mixture distributions.
#'
#' Plot the density of the sample(s) with optional arguments to improve the visualization.
#'
#' @param samples Observed samples (mixture distributions) from which the density will be plotted.
#' @param user.bounds (default to NULL) Bounds to limit the range of x-axis when plotting.
#' @param support Support of the distributions, to know whether density plot or histogram should be displayed.
#' @param main Title for the plot.
#'
#' @return a plot with the densities of the samples provided as inputs.
#'
#' @examples
#' ##### Continuous support:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm",
#'                   f3 = "norm", g3 = "norm")
#' list.param <- list(f1 = list(mean = 5, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 0, sd = 1), g2 = list(mean = -3, sd = 1.1),
#'                    f3 = list(mean = 9, sd = 1), g3 = list(mean = 6, sd = 2))
#' ## Simulate data:
#' sim1 <- rsimmix(n = 300, unknownComp_weight = 0.8, comp.dist = list(list.comp$f1,list.comp$g1),
#'                 comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' sim2 <- rsimmix(n= 250, unknownComp_weight = 0.85, comp.dist = list(list.comp$f2,list.comp$g2),
#'                 comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' sim3 <- rsimmix(n= 400, unknownComp_weight = 0.6, comp.dist = list(list.comp$f3,list.comp$g3),
#'                 comp.param = list(list.param$f3, list.param$g3))$mixt.data
#' plot_mixt_density(samples = list(sim1,sim2,sim3), user.bounds = NULL, support = "continuous")
#'
#' ####### Countable discrete support:
#' list.comp <- list(f1 = 'pois', g1 = 'pois',
#'                   f2 = 'pois', g2 = 'pois')
#' list.param <- list(f1 = list(lambda = 7), g1 = list(lambda = 1),
#'                    f2 = list(lambda = 2), g2 = list(lambda = 15))
#' sim1 <- rsimmix(n=4000, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                 comp.param=list(list.param$f1,list.param$g1))$mixt.data
#' sim2 <- rsimmix(n=3500, unknownComp_weight=0.3, comp.dist = list(list.comp$f2,list.comp$g2),
#'                 comp.param=list(list.param$f2,list.param$g2))$mixt.data
#' plot_mixt_density(samples = list(sim1,sim2), user.bounds = NULL, support = "discrete")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot_mixt_density <- function(samples, user.bounds = NULL, support = c("continuous","discrete"), main = "")
{
  n <- max(sapply(X = samples, FUN = length))
  n_pop <- length(samples)
  leg <- paste("Sample ", 1:n_pop, sep="")
  #leg <- eval(parse(text = paste("latex2exp::TeX('X_", 1:n_pop, "')", sep="")))

  densities <- lapply(X = samples, FUN = stats::density)
  x.axe.lim.old <- x.axe.lim <- c(min(apply(sapply(densities, "[[", "x"), 2, min)), max(apply(sapply(densities, "[[", "x"), 2, max)))

  if (!is.null(user.bounds)) x.axe.lim <- c(user.bounds[1],user.bounds[2])

  old_par <- graphics::par()$new
  if (support != "discrete") {
    y.axe.lim <- max(apply(sapply(densities, "[[", "y"), 2, max))
    for (i in 1:n_pop) {
      base::plot(densities[[i]], xlim = x.axe.lim, ylim = c(0,y.axe.lim), xlab = "", main = "", col = i, lty = i)
      graphics::par(new = TRUE)
    }
  } else {
    ## FIXME: pour l'argument 'breaks' defini comme ci-dessous, cela ne marche que si la loi discrete est a support positif!
    #?barplot
    for (i in 1:n_pop) {
      graphics::hist(samples[[i]], freq = TRUE, breaks = 0:ceiling(x.axe.lim.old[2]), xlim = x.axe.lim, ylim = c(0,n/2.3),
                     xlab = "", main = "", border = i, lty = i)
      graphics::par(new = TRUE)
    }
  }
  graphics::title(main = main)
  graphics::legend("topright", legend = leg, lty = 1:n_pop, col = 1:n_pop, bty = "n")

  on.exit(graphics::par(new = old_par))
}
