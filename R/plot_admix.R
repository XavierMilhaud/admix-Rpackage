#' Plot the density of some given sample(s)
#'
#' Plot the density of the sample(s) with optional arguments to improve the visualization.
#'
#' @param sim.X First sample from which the density will be plotted.
#' @param sim.Y (default to NULL) Second sample from which the density will be plotted.
#' @param user.bounds (default to NULL) Bounds to limit the range of x-axis when plotting.
#' @param support Support of the distributions, to know whether density plot or histogram should be displayed.
#' @param case Used for titles.
#'
#' @return a plot with the densities of the samples provided as inputs.
#'
#' @examples
#' comp.dist <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' comp.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = -2, sd = 0.8), g2 = list(mean = 2, sd = 0.9))
#' sim.X <- rsimmix(n=2000, unknownComp_weight = 0.7, comp.dist = list(comp.dist$f1,comp.dist$g1),
#'                  comp.param = list(comp.param$f1,comp.param$g1))
#' sim.Y <- rsimmix(n=2000, unknownComp_weight = 0.4, comp.dist = list(comp.dist$f2,comp.dist$g2),
#'                  comp.param = list(comp.param$f2,comp.param$g2))
#' plot_admix(sim.X[['mixt.data']], sim.Y[['mixt.data']],
#'            user.bounds = c(-6,6), support = 'continuous')
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot_admix <- function(sim.X, sim.Y = NULL, user.bounds = NULL, support = c("continuous","discrete"), case = "")
{
  n <- max(length(sim.X),length(sim.Y))
  densite.X <- stats::density(sim.X)
  if (!is.null(sim.Y)) {
    densite.Y <- stats::density(sim.Y)
    x.axe.lim.old <- x.axe.lim <- c(min(c(densite.X$x,densite.Y$x)) * 0.95, max(c(densite.X$x,densite.Y$x)) * 1.05)
    leg <- c(latex2exp::TeX("X_1"), latex2exp::TeX("X_2"))
  } else {
    densite.Y <- NA
    x.axe.lim.old <- x.axe.lim <- c(min(densite.X$x) * 0.95, max(densite.X$x) * 1.05)
    leg <- latex2exp::TeX("X_1")
  }
  if (!is.null(user.bounds)) x.axe.lim <- c(user.bounds[1],user.bounds[2])

  if (support != "discrete") {
    if (!is.null(sim.Y)) {
      y.axe.lim <- max(c(densite.X$y,densite.Y$y)) * 1.05
      base::plot(densite.X, xlim = x.axe.lim, ylim = c(0,y.axe.lim), xlab = "", main = "")
      old_par <- graphics::par()$new
      graphics::par(new = TRUE)
      base::plot(densite.Y, xlim = x.axe.lim, ylim = c(0,y.axe.lim), xlab = "", main = "", col = "red", lty = 2)
      on.exit(graphics::par(old_par))
    } else {
      y.axe.lim <- max(densite.X$y) * 1.05
      base::plot(densite.X, xlim = x.axe.lim, ylim = c(0,y.axe.lim), xlab = "", main = "")
    }
  } else {
    if (!is.null(sim.Y)) {
      ## FIXME: pour l'argument 'breaks' defini comme ci-dessous, cela ne marche que si la loi discrete est a support positif!
#      ?barplot
      graphics::hist(sim.X, freq = TRUE, breaks = 0:ceiling(x.axe.lim.old[2]), xlim = x.axe.lim, ylim = c(0,n/2.3),
                     xlab = "", main = "")
      old_par <- graphics::par()$new
      graphics::par(new = TRUE)
      graphics::hist(sim.Y, freq = TRUE, breaks = 0:ceiling(x.axe.lim.old[2]), xlim = x.axe.lim, ylim = c(0,n/2.3),
                     xlab = "", main = "", border = "red", lty = 2)
      on.exit(graphics::par(old_par))
    } else {
      graphics::hist(sim.X, freq = TRUE, breaks = 0:ceiling(x.axe.lim.old[2]), xlim = x.axe.lim, ylim = c(0,n/2.3),
                     xlab = "", main = "")
    }
  }
  graphics::title(main = paste("Density of the two-component mixture model, case ", case, sep = ""))
  graphics::legend("topright", legend = leg, lty = 1:2, col = c("black","red"), bty = "n")
}

