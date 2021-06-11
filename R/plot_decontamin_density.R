#' Plot the decontaminated density of the unknown component for an estimated admixture model
#'
#' Plot the decontaminated density of the unknown component in the admixture model under study, after inversion of the admixture
#' cumulative distribution function. Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with g a known CDF and p and f unknown quantities.
#'
#' @param x_val A vector of X-axis values at which to plot the decontaminated density f.
#' @param decontamin_obj An object from function 'decontamin_density_unknownComp', containing the X-axis values and range for plotting,
#'                       the corresponding Y-axis values, and the support over which to plot the results.
#' @param add_plot (default to FALSE) A boolean specifying if one plots the decontaminated density over an existing plot. Used for visual
#'                 comparison purpose.
#'
#' @details The decontaminated density is obtained by inverting the admixture density, given by l = p*f + (1-p)*g, to isolate the
#'          unknown component f after having estimated p.
#'
#' @return The plot of the decontaminated density.
#'
#' @examples
#' ####### Continuous support:
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=3000, unknownComp_weight=0.7, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2500, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' res1 <- decontamin_density_unknownComp(sample1 = sample1[['mixt.data']],
#'                                       comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                                       estim.p = estimate$prop.estim[1])
#' res2 <- decontamin_density_unknownComp(sample1 = sample2[['mixt.data']],
#'                                        comp.dist = list.comp[3:4], comp.param = list.param[3:4],
#'                                        estim.p = estimate$prop.estim[2])
#' ## Use appropriate sequence of x values:
#' plot_decontamin_density(x_val = seq(from=0, to=6, length.out=100), decontamin_obj = res1,
#'                         add_plot = FALSE)
#' plot_decontamin_density(x_val = seq(from=0, to=6, length.out=100), decontamin_obj = res2,
#'                         add_plot = TRUE)
#' ####### Countable discrete support:
#' list.comp <- list(f1 = 'pois', g1 = 'pois',
#'                   f2 = 'pois', g2 = 'pois')
#' list.param <- list(f1 = list(lambda = 3), g1 = list(lambda = 2),
#'                    f2 = list(lambda = 3), g2 = list(lambda = 4))
#' sample1 <- rsimmix(n=4000, unknownComp_weight=0.7, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=3500, unknownComp_weight=0.85, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'pois',
#'                   f2 = NULL, g2 = 'pois')
#' list.param <- list(f1 = NULL, g1 = list(lambda = 2),
#'                    f2 = NULL, g2 = list(lambda = 4))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' res1 <- decontamin_density_unknownComp(sample1 = sample1[['mixt.data']],
#'                                       comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                                       estim.p = estimate$prop.estim[1])
#' res2 <- decontamin_density_unknownComp(sample1 = sample2[['mixt.data']],
#'                                        comp.dist = list.comp[3:4], comp.param = list.param[3:4],
#'                                        estim.p = estimate$prop.estim[2])
#' ## Use appropriate sequence of x values:
#' plot_decontamin_density(x_val = seq(from=0, to=15, by=1), decontamin_obj = res1,
#'                         add_plot = FALSE)
#' plot_decontamin_density(x_val = seq(from=0, to=15, by=1), decontamin_obj = res2,
#'                         add_plot = TRUE)
#' ####### Finite discrete support:
#' list.comp <- list(f1 = 'multinom', g1 = 'multinom',
#'                   f2 = 'multinom', g2 = 'multinom')
#' list.param <- list(f1 = list(size=1, prob=c(0.3,0.4,0.3)), g1 = list(size=1, prob=c(0.6,0.3,0.1)),
#'                    f2 = list(size=1, prob=c(0.3,0.4,0.3)), g2 = list(size=1, prob=c(0.2,0.6,0.2)))
#' sample1 <- rsimmix(n=4000, unknownComp_weight=0.8, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=3500, unknownComp_weight=0.9, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'multinom',
#'                   f2 = NULL, g2 = 'multinom')
#' list.param <- list(f1 = NULL, g1 = list(size=1, prob=c(0.6,0.3,0.1)),
#'                    f2 = NULL, g2 = list(size=1, prob=c(0.2,0.6,0.2)))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' res1 <- decontamin_density_unknownComp(sample1 = sample1[['mixt.data']],
#'                                       comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                                       estim.p = estimate$prop.estim[1])
#' res2 <- decontamin_density_unknownComp(sample1 = sample2[['mixt.data']],
#'                                        comp.dist = list.comp[3:4], comp.param = list.param[3:4],
#'                                        estim.p = estimate$prop.estim[2])
#' ## Use appropriate sequence of x values:
#' plot_decontamin_density(x_val = seq(from=0, to=6, by=1), decontamin_obj = res1,
#'                         add_plot = FALSE)
#' plot_decontamin_density(x_val = seq(from=0, to=6, by=1), decontamin_obj = res2,
#'                         add_plot = TRUE)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot_decontamin_density <- function(x_val, decontamin_obj, add_plot = FALSE)
{
  support <- decontamin_obj$supp
  decontamin_dens_values <- decontamin_obj$decontamin_f(x_val)
  decontamin_dens_values[which(is.na(decontamin_dens_values))] <- 0
  decontamin_dens_values <- pmax(decontamin_dens_values, 0)
  x.range <- c(min(x_val), max(x_val))
  y.range <- c(min(decontamin_dens_values), max(decontamin_dens_values)*1.1)

  if (!add_plot) {
    if (support == "discrete") {
      graphics::barplot(height = decontamin_dens_values, names = as.character(x_val),
                        xlim = x.range, ylim = y.range, col = "grey")
    } else {
      plot(x = x_val, y = decontamin_dens_values, xlim = x.range, ylim = y.range, type = "l")
    }

  } else {

    old_par <- graphics::par()$new
    graphics::par(new = TRUE)
    if (support == "discrete") {
      graphics::barplot(height = decontamin_dens_values, names = as.character(x_val),
                        add = TRUE, col = "blue")
    } else {
      graphics::lines(x = x_val, y = decontamin_dens_values, col = "blue")
    }
    on.exit(graphics::par(old_par))
  }

}
