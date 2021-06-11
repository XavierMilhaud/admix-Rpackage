#' Plot the decontaminated cumulative distribution function (CDF) of the unknown component for an estimated admixture model
#'
#' Plot the decontaminated CDF of the unknown component of the admixture model under study, after inversion of the admixture
#' CDF. Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with g a known CDF and p and f unknown quantities.
#'
#' @param x_val A vector of x-axis values at which to plot the decontaminated cumulative distribution function F.
#' @param decontamin_cdf An object from function 'decontamin_cdf_unknownComp', containing the unknown decontaminated
#'                       distribution function, as well as the support of the distribution (either discrete or continuous).
#' @param add_plot A boolean (TRUE by default) specifying if one plots the decontaminated density over an existing plot. Used for visual
#'                 comparison purpose.
#'
#' @details The decontaminated CDF is obtained by inverting the admixture CDF, given by L = p*F + (1-p)*G, to isolate the
#'          unknown component F after having estimated p and L.
#'
#' @return The plot of the decontaminated cumulative distribution function.
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
#' ## Determine the decontaminated version of the unknown CDF by inversion:
#' res1 <- decontamin_cdf_unknownComp(sample1 = sample1[['mixt.data']],
#'                            comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                            estim.p = estimate$prop.estim[1])
#' res2 <- decontamin_cdf_unknownComp(sample1 = sample2[['mixt.data']],
#'                            comp.dist = list.comp[3:4], comp.param = list.param[3:4],
#'                            estim.p = estimate$prop.estim[2])
#' plot_decontamin_cdf(x_val = seq(from=-1, to=5, length.out=30), decontamin_cdf = res1,
#'                     add_plot = FALSE)
#' plot_decontamin_cdf(x_val = seq(from=-1, to=5, length.out=30), decontamin_cdf = res2,
#'                     add_plot = TRUE)
#' ####### Countable discrete support:
#' list.comp <- list(f1 = 'pois', g1 = 'pois',
#'                   f2 = 'pois', g2 = 'pois')
#' list.param <- list(f1 = list(lambda = 3), g1 = list(lambda = 2),
#'                    f2 = list(lambda = 3), g2 = list(lambda = 4))
#' sample1 <- rsimmix(n=4000, unknownComp_weight=0.7, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=3000, unknownComp_weight=0.9, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'pois',
#'                   f2 = NULL, g2 = 'pois')
#' list.param <- list(f1 = NULL, g1 = list(lambda = 2),
#'                    f2 = NULL, g2 = list(lambda = 4))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' res1 <- decontamin_cdf_unknownComp(sample1 = sample1[['mixt.data']],
#'                            comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                            estim.p = estimate$prop.estim[1])
#' res2 <- decontamin_cdf_unknownComp(sample1 = sample2[['mixt.data']],
#'                            comp.dist = list.comp[3:4], comp.param = list.param[3:4],
#'                            estim.p = estimate$prop.estim[2])
#' plot_decontamin_cdf(x_val = seq(from=0, to=15, by=1), decontamin_cdf = res1,
#'                     add_plot = FALSE)
#' plot_decontamin_cdf(x_val = seq(from=0, to=15, by=1), decontamin_cdf = res2,
#'                     add_plot = TRUE)
#' ####### Finite discrete support:
#' list.comp <- list(f1 = 'multinom', g1 = 'multinom',
#'                   f2 = 'multinom', g2 = 'multinom')
#' list.param <- list(f1 = list(size=1, prob=c(0.3,0.4,0.3)), g1 = list(size=1, prob=c(0.6,0.3,0.1)),
#'                    f2 = list(size=1, prob=c(0.3,0.4,0.3)), g2 = list(size=1, prob=c(0.2,0.6,0.2)))
#' sample1 <- rsimmix(n=5000, unknownComp_weight=0.7, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=4000, unknownComp_weight=0.9, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' list.comp <- list(f1 = NULL, g1 = 'multinom',
#'                   f2 = NULL, g2 = 'multinom')
#' list.param <- list(f1 = NULL, g1 = list(size=1, prob=c(0.6,0.3,0.1)),
#'                    f2 = NULL, g2 = list(size=1, prob=c(0.2,0.6,0.2)))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' res1 <- decontamin_cdf_unknownComp(sample1 = sample1[['mixt.data']],
#'                            comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                            estim.p = estimate$prop.estim[1])
#' res2 <- decontamin_cdf_unknownComp(sample1 = sample2[['mixt.data']],
#'                            comp.dist = list.comp[3:4], comp.param = list.param[3:4],
#'                            estim.p = estimate$prop.estim[2])
#' plot_decontamin_cdf(x_val = seq(from=0, to=5, by=1), decontamin_cdf = res1,
#'                     add_plot = FALSE)
#' plot_decontamin_cdf(x_val = seq(from=0, to=5, by=1), decontamin_cdf = res2,
#'                     add_plot = TRUE)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot_decontamin_cdf <- function(x_val, decontamin_cdf, add_plot = FALSE)
{
  if (!add_plot) {
    stats::plot.stepfun(x = decontamin_cdf, xval = x_val)
  } else {
    old_par <- graphics::par()$new
    graphics::par(new = TRUE)
    stats::plot.stepfun(x = decontamin_cdf, xval = x_val, add = TRUE, col = "blue")
    on.exit(graphics::par(old_par))
  }
}
