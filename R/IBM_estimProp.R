#' Estimate the weights related to the proportions of the unknown components of the two admixture models
#'
#' Estimate the component weights from the Inversion - Best Matching (IBM) method, related to the two admixture models
#' with respective probability density function (pdf) l1 and l2, such that:
#'    l1 = p1*f1 + (1-p1)*g1 and l2 = p2*f2 + (1-p2)*g2, where g1 and g2 are the known component densities.
#' For further details about IBM approach, see 'Details' below.
#'
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param known.prop (optional) Numeric vector with two elements, respectively the component weight for the unknown component in the first and in the second samples.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param with.correction Boolean indicating whether the solution (estimated proportions) should be adjusted or not
#'                       (with the constant determined thanks to the exact proportion, usually unknown except in case of simulations).
#' @param n.integ Number of data points generated for the distribution on which to integrate.
#'
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03201760
#'
#'
#' @return A list with the two estimates of the component weights for each of the admixture model, plus that of the theoretical model if specified.
#'
#' @examples
#' ##### On a simulated example to see whether the true parameters are well estimated.
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                  comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                  comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weights of the two admixture models (provide hat(theta)_n and theta^c):
#' estim <- IBM_estimProp(sample1 = sample1[['mixt.data']], sample2 = sample2[['mixt.data']],
#'                        known.prop = c(0.5,0.7), comp.dist = list.comp, comp.param = list.param,
#'                        with.correction = FALSE, n.integ = 1000)
#' estim[['prop.estim']]
#' estim[['theo.prop.estim']]
#' ##### On a real-life example (unknown component densities, unknown mixture weights).
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' ## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
#' estim <- IBM_estimProp(sample1 = sample1[['mixt.data']], sample2 = sample2[['mixt.data']],
#'                        known.prop = NULL, comp.dist = list.comp, comp.param = list.param,
#'                        with.correction = FALSE, n.integ = 1000)
#' estim[['prop.estim']]
#' estim[['theo.prop.estim']]
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_estimProp <- function(sample1, sample2, known.prop = NULL, comp.dist = NULL, comp.param = NULL, with.correction = TRUE, n.integ = 1000)
{
  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  ##------- Defines the support for integration by simulation --------##
  ## Allows to integrate the gap between F1 and F2 in the contrast computation. span(G) must contain span(X) and span(Y).
  ## Ideally, this distribution should put more weight to locations where differences between F1 and F2 are expected.
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    fit.allObs <- stats::density(c(sample1, sample2))
    G <- stats::rnorm(n.integ, sample(c(sample1, sample2), size = n.integ, replace = TRUE), fit.allObs$bw)
  } else {
    G <- unique(sort(c(unique(sample1), unique(sample2))))
  }

  if (G1equalG2) {
    ## leads to a one-dimensional optimization because known components of mixture distributions are the same:
    if (!is.null(known.prop)) {
      fixed.p.X <- known.prop[1]
    } else {
      fixed.p.X <- 0.3
    }
    par.init <- 0.5
  } else {
    ## Otherwise classical optimization to get the two component weights:
    fixed.p.X <- NULL
    par.init <- c(0.5,0.5)
  }

  ##------- Solution search by optimization --------##
  ## Use method 'L-BFGS-B' to optimize under constraints (on parameter space), and use 'try' to catch errors coming from divergent integrals (numerical issues):
  expr1_NM <- expression(stats::optim(par = par.init, fn = IBM_empirical_contrast, gr = NULL, fixed.p.X = fixed.p.X, sample1 = sample1,
                                      sample2 = sample2, G = G, comp.dist = comp.dist, comp.param = comp.param, method = "Nelder-Mead",
                                      control = list(trace = 0, maxit = 10000)))
  sol <- try(suppressWarnings(eval(expr1_NM)), silent = TRUE)
  count_error <- 0
  while ((class(sol) == "try-error") & (count_error < 3)) {
    sol <- NULL
    sol <- try(suppressWarnings(eval(expr1_NM)), silent = TRUE)
    count_error <- count_error + 1
  }

  ## To deal with extreme values that can be found and that cause numerical issues afterwards:
  if (any(abs(sol$par) > 5) | (class(sol) == "try-error")) {
    warning("In 'IBM_estimProp': optimization algorithm was changed (in 'optim') from 'Nelder-Mead' to 'BFGS' to avoid the solution to explose.")
    expr1_BFGS <- expression(stats::optim(par = par.init, fn = IBM_empirical_contrast, gr = NULL, fixed.p.X = fixed.p.X, sample1 = sample1,
                                          sample2 = sample2, G = G, comp.dist = comp.dist, comp.param = comp.param, method = "L-BFGS-B",
                                          lower = c(0.001,0.001), upper = c(10,10), control = list(trace = 0, maxit = 10000)))
    sol <- try(suppressWarnings(eval(expr1_BFGS)), silent = TRUE)
    count_error <- 0
    while ((class(sol) == "try-error") & (count_error < 3)) {
      sol <- NULL
      sol <- try(suppressWarnings(eval(expr1_BFGS)), silent = TRUE)
      count_error <- count_error + 1
    }
    if (class(sol) == "try-error") {
      warning("In 'IBM_estimProp': impossible to estimate the component weights with BFGS method. Switch back to Nelder-Mead algorithm to obtain a solution")
      sol <- eval(expr1_NM)
    }
  }

  estim.weights <- sol$par

  ## In case the underlying theoretical model is known, the theoretical contrast can be computed. Here, we need to specify all the component weights
  ## and all the component distributions, but there is no need to provide any observations.
  if (!is.null(known.prop)) {
    expr2_NM <- expression(stats::optim(par = par.init, fn = IBM_theoretical_contrast, gr = NULL, theo.par = known.prop, fixed.p.X = fixed.p.X,
                                        G = G, comp.dist = comp.dist, comp.param = comp.param, sample1 = sample1, sample2 = sample2,
                                        method = "Nelder-Mead", control = list(trace = 0, maxit = 10000)))
    sol.theo <- try(suppressWarnings(eval(expr2_NM)), silent = TRUE)
    count_error <- 0
    while ((class(sol.theo) == "try-error") & (count_error < 3)) {
      sol.theo <- NULL
      sol.theo <- try(suppressWarnings(eval(expr2_NM)), silent = TRUE)
      count_error <- count_error + 1
    }

    if (any(abs(sol.theo$par) > 5) | (class(sol.theo) == "try-error")) {
      warning("In 'IBM_estimProp': optimization algorithm was changed (in 'optim') from 'Nelder-Mead' to 'BFGS' to avoid the solution to explose.")
      expr2_BFGS <- expression(stats::optim(par = par.init, fn = IBM_theoretical_contrast, gr = NULL, theo.par = known.prop, fixed.p.X = fixed.p.X,
                                            G = G, comp.dist = comp.dist, comp.param = comp.param, sample1 = sample1, sample2 = sample2,
                                            method = "L-BFGS-B", lower = c(0.001,0.001), upper = c(10,10), control = list(trace = 0, maxit = 10000)))
      sol.theo <- try(suppressWarnings(eval(expr2_BFGS)), silent = TRUE)
      count_error <- 0
      while ((class(sol) == "try-error") & (count_error < 3)) {
        sol.theo <- NULL
        sol.theo <- try(suppressWarnings(eval(expr2_BFGS)), silent = TRUE)
        count_error <- count_error + 1
      }
      if (class(sol.theo) == "try-error") {
        warning("In 'IBM_estimProp': impossible to estimate the component weights with BFGS method. Switch back to Nelder-Mead algorithm to obtain a solution")
        sol.theo <- eval(expr2_NM)
      }
    }
    theo.weights <- sol.theo$par
  } else {
    theo.weights <- NULL
  }

  ##------- ## Adjustment of the fitted parameters when G1 = G2 --------##
  if (G1equalG2 & with.correction) {
    stopifnot(!is.null(known.prop))
    estim.weights <- estim.weights * (known.prop[1] / fixed.p.X)
    theo.weights <- theo.weights * (known.prop[1] / fixed.p.X)
  }

  ##------- Returns solution --------##
  solution <- list(integ.supp = sort(G), p.X.fixed = fixed.p.X, prop.estim = estim.weights, theo.prop.estim = theo.weights)
  return(solution)
}
