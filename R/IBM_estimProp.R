#' Estimates weights of unknown components from 2 admixtures using IBM
#'
#' Estimation of the component weights from the Inversion - Best Matching (IBM) method, related to two admixture models
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
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
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

IBM_estimProp <- function(sample1, sample2, known.prop = NULL, comp.dist = NULL, comp.param = NULL,
                          with.correction = TRUE, n.integ = 1000)
{
  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  ##------- Defines the support for integration by simulation --------##
  ## Allows to integrate the gap between F1 and F2 in the contrast computation. span(G) must contain span(X) and span(Y).
  ## Ideally, this distribution should put more weight to locations where differences between F1 and F2 are expected.
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
#    if (contrast.weighting.scheme == "uniform") {
      G <- stats::runif(n.integ, min = min(c(sample1, sample2)), max = max(c(sample1, sample2)))
#    } else if (contrast.weighting.scheme == "empirical") {
#      fit.allObs <- stats::density(c(sample1, sample2))
#      G <- stats::rnorm(n.integ, sample(c(sample1, sample2), size = n.integ, replace = TRUE), fit.allObs$bw)
#    } else if (contrast.weighting.scheme == "unknown.comp")  {
##      warning("Implemented only under H0 (F1=F2), with centered gaussian unknown distrib")
#      G <- rnorm(n.integ, mean = 0, sd = 1)
#    } else if (contrast.weighting.scheme == "instrumental") { # Instrumental parametric disitribution
#      stopifnot(!is.null(param.instru))
#      G <- rnorm(n.integ, mean = param.instru[1], sd = param.instru[2])
#    } else {
#      a.quantile <- qnorm(p = 0.025, mean = 0, sd = 1)
#      b.quantile <- qnorm(p = 0.975, mean = 0, sd = 1)
#      G <- stats::runif(n.integ, min = a.quantile, max = b.quantile)
#    }

  } else {
    G <- unique(sort(c(unique(sample1), unique(sample2))))
  }

  if (G1equalG2) {
    ## leads to a one-dimensional optimization because known components of mixture distributions are the same:
    if (!is.null(known.prop)) {
      fixed.p.X <- known.prop[1]
    } else {
      ## set arbitrary the proportion to 1/2, could be any other value belonging to ]0,1[:
      fixed.p.X <- 0.5
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
  while ((inherits(x = sol, what = "try-error", which = FALSE)) & (count_error < 3)) {
    sol <- NULL
    sol <- try(suppressWarnings(eval(expr1_NM)), silent = TRUE)
    count_error <- count_error + 1
  }
  if (inherits(x = sol, what = "try-error", which = FALSE)) { sol <- list(par = 100) }

  ## To deal with extreme values that can be found and that cause numerical issues afterwards:
  if (any(abs(sol[['par']]) > 5)) {
    #message("In 'IBM_estimProp': optimization algorithm was changed (in 'optim') from 'Nelder-Mead' to 'BFGS' to avoid the solution to explose.")
    expr1_BFGS <- expression(stats::optim(par = par.init, fn = IBM_empirical_contrast, gr = NULL, fixed.p.X = fixed.p.X, sample1 = sample1,
                                          sample2 = sample2, G = G, comp.dist = comp.dist, comp.param = comp.param, method = "L-BFGS-B",
                                          lower = c(0.001,0.001), upper = c(5,5), control = list(trace = 0, maxit = 10000)))
    sol <- try(suppressWarnings(eval(expr1_BFGS)), silent = TRUE)
    count_error <- 0
    while (inherits(x = sol, what = "try-error", which = FALSE) & (count_error < 7)) {
      sol <- NULL
      sol <- try(suppressWarnings(eval(expr1_BFGS)), silent = TRUE)
      count_error <- count_error + 1
    }
    if (inherits(x = sol, what = "try-error", which = FALSE)) {
      #message("In 'IBM_estimProp': impossible to estimate the component weights with BFGS method. Switch back to Nelder-Mead algorithm to obtain a solution")
      sol <- try(suppressWarnings(eval(expr1_NM)), silent = TRUE)
    }
  }

  if (inherits(x = sol, what = "try-error", which = FALSE)) {
    stop("In 'IBM_estimProp': whatever the optimization algorithm, the solution has not been found.")
  } else {
    estim.weights <- sol[['par']]
  }

  ## In case the underlying theoretical model is known, the theoretical contrast can be computed. Here, we need to specify all the component weights
  ## and all the component distributions, but there is no need to provide any observations.
  if (!is.null(known.prop)) {
    expr2_NM <- expression(stats::optim(par = par.init, fn = IBM_theoretical_contrast, gr = NULL, theo.par = known.prop, fixed.p.X = fixed.p.X,
                                        G = G, comp.dist = comp.dist, comp.param = comp.param, sample1 = sample1, sample2 = sample2,
                                        method = "Nelder-Mead", control = list(trace = 0, maxit = 10000)))
    sol.theo <- try(suppressWarnings(eval(expr2_NM)), silent = TRUE)
    count_error <- 0
    while (inherits(x = sol, what = "try-error", which = FALSE) & (count_error < 3)) {
      sol.theo <- NULL
      sol.theo <- try(suppressWarnings(eval(expr2_NM)), silent = TRUE)
      count_error <- count_error + 1
    }
    if (inherits(x = sol, what = "try-error", which = FALSE)) { sol.theo <- list(par = 100) }

    ## To deal with extreme values that can be found and that cause numerical issues afterwards:
    if (any(abs(sol.theo[['par']]) > 5)) {
      #message("In 'IBM_estimProp': optimization algorithm was changed (in 'optim') from 'Nelder-Mead' to 'BFGS' to avoid the solution to explose.")
      expr2_BFGS <- expression(stats::optim(par = par.init, fn = IBM_theoretical_contrast, gr = NULL, theo.par = known.prop, fixed.p.X = fixed.p.X,
                                            G = G, comp.dist = comp.dist, comp.param = comp.param, sample1 = sample1, sample2 = sample2,
                                            method = "L-BFGS-B", lower = c(0.001,0.001), upper = c(5,5), control = list(trace = 0, maxit = 10000)))
      sol.theo <- try(suppressWarnings(eval(expr2_BFGS)), silent = TRUE)
      count_error <- 0
      while (inherits(x = sol, what = "try-error", which = FALSE) & (count_error < 7)) {
        sol.theo <- NULL
        sol.theo <- try(suppressWarnings(eval(expr2_BFGS)), silent = TRUE)
        count_error <- count_error + 1
      }
      if (inherits(x = sol, what = "try-error", which = FALSE)) {
        #message("In 'IBM_estimProp': impossible to estimate the component weights with BFGS method. Switch back to Nelder-Mead algorithm to obtain a solution")
        sol.theo <- try(suppressWarnings(eval(expr2_NM)), silent = TRUE)
      }
    }
    if (inherits(x = sol, what = "try-error", which = FALSE)) {
      stop("In 'IBM_estimProp': whatever the optimization algorithm, the theoretical solution has not been found.")
    } else {
      theo.weights <- sol.theo[['par']]
    }
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


#' Empirical contrast in the IBM method
#'
#' Computes the empirical version of the contrast in the Inversion - Best Matching (IBM) method, to be minimized in the optimization process.
#' For further details about the contrast definition, see 'Details' below.
#'
#' @param par Numeric vector with two elements, corresponding to the two parameter values at which to compute the contrast. In practice
#'            the component weights for the two admixture models.
#' @param fixed.p.X Arbitrary value chosen by the user for the component weight related to the unknown component of the first
#'                  admixture model. Only useful for optimization when the known components of the two models are identical
#'                  (G1=G2, leading to unidimensional optimization).
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param G Distribution on which to integrate when calculating the contrast.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return The empirical contrast value evaluated at parameter values.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                comp.param = list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                comp.param = list(list.param$f2,list.param$g2))
#' ## Create the distribution on which the contrast will be integrated:
#' G <- stats::rnorm(n = 1000, mean = sample(c(sample1[['mixt.data']], sample2[['mixt.data']]),
#'                                           size = 1000, replace = TRUE),
#'                   sd = density(c(sample1[['mixt.data']], sample2[['mixt.data']]))$bw)
#' ## Compute the empirical contrast at parameters (p1,p2) = (0.2,0.7) in a real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' IBM_empirical_contrast(par = c(0.2,0.7), fixed.p.X = NULL, sample1 = sample1[['mixt.data']],
#'            sample2= sample2[['mixt.data']], G=G, comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

IBM_empirical_contrast <- function(par, fixed.p.X = NULL, sample1, sample2, G, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
    stop("Known components must be specified in the two admixture models.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
  }

  ## Extract the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
  #comp_emp <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_emp <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_emp)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                                                                                        paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_emp)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Differentiates when the 4 components were specified (simulations) and real-life cases (only known components are given).
  ## First, real-life case:
  if (all(unlist(sapply(comp.param, is.na)[c(1,3)]))) {
    if (any(exp.comp.dist == "stepfun")) {
      G1.fun <- eval(parse(text = expr[2]))
      G2.fun <- eval(parse(text = expr[4]))
      G1 <- function(z) G1.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      G1 <- function(z) { eval(parse(text = expr[2])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  } else {
    ## case of simulated data:
    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      G1.fun <- eval(parse(text = expr[2]))
      F2.fun <- eval(parse(text = expr[3]))
      G2.fun <- eval(parse(text = expr[4]))
      F1 <- function(z) F1.fun(z)
      G1 <- function(z) G1.fun(z)
      F2 <- function(z) F2.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
      F2 <- function(z) { eval(parse(text = expr[3])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  }

  ## Empirical cumulative distribution function from the two observed samples:
  L1 <- stats::ecdf(sample1)						# hat(L1)
  L2 <- stats::ecdf(sample2)						# hat(L2)

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  ##------- Computes the empirical contrast --------##
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  if (is.null(fixed.p.X)) {
    integrand <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
      F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
      weighted.difference <- (F.X.dataPoint - F.Y.dataPoint)^2 * densite.G.dataPoint
      weighted.difference
    }
  } else {   # for one-dimensional optimization
    if (G1equalG2) stopifnot(!is.null(fixed.p.X))
    integrand <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/fixed.p.X) * (L1(z) - (1-fixed.p.X) * G1(z))
      F.Y.dataPoint <- (1/par) * (L2(z) - (1-par) * G2(z))
      weighted.difference <- (F.X.dataPoint - F.Y.dataPoint)^2 * densite.G.dataPoint
      weighted.difference
    }
  }

  if (support == "continuous") {
    res <- stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L,
                            rel.tol = 1e-03)$value
  } else {
    res <- sum( unlist(sapply(supp.integration, integrand, par)) )
  }

  return(res)
}

#' Theoretical contrast with IBM approach
#'
#' Defines the theoretical contrast used in the Inversion - Best Matching (IBM) approach. Useful in case of simulation studies,
#' since all parameters are known. For further information about the considered contrast in the IBM approach, see 'Details' below.
#'
#' @param par Numeric vector with two elements, corresponding to the two parameter values at which to compute the contrast. In practice
#'            the component weights for the two admixture models.
#' @param theo.par Numeric vector with two elements, the known (true) mixture weights.
#' @param fixed.p.X Arbitrary value chosen by the user for the component weight related to the unknown component of the first
#'                  admixture model. Only useful for optimization when the known components of the two models are identical
#'                  (G1=G2, leading to unidimensional optimization).
#' @param G Distribution on which to integrate when calculating the contrast.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. No unknown elements permitted.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1='rnorm', g1='norm', f2='rnorm', g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. No unknown elements permitted. For instance, 'comp.param' could be specified
#'                   as follows: : list(f1 = list(mean=2,sd=0.3), g1 = list(mean=0,sd=1), f2 = list(mean=2,sd=0.3), g2 = list(mean=3,sd=1.1)).
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return The theoretical contrast value evaluated at parameter values.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Create the distribution on which the contrast will be integrated:
#' G <- stats::rnorm(n = 1000, mean = sample(c(sample1[['mixt.data']],sample2[['mixt.data']]),
#'                                           size = 1000, replace = TRUE),
#'                   sd = stats::density(c(sample1[['mixt.data']],sample2[['mixt.data']]))$bw)
#' ## Compute the theoretical contrast at parameters (p1,p2) = (0.2,0.7):
#' IBM_theoretical_contrast(par = c(0.2,0.7), theo.par = c(0.5,0.7), fixed.p.X = NULL, G = G,
#'                          comp.dist = list.comp, comp.param = list.param,
#'                          sample1 = sample1[['mixt.data']], sample2 = sample2[['mixt.data']])
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

IBM_theoretical_contrast <- function(par, theo.par, fixed.p.X = NULL, G = NULL, comp.dist, comp.param, sample1, sample2)
{
  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null))) {
    stop("All components must be specified in the two admixture models to be able to compute the theoretical contrast.")
  }

  ## Extracts the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
  #  comp_theo <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_theo <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_theo)) assign(x = names(comp_theo)[i], value = comp_theo[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_theo)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                                                                                         paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_theo)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  if (any(exp.comp.dist == "stepfun")) {
    F1.fun <- eval(parse(text = expr[1]))
    G1.fun <- eval(parse(text = expr[2]))
    F2.fun <- eval(parse(text = expr[3]))
    G2.fun <- eval(parse(text = expr[4]))
    F1 <- function(z) F1.fun(z)
    G1 <- function(z) G1.fun(z)
    F2 <- function(z) F2.fun(z)
    G2 <- function(z) G2.fun(z)
  } else {
    F1 <- function(z) { eval(parse(text = expr[1])) }
    G1 <- function(z) { eval(parse(text = expr[2])) }
    F2 <- function(z) { eval(parse(text = expr[3])) }
    G2 <- function(z) { eval(parse(text = expr[4])) }
  }

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  ## Defines bounds on which to integrate:
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  if (is.null(fixed.p.X)) {
    integrand <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/par[1]) * ((theo.par[1]*F1(z) + (1-theo.par[1])*G1(z)) - (1-par[1])*G1(z))
      F.Y.dataPoint <- (1/par[2]) * ((theo.par[2]*F2(z) + (1-theo.par[2])*G2(z)) - (1-par[2])*G2(z))
      weighted.difference <- (F.X.dataPoint - F.Y.dataPoint)^2 * densite.G.dataPoint
      weighted.difference
    }
  } else {
    ## for one-dimensional optimization
    if (G1equalG2) stopifnot(!is.null(fixed.p.X))
    integrand <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/fixed.p.X) * ((theo.par[1]*F1(z) + (1-theo.par[1])*G1(z)) - (1-fixed.p.X)*G1(z))
      F.Y.dataPoint <- (1/par) * ((theo.par[2]*F2(z) + (1-theo.par[2])*G2(z)) - (1-par)*G2(z))
      weighted.difference <- (F.X.dataPoint - F.Y.dataPoint)^2 * densite.G.dataPoint
      weighted.difference
    }
  }

  if (support == "continuous") {
    res <- stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = 1e-03)$value
  } else {
    res <- sum( unlist(sapply(supp.integration, integrand, par)) )
  }

  return(res)
}


#' Hessian matrix of the contrast function when using IBM
#'
#' Compute the hessian matrix of the contrast as defined in the IBM approach, at point (p1,p2). Here, based on
#' two samples following admixture models, where we recall that admixture models have probability distribution function
#' (pdf) given by l where l = p*f + (1-p)*g, where g represents the only known quantity and l is the pdf of the observed sample.
#' See 'Details' below for further information about the definition of the contrast.
#'
#' @param par Numeric vector with two elements (corresponding to the two unknown component weights) at which the hessian is computed.
#' @param fixed.p1 (optional, NULL by default) Arbitrary value chosen by the user for the component weight related to the unknown
#'                  component of the first admixture model. Only useful for optimization when the known components of the two
#'                  models are identical (G1=G2, leading to unidimensional optimization).
#' @param known.p (optional, NULL by default) Numeric vector with two elements, the known (true) mixture weights.
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param G Distribution on which to integrate when calculating the contrast.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return the hessian matrix of the contrast.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Define the distribution over which to integrate:
#' fit.all <- stats::density(x = c(sample1[['mixt.data']],sample2[['mixt.data']]))
#' G <- stats::rnorm(n = 1000, mean = sample(c(sample1[['mixt.data']], sample2[['mixt.data']]),
#'                                           size = 1000, replace = TRUE), sd = fit.all$bw)
#' ## Evaluate the hessian matrix at point (p1,p2) = (0.3,0.6):
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' IBM_hessian_contrast(par = c(0.3,0.6), fixed.p1 = NULL, known.p = NULL,
#'                      sample1 = sample1[['mixt.data']],  sample2 = sample2[['mixt.data']], G = G,
#'                      comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

IBM_hessian_contrast <- function(par, fixed.p1 = NULL, known.p = NULL, sample1, sample2, G, comp.dist = NULL, comp.param = NULL)
{
  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
    stop("Known components must be specified in the two admixture models.")
  }

  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
  }

  ## Extracts the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
  #  comp_hessian <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_hessian <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_hessian)) assign(x = names(comp_hessian)[i], value = comp_hessian[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_hessian)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                                                                                            paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_hessian)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Differentiates when the 4 components were specified (simulations) and real-life cases (only known components are given).
  ## First, real-life case:
  if (all(unlist(sapply(comp.param, is.na)[c(1,3)]))) {
    if (any(exp.comp.dist == "stepfun")) {
      G1.fun <- eval(parse(text = expr[2]))
      G2.fun <- eval(parse(text = expr[4]))
      G1 <- function(z) G1.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      G1 <- function(z) { eval(parse(text = expr[2])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  } else {
    ## case of simulated data:
    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      G1.fun <- eval(parse(text = expr[2]))
      F2.fun <- eval(parse(text = expr[3]))
      G2.fun <- eval(parse(text = expr[4]))
      F1 <- function(z) F1.fun(z)
      G1 <- function(z) G1.fun(z)
      F2 <- function(z) F2.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
      F2 <- function(z) { eval(parse(text = expr[3])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  }

  ## Calcul des fonctions de repartition (empiriques provenant de l'echantillon X, ou exactes) au point 'z':
  if ( !is.null(known.p) ) {
    L1 <- function(z) { known.p[1] * F1(z) + (1-known.p[1]) * G1(z) }
    L2 <- function(z) { known.p[2] * F2(z) + (1-known.p[2]) * G2(z) }
  } else {
    L1 <- stats::ecdf(sample1)  # hat(L1)
    L2 <- stats::ecdf(sample2)  # hat(L2)
  }

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  ## Loi de G et support d'integration:
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  ## Differentiates cases where G1 = G2 and G1 != G2 :
  if (G1equalG2) {

    stopifnot(length(par) == 1)
    integrand <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/fixed.p1) * (L1(z) - (1-fixed.p1) * G1(z))
      F.Y.dataPoint <- (1/par) * (L2(z) - (1-par) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.dataPoint <- (1/par^2) * (L2(z) - G2(z))
      hessianD.dataPoint <- -(2/par^3) * (L2(z) - G2(z))
      (hessianD.dataPoint * D + gradD.dataPoint^2) * densite.G.dataPoint
    }

    if (support == "continuous") {
      err.tol <- 1e-03
      hessian.contr <- try(suppressWarnings(
        2 * stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L,
                             rel.tol = err.tol)$value), silent = TRUE)
      while ( (class(hessian.contr) == "character") & (err.tol < 1) ) {
        err.tol <- err.tol * 10
        hessian.contr <- try(suppressWarnings(
          2 * stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L,
                               rel.tol = err.tol)$value), silent = TRUE)
      }
    } else {
      hessian.contr <- 2 * sum( unlist(sapply(supp.integration, integrand, par)) )
    }

  } else {

    integrand.1.1 <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        #densite.G.dataPoint <- (table(c(sample1,sample2)) / sum(table(c(sample1,sample2))))[which(names(table(c(sample1,sample2))) == z)]
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
      F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.pX.dataPoint <- (1/par[1]^2) * (G1(z) - L1(z))
      hessianD.dataPoint <- -(2/par[1]^3) * (G1(z) - L1(z))
      (gradD.pX.dataPoint^2 + D * hessianD.dataPoint) * densite.G.dataPoint
    }
    integrand.1.2 <- integrand.2.1 <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      gradD.pX.dataPoint <- (1/par[1]^2) * (G1(z) - L1(z))
      gradD.pY.dataPoint <- -(1/par[2]^2) * (G2(z) - L2(z))
      (gradD.pX.dataPoint * gradD.pY.dataPoint) * densite.G.dataPoint
    }
    integrand.2.2 <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
      F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.pY.dataPoint <- -(1/par[2]^2) * (G2(z) - L2(z))
      hessianD.dataPoint <- (2/par[1]^3) * (G2(z) - L2(z))
      (gradD.pY.dataPoint^2 + D * hessianD.dataPoint) * densite.G.dataPoint
    }

    hessian.contr <- matrix(NA, nrow = 2, ncol = 2)
    if (support == "continuous") {

      err.tol <- 1e-03
      term1 <- try(2 * stats::integrate(integrand.1.1, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term1) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term1 <- try(2 * stats::integrate(integrand.1.1, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      }
      hessian.contr[1,1] <- term1

      err.tol <- 1e-03
      term2 <- try(2 * stats::integrate(integrand.1.2, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term2) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term2 <- try(2 * stats::integrate(integrand.1.2, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      }
      hessian.contr[1,2] <- term2
      hessian.contr[2,1] <- hessian.contr[1,2]

      err.tol <- 1e-03
      term3 <- try(2 * stats::integrate(integrand.2.2, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term3) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term3 <- try(2 * stats::integrate(integrand.2.2, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      }
      hessian.contr[2,2] <- term3

    } else {

      hessian.contr[1,1] <- 2 * sum( unlist(sapply(supp.integration, integrand.1.1, par)) )
      hessian.contr[1,2] <- 2 * sum( unlist(sapply(supp.integration, integrand.1.2, par)) )
      hessian.contr[2,1] <- hessian.contr[1,2]
      hessian.contr[2,2] <- 2 * sum( unlist(sapply(supp.integration, integrand.2.2, par)) )
    }
  }

  return(hessian.contr)
}


#' Gap b/w estimated ECDF of unknown components from two admixture models
#'
#' Compute the 'gap' between two unknown empirical cumulative distribution functions (ECDF) at some given point, in admixture models
#' with probability distribution function (pdf) given by l where l = p*f + (1-p)*g.
#' Uses the inversion method to do so, i.e. f = (1/p) (l - (1-p)*g), where g represents the known component of the admixture
#' model and p is the unknown proportion of the unknown component. Therefore, compute:
#'    D(z,L1,L2,p1,p2) = F1(z,L1,p1) - F2(z,L2,p2)
#' This measure should be integrated over some domain to compute the global discrepancy in the IBM approach, see further information in 'Details' below.
#'
#' @param z the point at which the difference between both unknown (estimated) component distributions is computed.
#' @param par Numeric vector with two elements, corresponding to the weights of the unknown component for the two admixture models.
#' @param fixed.p1 (optional, NULL by default) Arbitrary value chosen by the user for the component weight related to the unknown
#'                  component of the first admixture model. Only useful for optimization when the known components of the two
#'                  models are identical (G1=G2, leading to unidimensional optimization).
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return the gap evaluated at the specified point between the unknown components of the two observed samples.
#'
#' @examples
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' IBM_gap(z = 2.8, par = c(0.3,0.6), fixed.p1 = NULL, sample1 = sample1[['mixt.data']],
#'         sample2 = sample2[['mixt.data']], comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

IBM_gap <- function(z, par, fixed.p1 = NULL, sample1, sample2, comp.dist, comp.param)
{
  if (is.null(fixed.p1)) stopifnot(length(par) == 2)
  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
    stop("Known components must be specified in the two admixture models.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
  }

  ## Extracts the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
  #  comp_IBM_gap <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_IBM_gap <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_IBM_gap)) assign(x = names(comp_IBM_gap)[i], value = comp_IBM_gap[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_IBM_gap)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                                                                                            paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_IBM_gap)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  if (any(exp.comp.dist == "stepfun")) {
    G1.fun <- eval(parse(text = expr[2]))
    G2.fun <- eval(parse(text = expr[4]))
    G1 <- function(z) G1.fun(z)
    G2 <- function(z) G2.fun(z)
  } else {
    G1 <- function(z) { eval(parse(text = expr[2])) }
    G2 <- function(z) { eval(parse(text = expr[4])) }
  }
  ## Calcul des fonctions de repartition empiriques utiles provenant des echantillons X et Y:
  L1 <- stats::ecdf(sample1)								    # L1
  L2 <- stats::ecdf(sample2)								    # L2

  ## Calcul explicite de la difference:
  if (is.null(fixed.p1)) {
    F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
    F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
  } else {
    F.X.dataPoint <- (1/fixed.p1) * (L1(z) - (1-fixed.p1) * G1(z))
    F.Y.dataPoint <- (1/par) * (L2(z) - (1-par) * G2(z))
  }
  difference <- F.X.dataPoint - F.Y.dataPoint

  return(difference)
}


#' Gap b/w CDF of unknown components from 2 admixture models
#'
#' Computes the difference of two cumulative distribution functions (CDF) corresponding to the two 'unknown' component
#' distributions coming from 2 admixture models, at some given point. Useful in case of simulation studies, since all parameters are known.
#' Remind that each admixture model has probability distribution function (pdf) given by l where  l = p*f + (1-p)*g.
#' Uses the inversion technique, i.e. f = (1/p) (l - (1-p)g), where g represents the known component of the admixture
#' model and p is the proportion of the unknown component. This difference must be integrated over some domain to compute the
#' global discrepancy, as introduced in the paper presenting the IBM approach (see 'Details' below).
#'
#' @param z Point at which the difference between the unknown component distributions of the two considered admixture models is computed.
#' @param par Numeric vector with two elements, corresponding to the two parameter values at which to compute the gap. In practice
#'            the component weights for the two admixture models.
#' @param known.p Numeric vector with two elements, the known (true) mixture weights.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. No unknown elements permitted.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1='rnorm', g1='norm', f2='rnorm', g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. No unknown elements permitted. For instance, 'comp.param' could be specified
#'                   as follows: : list(f1 = list(mean=2,sd=0.3), g1 = list(mean=0,sd=1), f2 = list(mean=2,sd=0.3), g2 = list(mean=3,sd=1.1)).
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return The gap between F1 and F2 (unknown components of the two admixture models), evaluated at the specified point.
#'
#' @examples
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' IBM_theoretical_gap(z = 2.8, par = c(0.3,0.6), known.p = c(0.5,0.5),
#'                     comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

IBM_theoretical_gap <- function(z, par, known.p = c(0.5,0.5), comp.dist, comp.param)
{
  if (is.null(known.p)) stop("In 'IBM_theoretical_gap': argument 'known.p' must be specified.")
  if ( (length(comp.dist) != 4) | (length(comp.param) != 4) | is.null(known.p) )
  { stop("The theoretical version of the gap can only be computed when all mixture components / weights are specified") }

  ## Extracts the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
  comp_gap <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  #  comp_gap <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  for (i in 1:length(comp_gap)) assign(x = names(comp_gap)[i], value = comp_gap[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_gap)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                                                                                        paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_gap)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  if (any(exp.comp.dist == "stepfun")) {
    F1.fun <- eval(parse(text=expr[1]))
    G1.fun <- eval(parse(text=expr[2]))
    F2.fun <- eval(parse(text=expr[3]))
    G2.fun <- eval(parse(text=expr[4]))
    F1 <- function(z) F1.fun(z)
    G1 <- function(z) G1.fun(z)
    F2 <- function(z) F2.fun(z)
    G2 <- function(z) G2.fun(z)
  } else {
    F1 <- function(z) { eval(parse(text=expr[1])) }
    G1 <- function(z) { eval(parse(text=expr[2])) }
    F2 <- function(z) { eval(parse(text=expr[3])) }
    G2 <- function(z) { eval(parse(text=expr[4])) }
  }

  L1.theo <- function(z) { known.p[1] * F1(z) + (1-known.p[1]) * G1(z) }
  L2.theo <- function(z) { known.p[2] * F2(z) + (1-known.p[2]) * G2(z) }

  ## Calcul explicite de la difference:
  if (length(par) == 2) {
    F.X.dataPoint <- (1/par[1]) * (L1.theo(z) - (1-par[1]) * G1(z))
    F.Y.dataPoint <- (1/par[2]) * (L2.theo(z) - (1-par[2]) * G2(z))
  } else {
    F.X.dataPoint <- (1/known.p[1]) * (L1.theo(z) - (1-known.p[1]) * G1(z))
    F.Y.dataPoint <- (1/par) * (L2.theo(z) - (1-par) * G2(z))
  }
  difference <- F.X.dataPoint - F.Y.dataPoint

  return(difference)
}
