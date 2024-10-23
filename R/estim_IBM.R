#' Estimates weights of unknown components from 2 admixtures using IBM
#'
#' Estimation of the component weights from the Inversion - Best Matching (IBM) method, related to two admixture models
#' with respective probability density function (pdf) l1 and l2, such that:
#'    l1 = p1*f1 + (1-p1)*g1 and l2 = p2*f2 + (1-p2)*g2, where g1 and g2 are the known component densities.
#' For further details about IBM approach, see 'Details' below.
#'
#' @param samples (List) List of the two considered samples.
#' @param admixMod (List) List of objects of class 'admix_model', one for each sample.
#' @param n.integ Number of data points generated for the distribution on which to integrate.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return An object of class 'estim_IBM', containing 7 attributes: 1) the number of samples under study; 2) the sizes of samples;
#'         3) the information about mixture components (distributions and parameters) for each sample; 4) the estimation
#'         method (Inversion Best Matching here, see the given reference); 5) the estimated mixing proportions (weights of the
#'         unknown component distributions in each sample); 6) the arbitrary value of the mixing weight in the first admixture sample
#'         (in case of equal known components, see the given reference); 7) the support of integration that was used in the computations.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1500, weight = 0.5,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' mixt2 <- twoComp_mixt(n = 2000, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 5, "sd" = 2)))
#' data2 <- getmixtData(mixt2)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#'
#' ## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
#' estim_IBM(samples = list(data1,data2),
#'           admixMod = list(admixMod1,admixMod2), n.integ = 1000)
#'
#' ## Example 2: multinomial distribution:
#' mixt1 <- twoComp_mixt(n = 1500, weight = 0.8,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size" = 1, "prob" = c(0.2,0.3,0.5)),
#'                                         list("size" = 1, "prob" = c(0.1,0.6,0.3))))
#' data1 <- getmixtData(mixt1)
#' mixt2 <- twoComp_mixt(n = 2000, weight = 0.3,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size" = 1, "prob" = c(0.2,0.3,0.5)),
#'                                         list("size" = 1, "prob" = c(0.7,0.1,0.2))))
#' data2 <- getmixtData(mixt2)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#'
#' ## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
#' estim_IBM(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

estim_IBM <- function(samples, admixMod, n.integ = 1000)
{
  warning(" IBM estimators of two unknown proportions are reliable only if the two
    corresponding unknown component distributions have been tested equal (see 'admix_test()').
    Furthermore, when both the known and unknown component distributions of the mixture
    models are identical, the IBM approach provides an estimation of the ratio of the
    actual mixing weights rather than an estimation of the unknown weights themselves.\n")

  ##------- Defines the support for integration by simulation --------##
  ## Allows to integrate the gap between F1 and F2 in the contrast computation. span(G) must contain span(X) and span(Y).
  ## Ideally, this distribution should put more weight to locations where differences between F1 and F2 are expected.
  support <- detect_support_type(samples[[1]], samples[[2]])
  if (support == "Continuous") {
    G <- stats::runif(n.integ, min = min(c(samples[[1]], samples[[2]])), max = max(c(samples[[1]], samples[[2]])))
  } else {
    G <- unique(sort(c(unique(samples[[1]]), unique(samples[[2]]))))
  }

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])
  if (G1equalG2) {
    ## Leads to a one-dimensional optimization because known components of mixture distributions are the same.
    ## Set arbitrarily the proportion to 1/2, could be any other value belonging to ]0,1[:
    fixed.p.X <- 0.5
    par.init <- 0.5
  } else {
    ## Otherwise classical optimization to get the two component weights:
    fixed.p.X <- NULL
    par.init <- c(0.5,0.5)
  }

  ##------- Solution search by optimization --------##
  ## Use method 'L-BFGS-B' to optimize under constraints (on parameter space), and use 'try' to catch errors coming from divergent integrals (numerical issues):
  expr1_NM <- expression(stats::optim(par = par.init, fn = IBM_empirical_contrast, gr = NULL, samples = samples, admixMod = admixMod,
                                      G = G, fixed.p.X = fixed.p.X, method = "Nelder-Mead", control = list(trace = 0, maxit = 10000)))
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
    #message("In 'IBM_estim': optimization algorithm was changed (in 'optim') from 'Nelder-Mead' to 'BFGS' to avoid the solution to explose.")
    expr1_BFGS <- expression(stats::optim(par = par.init, fn = IBM_empirical_contrast, gr = NULL, samples = samples, admixMod = admixMod,
                                          G = G, fixed.p.X = fixed.p.X, method = "L-BFGS-B", lower = c(0.001,0.001), upper = c(5,5),
                                          control = list(trace = 0, maxit = 10000)))
    sol <- try(suppressWarnings(eval(expr1_BFGS)), silent = TRUE)
    count_error <- 0
    while (inherits(x = sol, what = "try-error", which = FALSE) & (count_error < 7)) {
      sol <- NULL
      sol <- try(suppressWarnings(eval(expr1_BFGS)), silent = TRUE)
      count_error <- count_error + 1
    }
    if (inherits(x = sol, what = "try-error", which = FALSE)) {
      #message("In 'IBM_estim': impossible to estimate the component weights with BFGS method. Switch back to Nelder-Mead algorithm to obtain a solution")
      sol <- try(suppressWarnings(eval(expr1_NM)), silent = TRUE)
    }
  }

  if (inherits(x = sol, what = "try-error", which = FALSE)) {
    stop("In 'IBM_estim': whatever the optimization algorithm, the solution has not been found.")
  } else {
    estim.weights <- sol[['par']]
  }

  res <- list(
    n_populations = length(samples),
    population_sizes = sapply(X = samples, FUN = length),
    admixture_models = admixMod,
    estimation_method = "Inversion Best Matching (IBM)",
    estimated_mixing_weights = estim.weights,
    p.X.fixed = fixed.p.X,
    integ.supp = sort(G)
    )
  class(res) <- c("estim_IBM", "admix_estim")
  res$call <- match.call()
  return(res)
}


#' Print method for objects of class 'estim_IBM'
#'
#' Print the results stored in an object of class 'estim_IBM'.
#'
#' @param x An object of class 'estim_IBM'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.estim_IBM <- function(x, ...)
{
  cat("Call:")
  print(x$call)
  cat("\n")
  cat("Estimated mixing proportion (of the unknown component) in the 1st sample: ", x$estimated_mixing_weights[1], "\n")
  cat("Estimated mixing proportion (of the unknown component) in the 2nd sample: ", x$estimated_mixing_weights[2], "\n")
  if (!is.null(x$p.X.fixed))
    cat("\nFixed value for the mixing weight in the 1st sample (case of identical known components in the 2 samples):", x$p.X.fixed, "\n\n")
}


#' Summary method for objects 'estim_IBM'
#'
#' Summarizes the results stored in an object of class 'estim_IBM'.
#'
#' @param object An object of class 'estim_IBM'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.estim_IBM <- function(object, ...)
{
  cat("Call:")
  print(object$call)
  cat("\n")
  cat("------- Samples -------\n")
  cat("Number of samples: ", object$n_populations, "\n")
  cat("Sample sizes: ", object$population_sizes, "\n")
  for (k in 1:object$n_populations) {
    cat("-> Distribution and parameters of the known component \n for admixture model #", k, ": ", sep="")
    cat(paste(sapply(object$admixture_models[[k]], "[[", "known")[1:2], collapse = " - "))
    cat("\n")
  }
  cat("\n------- Estimation results -------\n")
  cat("Estimated mixing proportion (of the unknown component) in the 1st sample: ", object$estimated_mixing_weights[1], "\n")
  cat("Estimated mixing proportion (of the unknown component) in the 2nd sample: ", object$estimated_mixing_weights[2], "\n\n")
  if (!is.null(object$p.X.fixed))
    cat("\nFixed value for the mixing weight in the 1st sample (case of identical known components in the 2 samples):", object$p.X.fixed, "\n\n")
  cat("------- Support -------\n")
  cat("Integration support: ", paste(utils::head(object$integ.supp,3), collapse=" "), "...",
      paste(utils::tail(object$integ.supp,3), collapse = " "), "\n\n", sep="")
}


#' Empirical contrast in the IBM method
#'
#' Computes the empirical version of the contrast in the Inversion - Best Matching (IBM) method, to be minimized in the optimization process.
#' For further details about the contrast definition, see 'Details' below.
#'
#' @param par Numeric vector with two elements, corresponding to the two parameter values at which to compute the contrast. In practice
#'            the component weights for the two admixture models.
#' @param samples (List) List of the two considered samples.
#' @param admixMod (List) List of objects of class 'admix_model', one for each sample.
#' @param G Distribution on which to integrate when calculating the contrast.
#' @param fixed.p.X (default to NULL) Arbitrary value chosen by the user for the component weight related to the unknown component of
#'                  the first admixture model. Only useful for optimization when the known components of the two models are identical
#'                  (G1=G2, leading to unidimensional optimization).
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return The empirical contrast value evaluated at parameter values.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1500, weight = 0.5,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' mixt2 <- twoComp_mixt(n = 2000, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 5, "sd" = 2)))
#' data2 <- getmixtData(mixt2)
#'
#' ## Create the distribution on which the contrast will be integrated:
#' G <- stats::rnorm(n = 1000, mean = sample(c(data1, data2), size = 1000, replace = TRUE),
#'                   sd = density(c(data1, data2))$bw)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#'
#' ## Compute the empirical contrast at parameters (p1,p2) = (0.2,0.7) in a real-life setting:
#' IBM_empirical_contrast(par = c(0.2,0.7), samples = list(data1, data2),
#'                        admixMod = list(admixMod1, admixMod2), G=G, fixed.p.X = NULL)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

IBM_empirical_contrast <- function(par, samples, admixMod, G, fixed.p.X = NULL)
{
  ## Extract the information on component distributions:
  knownCDF_comp.dist <- paste0("p", unlist(sapply(admixMod, '[[', 'comp.dist')["known", ]))
  if (any(knownCDF_comp.dist == "pmultinom")) knownCDF_comp.dist[which(knownCDF_comp.dist == "pmultinom")] <- "stepfun"
  comp_emp <- sapply(X = knownCDF_comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])

  ## Create the expression involved in future assessments of the CDF:
#  make.expr.step <- function(i) paste(names(comp_emp)[i],"(x = 1:", length(admixMod[[i]]$comp.param$known[[2]]$prob),
#      paste(", y = ", paste("cumsum(c(0,", paste(admixMod[[i]]$comp.param$known[[2]]$prob, collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr.step <- function(i) paste(names(comp_emp)[i],"(x = 1:", length(admixMod[[i]]$comp.param$known$prob),
                                      paste(", y = ", paste("cumsum(c(0,", paste(admixMod[[i]]$comp.param$known$prob, collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_emp)[i],"(z,", paste(names(admixMod[[i]]$comp.param$known),
                                    "=", admixMod[[i]]$comp.param$known, sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(knownCDF_comp.dist))
  expr[which(knownCDF_comp.dist == "stepfun")] <- sapply(which(knownCDF_comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  if (any(knownCDF_comp.dist == "stepfun")) {
    G1.fun <- eval(parse(text = expr[1]))
    G2.fun <- eval(parse(text = expr[2]))
    G1 <- function(z) G1.fun(z)
    G2 <- function(z) G2.fun(z)
  } else {
    G1 <- function(z) { eval(parse(text = expr[1])) }
    G2 <- function(z) { eval(parse(text = expr[2])) }
  }

  ## Empirical cumulative distribution function from the two observed samples:
  L1 <- stats::ecdf(samples[[1]])						# hat(L1)
  L2 <- stats::ecdf(samples[[2]])						# hat(L2)

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])

  ##------- Computes the empirical contrast --------##
  support <- detect_support_type(samples[[1]], samples[[2]])
  if (support == "Continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  if (!G1equalG2) {
    integrand <- function(z, par) {
      if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else { densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]]))) }
      F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
      F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
      weighted.difference <- (F.X.dataPoint - F.Y.dataPoint)^2 * densite.G.dataPoint
      weighted.difference
    }
  } else {   # for one-dimensional optimization
    stopifnot("Specify argument 'fixed.p.X' since known components are identical " = !is.null(fixed.p.X))
    integrand <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else { densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]]))) }
      F.X.dataPoint <- (1/fixed.p.X) * (L1(z) - (1-fixed.p.X) * G1(z))
      F.Y.dataPoint <- (1/par) * (L2(z) - (1-par) * G2(z))
      weighted.difference <- (F.X.dataPoint - F.Y.dataPoint)^2 * densite.G.dataPoint
      weighted.difference
    }
  }

  if (support == "Continuous") {
    res <- stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2],
                            par, subdivisions = 10000L, rel.tol = 1e-03)$value
  } else {
    res <- sum( unlist(sapply(supp.integration, integrand, par)) )
  }
  return(res)
}
