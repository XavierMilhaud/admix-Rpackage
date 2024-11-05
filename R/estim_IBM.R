#' Estimates weights of unknown components from 2 admixtures using IBM
#'
#' Estimation of the component weights from the Inversion - Best Matching (IBM) method, related to two admixture models
#' with respective probability density function (pdf) l1 and l2, such that:
#'    l1 = p1*f1 + (1-p1)*g1 and l2 = p2*f2 + (1-p2)*g2, where g1 and g2 are the known component densities.
#' For further details about IBM approach, see 'Details' below.
#'
#' @param samples A list of the two considered samples.
#' @param admixMod A list of two objects of class \link[admix]{admix_model}, one for each sample.
#' @param n.integ Number of data points generated for the distribution on which to integrate.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return An object of class \link[admix]{estim_IBM}, containing 7 attributes: 1) the number of samples under study; 2) the sizes of samples;
#'         3) the information about mixture components (distributions and parameters) for each sample; 4) the estimation
#'         method (Inversion Best Matching here, see the given reference); 5) the estimated mixing proportions (weights of the
#'         unknown component distributions in each sample); 6) the arbitrary value of the mixing weight in the first admixture sample
#'         (in case of equal known components, see the given reference); 7) the support of integration that was used in the computations.
#'
#' @examples
#' ## Continuous support: simulate mixture data.
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
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
#' estim_IBM(samples = list(data1,data2),
#'           admixMod = list(admixMod1,admixMod2), n.integ = 1000)
#'
#' ## Example 2: multinomial distribution (discrete case).
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
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
#' estim_IBM(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

estim_IBM <- function(samples, admixMod, n.integ = 1000)
{
  stopifnot("Wrong number of samples... Must be 2!" = length(samples) == 2)

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
    equal.knownComp = G1equalG2,
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
  if (x$equal.knownComp) {
    cat("Fixed proportion (of the unknown component) in the 1st sample (since equal unknown components): ", x$p.X.fixed, "\n")
    cat("Estimated mixing proportion (of the unknown component) in the 2nd sample: ", x$estimated_mixing_weights, "\n")
  } else {
    cat("Estimated mixing proportion (of the unknown component) in the 1st sample: ", x$estimated_mixing_weights[1], "\n")
    cat("Estimated mixing proportion (of the unknown component) in the 2nd sample: ", x$estimated_mixing_weights[2], "\n")
  }
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
  cat("Are the known component distributions equal? ", object$equal.knownComp,"\n")
  cat("\n------- Estimation results -------\n")
  if (object$equal.knownComp) {
    cat("Fixed proportion (of the unknown component) in the 1st sample (since equal unknown components): ", object$p.X.fixed, "\n")
    cat("Estimated mixing proportion (of the unknown component) in the 2nd sample: ", object$estimated_mixing_weights, "\n")
  } else {
    cat("Estimated mixing proportion (of the unknown component) in the 1st sample: ", object$estimated_mixing_weights[1], "\n")
    cat("Estimated mixing proportion (of the unknown component) in the 2nd sample: ", object$estimated_mixing_weights[2], "\n")
  }
  cat("\n------- Support -------\n")
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



#' Estimates the covariance matrix of the gaussian vector (IBM)
#'
#' Nonparametric estimation of the covariance matrix of the gaussian vector at point 'z', considering the use of Inversion - Best
#' Matching (IBM) method to estimate the model parameters in two-sample admixture models.
#' Recall that the two admixture models have respective probability density functions (pdf) l1 and l2, such that:
#'   l1 = p1*f1 + (1-p1)*g1 and l2 = p2*f2 + (1-p2)*g2, where g1 and g2 are the known component densities.
#' Further information for the IBM approach are given in 'Details' below.
#'
#' @param x Time point at which the 1st (related to the 1st parameter) underlying empirical process is looked through.
#' @param y Time point at which the 2nd (related to the 2nd parameter) underlying empirical process is looked through.
#' @param IBMestim.obj An object of class 'estim_IBM'.
#' @param samples (List) List of the two considered samples.
#' @param admixMod (List) List of objects of class 'admix_model', one for each sample.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return The estimated variance-covariance matrix of the gaussian vector Z = (hat(p1),(hat(p2),Dn(z)), at location '(x,y)'.
#'
#' @examples
#' \donttest{
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 2000, weight = 0.2,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 2, "sd" = 3),
#'                                         list("mean" = -2, "sd" = 1)))
#' plot(mixt1, xlim = c(0,30), ylim = c(0,0.15))
#' data1 <- getmixtData(mixt1)
#' mixt2 <- twoComp_mixt(n = 1500, weight = 0.5,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 2, "sd" = 3),
#'                                         list("mean" = 6, "sd" = 1)))
#' plot(mixt2, add.plot = TRUE, xlim = c(0,30), ylim = c(0,0.15), col = "blue")
#' data2 <- getmixtData(mixt2)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#'
#' ## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
#' est <- estim_IBM(samples = list(data1,data2),
#'                  admixMod = list(admixMod1,admixMod2), n.integ = 1000)
#'
#' IBM_estimVarCov_gaussVect(x = mean(data1), y = mean(data2), IBMestim.obj = est,
#'                           samples=list(data1,data2), admixMod = list(admixMod1,admixMod2))
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

IBM_estimVarCov_gaussVect <- function(x, y, IBMestim.obj, samples, admixMod)
{
  estimators <- IBMestim.obj$estimated_mixing_weights
  ## Location at which the variance-covariance matrix is evaluated:
  varCov.gaussianVect <- IBM_mat_Sigma(x = x, y = y, par = estimators, samples = samples, admixMod = admixMod, G = IBMestim.obj$integ.supp)
  ## Adjusts by the normalization factor to get the distribution of the gaussian vector Z=(hat(p1), hat(p2), Dn(z)):
  normalization.factor <- IBM_normalization_term(z = x, IBMestim.obj = IBMestim.obj, samples = samples, admixMod = admixMod)
  varCov.transfo_gaussVect <- normalization.factor %*% varCov.gaussianVect %*% t(normalization.factor)

  return(varCov.transfo_gaussVect)
}


## Normalization of the variance-covariance by matrix 'M' in the paper
IBM_normalization_term <- function(z, IBMestim.obj, samples, admixMod)
{
  ## Adjusts by the normalization factor to get the distribution of the gaussian vector Z=(hat(p1), hat(p2), Dn(z)):
  L <- IBM_mat_L(z = z, par = IBMestim.obj$estimated_mixing_weights, IBMestim.obj = IBMestim.obj, samples = samples,
                 admixMod = admixMod)
  inv_J <- solve(IBM_mat_J(par = IBMestim.obj$estimated_mixing_weights, IBMestim.obj = IBMestim.obj, samples = samples,
                           admixMod = admixMod, G = IBMestim.obj$integ.supp))

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])
  n1 <- length(samples[[1]])
  n2 <- length(samples[[2]])
  xi <- 1 / sqrt(max(n1,n2)/min(n1,n2))

  if (G1equalG2) {
    if (n1 <= n2) {
      C <- matrix(c(-1,0,-xi,0,
                    0,1,0,0,
                    0,0,0,1), nrow = 3, ncol = 4, byrow = T)
    } else {
      C <- matrix(c(-xi,0,-1,0,
                    0,1,0,0,
                    0,0,0,1), nrow = 3, ncol = 4, byrow = T)
    }
  } else {
    if (n1 <= n2) {
      C <- matrix(c(-1,0,0,0,-xi,0,
                    0,-1,0,-xi,0,0,
                    0,0,1,0,0,0,
                    0,0,0,0,0,1), nrow = 4, ncol = 6, byrow = T)
    } else {
      C <- matrix(c(-xi,0,0,0,-1,0,
                    0,-xi,0,-1,0,0,
                    0,0,1,0,0,0,
                    0,0,0,0,0,1), nrow = 4, ncol = 6, byrow = T)
    }
  }

  return(L %*% inv_J %*% C)
}


## Matrice J de l'article (normalisation par la variance):
IBM_mat_J <- function(par, IBMestim.obj, samples, admixMod, G)
{
  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])

  if (G1equalG2) {
    stopifnot("Error in 'IBM_mat_J: 'p.X.fixed' must have a value" = !is.null(IBMestim.obj$p.X.fixed))
    J <- diag(1, nrow = 3, ncol = 3)
    J[1,1] <- IBM_hessian_contrast(par = par, IBMestim.obj = IBMestim.obj, samples = samples, admixMod = admixMod, G = G)
  } else {
    J <- diag(1, nrow = 4, ncol = 4)
    J[1:2,1:2] <- IBM_hessian_contrast(par = par, IBMestim.obj = IBMestim.obj, samples = samples, admixMod = admixMod, G = G)
  }
  return(J)
}

## Matrice L de l'article (matrice muette pour les estimateurs de p1 et p2):
IBM_mat_L <- function(z, par, IBMestim.obj, samples, admixMod)
{
  ## Extract the information on component distributions:
  knownCDF_comp.dist <- paste0("p", unlist(sapply(admixMod, '[[', 'comp.dist')["known", ]))
  if (any(knownCDF_comp.dist == "pmultinom")) knownCDF_comp.dist[which(knownCDF_comp.dist == "pmultinom")] <- "stepfun"
  comp_emp <- sapply(X = knownCDF_comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])

  ## Create the expression involved in future assessments of the CDF:
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

  ## Cumulative distribution functions (empirical) at point 'z':
  L1 <- stats::ecdf(samples[[1]])
  L2 <- stats::ecdf(samples[[1]])

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])
  n1 <- length(samples[[1]])
  n2 <- length(samples[[2]])
  xi <- 1 / sqrt(max(n1,n2)/min(n1,n2))

  if (G1equalG2) {
    stopifnot("Error in 'IBM_mat_L: 'p.X.fixed' must have a value" = !is.null(IBMestim.obj$p.X.fixed))
    L <- matrix(0, nrow = 2, ncol = 3)
    if (n1 <= n2) {
      L[1,1] <- 1
      L[2,1] <- (1/par^2) * (L2(z) - G2(z))
      L[2,2] <- 1 / IBMestim.obj$p.X.fixed
      L[2,3] <- -xi / par
    } else {
      L[1,1] <- 1
      L[2,1] <- (1/par^2) * (L2(z) - G2(z))
      L[2,2] <- xi / IBMestim.obj$p.X.fixed
      L[2,3] <- -1 / par
    }
  } else {
    L <- matrix(0, nrow = 3, ncol = 4)
    if (n1 <= n2) {
      L[1,1] <- L[2,2] <- 1
      L[3,1] <- -(1/par[1]^2) * (L1(z) - G1(z))
      L[3,2] <- (1/par[2]^2) * (L2(z) - G2(z))
      L[3,3] <- 1 / par[1]
      L[3,4] <- -xi / par[2]
    } else {
      L[1,1] <- L[2,2] <- 1
      L[3,1] <- -(1/par[1]^2) * (L1(z) - G1(z))
      L[3,2] <- (1/par[2]^2) * (L2(z) - G2(z))
      L[3,3] <- xi / par[1]
      L[3,4] <- -1 / par[2]
    }
  }

  return(L)
}


## Variance-covariance matrix of the Gaussian random vector. Arguments are defined as follows:
## This function returns the non-normalized variance-covariance matrix of the gaussian random vector.
IBM_mat_Sigma <- function(x, y, par, IBMestim.obj, samples, admixMod, G)
{
  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])

  if (G1equalG2) {
    Sigma <- matrix(0, nrow = 4, ncol = 4)
    ## Upper block of the matrix :
    Sigma[1:2,1:2] <- IBM_Sigma1(x = x, y = y, par = par, IBMestim.obj = IBMestim.obj, samples = samples, admixMod = admixMod, G = G)
    ## Lower block of the matrix :
    Sigma[3:4,3:4] <- IBM_Sigma2(x = x, y = y, par = par, IBMestim.obj = IBMestim.obj, samples = samples, admixMod = admixMod, G = G)
  } else {
    Sigma <- matrix(0, nrow = 6, ncol = 6)
    ## Upper block of the matrix :
    Sigma[1:3,1:3] <- IBM_Sigma1(x = x, y = y, par = par, IBMestim.obj = IBMestim.obj, samples = samples, admixMod = admixMod, G = G)
    ## Lower block of the matrix :
    Sigma[4:6,4:6] <- IBM_Sigma2(x = x, y = y, par = par, IBMestim.obj = IBMestim.obj, samples = samples, admixMod = admixMod, G = G)
  }
  return(Sigma)
}


## Upper block of the variance-covariance matrix:
IBM_Sigma1 <- function(x, y, par, IBMestim.obj, samples, admixMod, G)
{
  ## Extract the information on component distributions:
  knownCDF_comp.dist <- paste0("p", unlist(sapply(admixMod, '[[', 'comp.dist')["known", ]))
  if (any(knownCDF_comp.dist == "pmultinom")) knownCDF_comp.dist[which(knownCDF_comp.dist == "pmultinom")] <- "stepfun"
  comp_emp <- sapply(X = knownCDF_comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])

  ## Create the expression involved in future assessments of the CDF:
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

  ## Empirical cumulative distribution functions:
  L1 <- stats::ecdf(samples[[1]])
  L2 <- stats::ecdf(samples[[2]])

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])

  ## Loi de G et support d'integration :
  support <- detect_support_type(samples[[1]], samples[[2]])
  if (support == "Continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  ## Numerical computation of each term included in the variance-covariance matrix of the gaussian vector:
  if (G1equalG2) {
    stopifnot(length(par) == 1)
    psi1 <- function(z) 2*( ((2-par)/par^3) * G1(z) - (2/par^3)*L2(z) + (1/(par^2*IBMestim.obj$p.X.fixed))*L1(z) -
                              ((1-IBMestim.obj$p.X.fixed)/(par^2*IBMestim.obj$p.X.fixed)) * G1(z) )
    psi2 <- function(z) 2*( (1/(par^2*IBMestim.obj$p.X.fixed)) * (L2(z) - G1(z)) )
    integrand.sigma1.11_g1eqg2 <- function(x,y) {
      if (support == "Continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2(x) * psi2(y) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma1.12_g1eqg2 <- function(x) {
      if (support == "Continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2(x) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.x
    }
    integrand.sigma1.21_g1eqg2 <- function(y) {
      if (support == "Continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2(y) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.y
    }

    sigma <- matrix(NA, nrow = 2, ncol = 2)
    if (support == "Continuous") {
      sigma[1,1] <- pracma::integral2(Vectorize(integrand.sigma1.11_g1eqg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,2] <- stats::integrate(Vectorize(integrand.sigma1.12_g1eqg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,1] <- stats::integrate(Vectorize(integrand.sigma1.21_g1eqg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,2] <- estimVarCov_empProcess(x, y, samples[[1]])
    } else {
      sigma[1,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma1.11_g1eqg2, x, supp.integration)) )
      sigma[1,2] <- sum( unlist(sapply(supp.integration, integrand.sigma1.12_g1eqg2)) )
      sigma[2,1] <- sum( unlist(sapply(supp.integration, integrand.sigma1.21_g1eqg2)) )
      sigma[2,2] <- estimVarCov_empProcess(x, y, samples[[1]])
    }

  } else {

    psi1.1 <- function(z) 2*( ((2-par[1])/par[1]^3) * G1(z) - (2/par[1]^3)*L1(z) + (1/(par[1]^2*par[2]))*L2(z) - ((1-par[2])/(par[1]^2*par[2])) * G2(z) )
    psi1.2 <- function(z) 2*( (1/(par[1]^2*par[2])) * (L1(z) - G1(z)) )
    psi2.1 <- function(z) 2*( ((2-par[2])/par[2]^3) * G2(z) - (2/par[2]^3)*L2(z) + (1/(par[2]^2*par[1]))*L1(z) - ((1-par[1])/(par[2]^2*par[1])) * G1(z) )
    psi2.2 <- function(z) 2*( (1/(par[2]^2*par[1])) * (L2(z) - G2(z)) )
    integrand.sigma1.11_g1diffg2 <- function(x,y) {
      if (support == "Continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1.1(x) * psi1.1(y) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma1.12_g1diffg2 <- function(x,y) {
      if (support == "Continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1.1(x) * psi2.2(y) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma1.13_g1diffg2 <- function(x) {
      if (support == "Continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1.1(x) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.x
    }
    integrand.sigma1.22_g1diffg2 <- function(x,y) {
      if (support == "Continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2.2(x) * psi2.2(y) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma1.23_g1diffg2 <- function(x) {
      if (support == "Continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2.2(x) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.x
    }
    integrand.sigma1.31_g1diffg2 <- function(y) {
      if (support == "Continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1.1(y) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.y
    }
    integrand.sigma1.32_g1diffg2 <- function(y) {
      if (support == "Continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2.2(y) * estimVarCov_empProcess(x,y,samples[[1]]) * densite.G.dataPoint.y
    }

    sigma <- matrix(NA, nrow = 3, ncol = 3)
    if (support == "Continuous") {
      sigma[1,1] <- pracma::integral2(Vectorize(integrand.sigma1.11_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,2] <- sigma[2,1] <- pracma::integral2(Vectorize(integrand.sigma1.12_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,3] <- stats::integrate(Vectorize(integrand.sigma1.13_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,2] <- pracma::integral2(Vectorize(integrand.sigma1.22_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[2,3] <- stats::integrate(Vectorize(integrand.sigma1.23_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,1] <- stats::integrate(Vectorize(integrand.sigma1.31_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,2] <- stats::integrate(Vectorize(integrand.sigma1.32_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,3] <- estimVarCov_empProcess(x,y,samples[[1]])
    } else {
      sigma[1,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma1.11_g1diffg2, x, supp.integration)) )
      sigma[1,2] <- sigma[2,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma1.12_g1diffg2, x, supp.integration)) )
      sigma[1,3] <- sum( unlist(sapply(supp.integration, integrand.sigma1.13_g1diffg2)) )
      sigma[2,2] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma1.22_g1diffg2, x, supp.integration)) )
      sigma[2,3] <- sum( unlist(sapply(supp.integration, integrand.sigma1.23_g1diffg2)) )
      sigma[3,1] <- sum( unlist(sapply(supp.integration, integrand.sigma1.31_g1diffg2)) )
      sigma[3,2] <- sum( unlist(sapply(supp.integration, integrand.sigma1.32_g1diffg2)) )
      sigma[3,3] <- estimVarCov_empProcess(x,y,samples[[1]])
    }
  }

  return(sigma)
}


## Lower block:
IBM_Sigma2 <- function(x, y, par, IBMestim.obj, samples, admixMod, G)
{
  ## Extract the information on component distributions:
  knownCDF_comp.dist <- paste0("p", unlist(sapply(admixMod, '[[', 'comp.dist')["known", ]))
  if (any(knownCDF_comp.dist == "pmultinom")) knownCDF_comp.dist[which(knownCDF_comp.dist == "pmultinom")] <- "stepfun"
  comp_emp <- sapply(X = knownCDF_comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])

  ## Create the expression involved in future assessments of the CDF:
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

  ## Cumulative distribution function (either empirical from the two observed samples, or theoretical):
  L1 <- stats::ecdf(samples[[1]])
  L2 <- stats::ecdf(samples[[2]])

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])

  ## Loi de G et support d'integration:
  support <- detect_support_type(samples[[1]], samples[[2]])
  if (support == "Continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  ## Numerical computation of each term included in the variance-covariance matrix of the gaussian vector:
  if (G1equalG2) {

    stopifnot(length(par) == 1)
    psi1 <- function(z) 2*( ((2-par)/par^3) * G1(z) - (2/par^3)*L2(z) + (1/(par^2*IBMestim.obj$p.X.fixed))*L1(z) -
                              ((1-IBMestim.obj$p.X.fixed)/(par^2*IBMestim.obj$p.X.fixed)) * G1(z) )
    psi2 <- function(z) 2*( (1/(par^2*IBMestim.obj$p.X.fixed)) * (L2(z) - G1(z)) )

    integrand.sigma2.11_g1eqg2 <- function(x,y) {
      if (support == "Continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1(x) * psi1(y) * estimVarCov_empProcess(x,y, samples[[2]]) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma2.12_g1eqg2 <- function(x) {
      if (support == "Continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1(x) * estimVarCov_empProcess(x,y, samples[[2]]) * densite.G.dataPoint.x
    }
    integrand.sigma2.21_g1eqg2 <- function(y) {
      if (support == "Continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1(y) * estimVarCov_empProcess(x,y, samples[[2]]) * densite.G.dataPoint.y
    }

    sigma <- matrix(NA, nrow = 2, ncol = 2)
    if (support == "Continuous") {
      sigma[1,1] <- pracma::integral2(Vectorize(integrand.sigma2.11_g1eqg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,2] <- stats::integrate(Vectorize(integrand.sigma2.12_g1eqg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,1] <- stats::integrate(Vectorize(integrand.sigma2.21_g1eqg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,2] <- estimVarCov_empProcess(x, y, samples[[2]])
    } else {
      sigma[1,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma2.11_g1eqg2, x, supp.integration)) )
      sigma[1,2] <- sum( unlist(sapply(supp.integration, integrand.sigma2.12_g1eqg2)) )
      sigma[2,1] <- sum( unlist(sapply(supp.integration, integrand.sigma2.21_g1eqg2)) )
      sigma[2,2] <- estimVarCov_empProcess(x, y, samples[[2]])
    }

  } else {

    psi1.1 <- function(z) 2*( ((2-par[1])/par[1]^3) * G1(z) - (2/par[1]^3)*L1(z) + (1/(par[1]^2*par[2]))*L2(z) - ((1-par[2])/(par[1]^2*par[2])) * G2(z) )
    psi1.2 <- function(z) 2*( (1/(par[1]^2*par[2])) * (L1(z) - G1(z)) )
    psi2.1 <- function(z) 2*( ((2-par[2])/par[2]^3) * G2(z) - (2/par[2]^3)*L2(z) + (1/(par[2]^2*par[1]))*L1(z) - ((1-par[1])/(par[2]^2*par[1])) * G1(z) )
    psi2.2 <- function(z) 2*( (1/(par[2]^2*par[1])) * (L2(z) - G2(z)) )

    integrand.sigma2.11_g1diffg2 <- function(x,y) {
      if (support == "Continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2.1(x) * psi2.1(y) * estimVarCov_empProcess(x,y,samples[[2]]) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma2.12_g1diffg2 <- function(x,y) {
      if (support == "Continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2.1(x) * psi1.2(y) * estimVarCov_empProcess(x,y,samples[[2]]) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma2.13_g1diffg2 <- function(x) {
      if (support == "Continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2.1(x) * estimVarCov_empProcess(x,y,samples[[2]]) * densite.G.dataPoint.x
    }
    integrand.sigma2.22_g1diffg2 <- function(x,y) {
      if (support == "Continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1.2(x) * psi1.2(y) * estimVarCov_empProcess(x,y, samples[[2]]) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma2.23_g1diffg2 <- function(x) {
      if (support == "Continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1.2(x) * estimVarCov_empProcess(x,y, samples[[2]]) * densite.G.dataPoint.x
    }
    integrand.sigma2.31_g1diffg2 <- function(y) {
      if (support == "Continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi2.1(y) * estimVarCov_empProcess(x,y, samples[[2]]) * densite.G.dataPoint.y
    }
    integrand.sigma2.32_g1diffg2 <- function(y) {
      if (support == "Continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      psi1.2(y) * estimVarCov_empProcess(x,y, samples[[2]]) * densite.G.dataPoint.y
    }

    sigma <- matrix(NA, nrow = 3, ncol = 3)
    if (support == "Continuous") {
      sigma[1,1] <- pracma::integral2(Vectorize(integrand.sigma2.11_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,2] <- sigma[2,1] <- pracma::integral2(Vectorize(integrand.sigma2.12_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,3] <- stats::integrate(Vectorize(integrand.sigma2.13_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,2] <- pracma::integral2(Vectorize(integrand.sigma2.22_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[2,3] <- stats::integrate(Vectorize(integrand.sigma2.23_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,1] <- stats::integrate(Vectorize(integrand.sigma2.31_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,2] <- stats::integrate(Vectorize(integrand.sigma2.32_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,3] <- estimVarCov_empProcess(x, y, samples[[2]])
    } else {
      sigma[1,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma2.11_g1diffg2, x, supp.integration)) )
      sigma[1,2] <- sigma[2,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma2.12_g1diffg2, x, supp.integration)) )
      sigma[1,3] <- sum( unlist(sapply(supp.integration, integrand.sigma2.13_g1diffg2)) )
      sigma[2,2] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma2.22_g1diffg2, x, supp.integration)) )
      sigma[2,3] <- sum( unlist(sapply(supp.integration, integrand.sigma2.23_g1diffg2)) )
      sigma[3,1] <- sum( unlist(sapply(supp.integration, integrand.sigma2.31_g1diffg2)) )
      sigma[3,2] <- sum( unlist(sapply(supp.integration, integrand.sigma2.32_g1diffg2)) )
      sigma[3,3] <- estimVarCov_empProcess(x, y, samples[[2]])
    }
  }

  return(sigma)
}


#' Hessian matrix of the contrast function when using IBM
#'
#' Compute the hessian matrix of the contrast as defined in the IBM approach, at point (p1,p2). Here, based on
#' two samples following admixture models, where we recall that admixture models have probability distribution function
#' (pdf) given by l where l = p*f + (1-p)*g, where g represents the only known quantity and l is the pdf of the observed sample.
#' See 'Details' below for further information about the definition of the contrast.
#'
#' @param par Numeric vector with two elements (corresponding to the two unknown component weights) at which the hessian is computed.
#' @param IBMestim.obj An object of class 'estim_IBM'.
#' @param samples (List) List of the two considered samples.
#' @param admixMod (List) List of objects of class 'admix_model', one for each sample.
#' @param G Distribution on which to integrate when calculating the contrast.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return the hessian matrix of the contrast.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 20000, weight = 0.2,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 2, "scale" = 3),
#'                                         list("rate" = 1/3)))
#' data1 <- getmixtData(mixt1)
#' mixt2 <- twoComp_mixt(n = 10000, weight = 0.5,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 2, "scale" = 3),
#'                                         list("rate" = 1/5)))
#' data2 <- getmixtData(mixt2)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#'
#' ## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
#' est <- estim_IBM(samples = list(data1,data2),
#'                  admixMod = list(admixMod1,admixMod2), n.integ = 1000)
#'
#' ## Define the distribution over which to integrate:
#' fit.all <- stats::density(x = c(data1, data2))
#' G <- stats::rnorm(n = 1000, mean = sample(c(data1, data2),
#'                                           size = 1000, replace = TRUE), sd = fit.all$bw)
#' ## Evaluate the hessian matrix at point (p1,p2) = (0.3,0.6):
#' IBM_hessian_contrast(par = c(0.3,0.6), IBMestim.obj = est, samples = list(data1,data2),
#'                      admixMod = list(admixMod1, admixMod2), G = G)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

IBM_hessian_contrast <- function(par, IBMestim.obj, samples, admixMod, G)
{
  ## Extract the information on component distributions:
  knownCDF_comp.dist <- paste0("p", unlist(sapply(admixMod, '[[', 'comp.dist')["known", ]))
  if (any(knownCDF_comp.dist == "pmultinom")) knownCDF_comp.dist[which(knownCDF_comp.dist == "pmultinom")] <- "stepfun"
  comp_emp <- sapply(X = knownCDF_comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])

  ## Create the expression involved in future assessments of the CDF:
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

  ## Cumulative distribution function (either empirical from the two observed samples, or theoretical):
  L1 <- stats::ecdf(samples[[1]])
  L2 <- stats::ecdf(samples[[2]])

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])

  ## Loi de G et support d'integration:
  support <- detect_support_type(samples[[1]], samples[[2]])
  if (support == "Continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  ## Differentiates cases where G1 = G2 and G1 != G2 :
  if (G1equalG2) {

    stopifnot(length(par) == 1)
    integrand <- function(z, par) {
      if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      F.X.dataPoint <- (1/IBMestim.obj$p.X.fixed) * (L1(z) - (1-IBMestim.obj$p.X.fixed) * G1(z))
      F.Y.dataPoint <- (1/par) * (L2(z) - (1-par) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.dataPoint <- (1/par^2) * (L2(z) - G2(z))
      hessianD.dataPoint <- -(2/par^3) * (L2(z) - G2(z))
      (hessianD.dataPoint * D + gradD.dataPoint^2) * densite.G.dataPoint
    }

    if (support == "Continuous") {
      err.tol <- 1e-03
      hessian.contr <- try(suppressWarnings(
        2 * stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par,
                             subdivisions = 10000L, rel.tol = err.tol)$value), silent = TRUE)
      while ( (class(hessian.contr) == "character") & (err.tol < 1) ) {
        err.tol <- err.tol * 10
        hessian.contr <- try(suppressWarnings(
          2 * stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par,
                               subdivisions = 10000L, rel.tol = err.tol)$value), silent = TRUE)
      }
    } else {
      hessian.contr <- 2 * sum( unlist(sapply(supp.integration, integrand, par)) )
    }

  } else {

    integrand.1.1 <- function(z, par) {
      if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
      F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.pX.dataPoint <- (1/par[1]^2) * (G1(z) - L1(z))
      hessianD.dataPoint <- -(2/par[1]^3) * (G1(z) - L1(z))
      (gradD.pX.dataPoint^2 + D * hessianD.dataPoint) * densite.G.dataPoint
    }
    integrand.1.2 <- integrand.2.1 <- function(z, par) {
      if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      gradD.pX.dataPoint <- (1/par[1]^2) * (G1(z) - L1(z))
      gradD.pY.dataPoint <- -(1/par[2]^2) * (G2(z) - L2(z))
      (gradD.pX.dataPoint * gradD.pY.dataPoint) * densite.G.dataPoint
    }
    integrand.2.2 <- function(z, par) {
      if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
      }
      F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
      F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.pY.dataPoint <- -(1/par[2]^2) * (G2(z) - L2(z))
      hessianD.dataPoint <- (2/par[1]^3) * (G2(z) - L2(z))
      (gradD.pY.dataPoint^2 + D * hessianD.dataPoint) * densite.G.dataPoint
    }

    hessian.contr <- matrix(NA, nrow = 2, ncol = 2)
    if (support == "Continuous") {

      err.tol <- 1e-03
      term1 <- try(2 * stats::integrate(integrand.1.1, lower = supp.integration[1], upper = supp.integration[2], par,
                                        subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term1) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term1 <- try(2 * stats::integrate(integrand.1.1, lower = supp.integration[1], upper = supp.integration[2], par,
                                          subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      }
      hessian.contr[1,1] <- term1

      err.tol <- 1e-03
      term2 <- try(2 * stats::integrate(integrand.1.2, lower = supp.integration[1], upper = supp.integration[2], par,
                                        subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term2) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term2 <- try(2 * stats::integrate(integrand.1.2, lower = supp.integration[1], upper = supp.integration[2], par,
                                          subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      }
      hessian.contr[1,2] <- term2
      hessian.contr[2,1] <- hessian.contr[1,2]

      err.tol <- 1e-03
      term3 <- try(2 * stats::integrate(integrand.2.2, lower = supp.integration[1], upper = supp.integration[2], par,
                                        subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term3) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term3 <- try(2 * stats::integrate(integrand.2.2, lower = supp.integration[1], upper = supp.integration[2], par,
                                          subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
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
#' @param samples (List) List of the two considered samples.
#' @param admixMod (List) List of objects of class 'admix_model', one for each sample.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return the gap evaluated at the specified point between the unknown components of the two observed samples.
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
#'                       comp.param = list(list("mean" = 1, "sd" = 0.1),
#'                                         list("mean" = 5, "sd" = 2)))
#' data2 <- getmixtData(mixt2)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#'
#' IBM_gap(z = 2.8, par = c(0.3,0.6), samples = list(data1,data2),
#'         admixMod = list(admixMod1,admixMod2))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

IBM_gap <- function(z, par, samples, admixMod)
{
  ## Extract the information on component distributions:
  knownCDF_comp.dist <- paste0("p", unlist(sapply(admixMod, '[[', 'comp.dist')["known", ]))
  if (any(knownCDF_comp.dist == "pmultinom")) knownCDF_comp.dist[which(knownCDF_comp.dist == "pmultinom")] <- "stepfun"
  comp_emp <- sapply(X = knownCDF_comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])

  ## Create the expression involved in future assessments of the CDF:
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

  ## Cumulative distribution function (either empirical from the two observed samples, or theoretical):
  L1 <- stats::ecdf(samples[[1]])
  L2 <- stats::ecdf(samples[[2]])

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])
  warning("In function 'IBM_gap': in case of identical known components, 'fixed.p1' was arbitrarily set to 0.5.")
  if (G1equalG2) { fixed.p1 <- 0.5
  } else { fixed.p1 <- NULL }

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


IBM_theoretical_gap <- function(z, par, known.p, mixtMod)
{
  if (is.null(known.p)) stop("In 'IBM_theoretical_gap': argument 'known.p' must be specified.")
  stopifnot("In 'IBM_theoretical_gap', wrong class for argument 'mixtMod'" =
              all(sapply(X = mixtMod, FUN = "class") == c("twoComp_mixt","twoComp_mixt")))

  ## Extract the information on component distributions:
  CDF_comp.dist <- paste0("p", unlist(sapply(mixtMod, '[[', 'comp.dist')))
  if (any(CDF_comp.dist == "pmultinom")) CDF_comp.dist[which(CDF_comp.dist == "pmultinom")] <- "stepfun"
  comp_emp <- sapply(X = CDF_comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])

  ## Create the expression involved in future assessments of the CDF:
  list.param_dist <- list(mixtMod[[1]]$comp.param[[1]], mixtMod[[1]]$comp.param[[2]],
                          mixtMod[[2]]$comp.param[[1]], mixtMod[[2]]$comp.param[[2]])
  make.expr.step <- function(i) paste(names(comp_emp)[i],"(x = 1:", length(list.param_dist[[i]]$prob),
                                      paste(", y = ", paste("cumsum(c(0,", paste(list.param_dist[[i]]$prob, collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_emp)[i],"(z,", paste(names(list.param_dist[[i]]),
                                                                 "=", list.param_dist[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(CDF_comp.dist))
  expr[which(CDF_comp.dist == "stepfun")] <- sapply(which(CDF_comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  if (any(CDF_comp.dist == "stepfun")) {
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
