#' Estimation of the admixture parameters by Bordes & Vandekerkhove (2010)
#'
#' Estimates parameters in an admixture model where the unknown component is assumed to have a symmetric density.
#' More precisely, estimates the two parameters (mixture weight and location shift) in the admixture model with pdf:
#'          l(x) = p*f(x-mu) + (1-p)*g(x), x in R,
#' where g is the known component, p is the proportion and f is the unknown component with symmetric density.
#' The localization shift parameter is denoted mu, and the component weight p.
#' See the reference below for further details.
#'
#' @param samples The observed sample under study.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#' @param method The method used throughout the optimization process, either 'L-BFGS-B' or 'Nelder-Mead' (see ?optim).
#'
#' @references
#' \insertRef{BordesVandekerkhove2010}{admix}
#'
#' @return An object of class 'estim_BVdk', containing 7 attributes: 1) the number of sample under study (set to 1 here);
#'         2) the sample size; 3) the information about mixture components (distributions and parameters); 4) the estimation
#'         method (Bordes and Vandekerkhove here, see the given reference); 5) the estimated mixing proportion (weight of the
#'         unknown component distribution); 6) the estimated location parameter of the unknown component distribution (with symetric
#'         density); 7) the optimization method that was used.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 200, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' ## Retrieves the mixture data:
#' data1 <- getmixtData(mixt1)
#' ## Define the admixture model:
#' admixMod <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                         knownComp_param = mixt1$comp.param[[2]])
#' ## Perform the estimation of parameters in real-life:
#' estim_BVdk(samples = data1, admixMod = admixMod, method = 'L-BFGS-B')
#'
#' ## Second example:
#' mixt2 <- twoComp_mixt(n = 200, weight = 0.65,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean" = -1, "sd" = 0.5),
#'                                         list("rate" = 1)))
#' data2 <- getmixtData(mixt2)
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                         knownComp_param = mixt2$comp.param[[2]])
#' ## Perform the estimation of parameters in real-life:
#' estim_BVdk(samples = data2, admixMod = admixMod2, method = 'L-BFGS-B')
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

estim_BVdk <- function(samples, admixMod, method = c("L-BFGS-B","Nelder-Mead"))
{
  warning("Estimation by 'BVdk' assumes the unknown component distribution
  to have a symmetric probability density function.\n")
  ## Extract useful information about known component distribution:
  comp.dist.sim <- paste0("r", admixMod$comp.dist$known)
  comp.sim <- sapply(X = comp.dist.sim, FUN = get, mode = "function")
  assign(x = names(comp.sim)[1], value = comp.sim[[1]])
  expr.sim <- paste(names(comp.sim)[1],"(n=100000,", paste(names(admixMod$comp.param$known),
                    "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  ## Initialization of the parameters: localization parameter is initialized depending on whether the global mean
  ## of the sample is lower than the mean of the known component (or not).
  if (mean(samples) > mean(eval(parse(text = expr.sim)))) {
    init.param <- c(0.5, 0.99 * max(samples))
  } else {
    init.param <- c(0.5, (min(samples) + 0.01 * abs(min(samples))))
  }

  ## Select the bandwith :
  bandw <- stats::density(samples)$bw

  method <- match.arg(method)
  if (method == "Nelder-Mead") {
    sol <- stats::optim(par = init.param, fn = BVdk_contrast, gr = BVdk_contrast_gradient, data = samples,
                        admixMod = admixMod, h = bandw, method = "Nelder-Mead")
  } else {
    sol <- stats::optim(par = init.param, fn = BVdk_contrast, gr = BVdk_contrast_gradient, data = samples, admixMod = admixMod,
                        h = bandw, method = "L-BFGS-B", lower = c(0.001,min(samples)), upper = c(0.999,max(samples)))
  }

  obj <- list(
    n_populations = 1,
    population_sizes = length(samples),
    admixture_models = admixMod,
    estimation_method = "Bordes and Vandekerkhove",
    estimated_mixing_weights = sol$par[1],
    estimated_locations = sol$par[2],
    optim_method = method
  )
  class(obj) <- c("estim_BVdk", "admix_estim")
  obj$call <- match.call()

  return(obj)
}


#' Print method for objects 'estim_BVdk'
#'
#' Print the results stored in an object of class 'estim_BVdk'.
#'
#' @param x An object of class 'estim_BVdk'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.estim_BVdk <- function(x, ...)
{
  cat("\nCall:")
  print(x$call)
  cat("\n")
  cat("Estimated mixing proportion: ", x$estimated_mixing_weights, "\n")
  cat("Estimated location parameter: ", x$estimated_locations, "\n\n")
}


#' Summary method for objects 'estim_BVdk'
#'
#' Summarizes the results stored in an object of class 'estim_BVdk'.
#'
#' @param object An object of class 'estim_BVdk'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.estim_BVdk <- function(object, ...)
{
  cat("Call:")
  print(object$call)
  cat("\n")
  cat("------- Sample -------\n")
  cat("Sample size: ", object$population_sizes, "\n")
  cat("-> Distribution and parameters of the known component \n in the admixture model: ", sep="")
  cat(object$admixture_models$comp.dist$known, "\n")
  print(unlist(object$admixture_models$comp.param$known, use.names = TRUE))
  cat("\n------- Estimation results -------\n")
  cat("Estimated mixing proportion: ", object$estimated_mixing_weights, "\n")
  cat("Estimated location parameter: ", object$estimated_locations, "\n\n")
  cat("------- Optimization -------\n")
  cat("Optimization method: ", object$optim_method, "\n\n")
}


#' Contrast as defined in Bordes & Vandekerkhove (2010)
#'
#' Compute the contrast as defined in Bordes & Vandekerkhove (2010) (see below in section 'Details'), needed for
#' optimization purpose. Remind that one considers an admixture model with symmetric unknown density, i.e.
#'          l(x) = p*f(x-mu) + (1-p)*g(x),
#' where l denotes the probability density function (pdf) of the mixture with known component pdf g, p is the unknown mixture
#' weight, f relates to the unknown symmetric component pdf f, and mu is the location shift parameter.

#' @param param Numeric vector of two elements, corresponding to the two parameters (first the unknown component weight, and
#'              then the location shift parameter of the symmetric unknown component distribution).
#' @param data Numeric vector of observations following the mixture model given by the pdf l.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#' @param h Width of the window used in the kernel estimations.
#'
#' @references
#' \insertRef{BordesVandekerkhove2010}{admix}
#'
#' @return The value of the contrast.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1000, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' ## Define the admixture model:
#' admixMod <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                         knownComp_param = mixt1$comp.param[[2]])
#' ## Compute the contrast value for some given parameter vector in real-life framework:
#' BVdk_contrast(param = c(0.3,2), data = data1,
#'               admixMod = admixMod, h = density(data1)$bw)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

BVdk_contrast <- function(param, data, admixMod, h)
{
  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("p", admixMod$comp.dist$known)
  comp_BVdk <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  assign(x = names(comp_BVdk)[1], value = comp_BVdk[[1]])

  ## Creates the expression allowing further to generate the right data:
  expr1 <- paste(names(comp_BVdk)[1],"(data[i] + mu,", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  expr2 <- paste(names(comp_BVdk)[1],"(mu - data[i],", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")

  p <- param[1]
  mu <- param[2]

  ## Here, 'G' is the cdf of the admixture model, Fo is the cdf of the known component, and H is the cdf of the unknown component.
  n <- length(data)
  H <- G <- Fo <- rep(0, n)
  for (i in 1:n) {
    G[i]  <- mean(kernel_cdf(data[i] + mu - data, h)) + mean(kernel_cdf(mu - data[i] - data, h))
    Fo[i] <- eval(parse(text = expr1)) + eval(parse(text = expr2))
  }

  ## cf formula (2.3) p.5 :
  H <- ( (G - (1-p) * Fo) / p ) - 1

  return( mean(H^2) )
}


#' Gradient of the contrast defined in Bordes & Vandekerkhove (2010)
#'
#' Compute the gradient of the contrast as defined in Bordes & Vandekerkhove (2010) (see below in section 'Details'), needed for optimization purpose. Remind
#' that one considers an admixture model, i.e. l = p*f + (1-p)*g ; where l denotes the probability density function (pdf) of
#' the mixture with known component pdf g, p is the unknown mixture weight, and f relates to the  unknown symmetric
#' component pdf f.
#'
#' @param param A numeric vector with two elements corresponding to the parameters to be estimated. First the unknown component
#'              weight, and second the location shift parameter of the symmetric unknown component distribution.
#' @param data A vector of observations following the admixture model given by the pdf l.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#' @param h The window width used in the kernel estimations.
#'
#' @references
#' \insertRef{BordesVandekerkhove2010}{admix}
#'
#' @return A numeric vector composed of the two partial derivatives w.r.t. the two parameters on which to optimize the contrast.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1000, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' ## Define the admixture model:
#' admixMod <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                         knownComp_param = mixt1$comp.param[[2]])
#' BVdk_contrast_gradient(param = c(0.3,2), data = data1,
#'                        admixMod = admixMod, h = density(data1)$bw)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

BVdk_contrast_gradient <- function(param, data, admixMod, h)
{
  ## Extracts the information on component distributions and stores in expressions:
  comp.dist.cdf <- paste0("p", admixMod$comp.dist$known)
  comp_cdf <- sapply(X = comp.dist.cdf, FUN = get, mode = "function")
  assign(x = names(comp_cdf)[1], value = comp_cdf[[1]])
  ## Creates the expression allowing further to generate the right data:
  expr1.cdf <- paste(names(comp_cdf)[1],"(data[i] + mu,", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  expr2.cdf <- paste(names(comp_cdf)[1],"(mu - data[i],", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")

  ## Same with density functions :
  comp.dist.dens <- paste0("d", admixMod$comp.dist$known)
  comp_dens <- sapply(X = comp.dist.dens, FUN = get, mode = "function")
  assign(x = names(comp_dens)[1], value = comp_dens[[1]])
  ## Creates the expression allowing further to generate the right data:
  expr1.dens <- paste(names(comp_dens)[1],"(data[i] + mu,", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  expr2.dens <- paste(names(comp_dens)[1],"(mu - data[i],", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")

  p <- param[1]
  mu <- param[2]

  ## Here, 'G' is the cdf of the admixture model, Fo is the cdf of the known component, and H is the cdf of the unknown component.
  ## Lowercase letters refer to the corresponding densities.
  n <- length(data)
  H <- G <- Fo <- g <- fo <- rep(0, n)

  for (i in 1:n) {
    G[i]  <- mean( kernel_cdf(data[i] + mu - data, h) + kernel_cdf(mu - data[i] - data, h) )
    g[i]  <- mean( kernel_density(data[i] + mu - data, h) + kernel_density(mu - data[i] - data, h) )
    Fo[i] <- eval(parse(text = expr1.cdf)) + eval(parse(text = expr2.cdf))
    fo[i] <- eval(parse(text = expr1.dens)) + eval(parse(text = expr2.dens))
  }

  ## Inversion formula to isolate the unknown component cdf:
  H <- ( (G - (1-p) * Fo) / p ) - 1
  ## Partial derivative with respect to the component weight 'p':
  d_p_H <- (Fo - G) / p^2
  ## Partial derivative with respect to the location parmaeter 'mu':
  d_mu_H <- (g - (1-p) * fo) / p

  return( c(mean(H * d_p_H), mean(H * d_mu_H)) )
}


#' Variance of estimators in an admixture model with symmetric unknown density.
#'
#' Semiparametric estimation of the variance of the estimators related to the mixture weight p and the location shift parameter mu,
#' considering the admixture model with probability density function l:
#'          l(x) = p*f(x-mu) + (1-p)*g(x), x in R,
#' where g is the known component of the two-component mixture, p is the unknown proportion, f is the unknown component density and
#' mu is the location shift. See 'Details' below for more information.
#'
#' @param estim An object of class 'estim_BVdk', containing the estimators of unknown quantities in the admixture model.
#' @param data The observed sample under study.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#'
#' @details See formulas pp.28--30 in Appendix of Bordes, L. and Vandekerkhove, P. (2010).
#'
#' @references
#' \insertRef{BordesVandekerkhove2010}{admix}
#'
#' @return A list containing 1) the variance-covariance matrix of the estimators (assessed at the specific time points 'u' and 'v'
#'         such that u = v = mean(data)); 2) the variance of the estimator of the unknown mixture weight; 3) the variance of the
#'         estimator of the location shift parameter; 4) the variance of the estimator of the unknown component cumulative
#'         distribution function at points 'u' and 'v'.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 200, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' ## Define the admixture model:
#' admixMod <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                         knownComp_param = mixt1$comp.param[[2]])
#'
#' ## Perform the estimation of parameters in real-life:
#' estim <- estim_BVdk(samples = data1, admixMod = admixMod, method = 'L-BFGS-B')
#' BVdk_varCov_estimators(estim = estim, data = data1, admixMod = admixMod)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

BVdk_varCov_estimators <- function(estim, data, admixMod)
{
  n <- length(data)
  bandw <- stats::density(data)$bw
  u <- v <- round(mean(data))

  inv_J <- solve(BVdk_mat_J(estim, data, admixMod, h = bandw))
  mat_L.u <- BVdk_mat_L(u = u, estim, data, admixMod, h = bandw)
  mat_L.v <- BVdk_mat_L(u = v, estim, data, admixMod, h = bandw)
  sigma.mat.uv <- BVdk_mat_Sigma(u = u, v = v, estim, data, admixMod, h = bandw)
  varCov.matrix <- (1/n) * (mat_L.u %*% inv_J %*% sigma.mat.uv %*% inv_J %*% t(mat_L.v))

  return( list(varCov.mat = varCov.matrix,
               var.estim_prop = varCov.matrix[1,1],
               var.estim_location = varCov.matrix[2,2],
               var.estim_cdf = varCov.matrix[3,3]) )
}


## Estimation of matrix Sigma(u,v) from empirical versions (cf p.28 et p.29), with arguments:
## - u the time point at which the first (related to the first parameter) underlying empirical process is looked through.
## - v the time point at which the second (related to the second parameter) underlying empirical process is looked through.
BVdk_mat_Sigma <- function(u, v, estim, data, admixMod, h)
{
  n <- length(data)
  aux1 <- aux2 <- rep(0, n)
  aux1 <- sapply(X = data, FUN = h1, estim, admixMod)
  aux2 <- sapply(X = data, FUN = h2, data, estim, admixMod, h)

  sorted_data <- sort(data)
  ul <- sapply(X = data, FUN = l_fun, u, sorted_data, estim)
  vl <- sapply(X = data, FUN = l_fun, v, sorted_data, estim)

  aux12 <- terms_sigma11 <- terms_sigma22 <- terms_sigma12 <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      aux12[i,j] <- ka(data[i], data[j], sorted_data, estim)
      terms_sigma11[i,j] <- aux1[i] * aux1[j] * aux12[i,j]
      terms_sigma22[i,j] <- aux2[i] * aux2[j] * aux12[i,j]
      terms_sigma12[i,j] <- (aux1[i] * aux2[j] + aux1[j] * aux2[i]) * aux12[i,j]
    }
  }

  sigma_mat <- matrix(NA, nrow = 3, ncol = 3)
  sigma_mat[3,3] <- .Call('_admix_Donsker_correl_cpp',
                          PACKAGE = 'admix',
                          u + estim$estimated_locations,
                          v + estim$estimated_locations, sorted_data)
  sigma_mat[1,3] <- (2 / (n * estim$estimated_mixing_weights)) * sum(aux1 * vl)
  sigma_mat[3,1] <- (2 / (n * estim$estimated_mixing_weights)) * sum(aux1 * ul)
  sigma_mat[2,3] <- (2 / (n * estim$estimated_mixing_weights)) * sum(aux2 * vl)
  sigma_mat[3,2] <- (2 / (n * estim$estimated_mixing_weights)) * sum(aux2 * ul)
  sigma_mat[1,1] <- (8 / (n * (n-1) * estim$estimated_mixing_weights^2)) * sum(terms_sigma11, na.rm = TRUE)
  sigma_mat[2,2] <- (8 / (n * (n-1) * estim$estimated_mixing_weights^2)) * sum(terms_sigma22, na.rm = TRUE)
  sigma_mat[1,2] <- sigma_mat[2,1] <- (4 / (n * (n-1) * estim$estimated_mixing_weights^2)) * sum(terms_sigma12, na.rm = TRUE)

  return(sigma_mat)
}


## Define the matrix J that is useful for the consistency theorem of estimators (cf p.29-30) .
BVdk_mat_J <- function(estim, data, admixMod, h)
{
  n <- length(data)
  aux1 <- aux2 <- rep(0, n)
  aux1 <- sapply(X = data, FUN = h1, estim, admixMod)
  aux2 <- sapply(X = data, FUN = h2, data, estim, admixMod, h)
  J <- matrix(0, nrow = 3, ncol = 3)
  J[1,1] <- -(2/n) * sum(aux1^2)
  J[2,2] <- -(2/n) * sum(aux2^2)
  J[1,2] <- -(2/n) * sum(aux1 * aux2)
  J[2,1] <- J[1,2]
  J[3,3] <- 1
  return(J)
}

## Define the matrix L, useful for the consistency theorem of estimators (cf p.29-30).
## Returns the evaluation of the matrix terms at point 'u'.
BVdk_mat_L <- function(u, estim, data, admixMod, h)
{
  L <- matrix(0, nrow = 3, ncol = 3)
  L[1,1] <- L[2,2] <- 1
  L[3,1] <- h3(u, data, estim, admixMod)
  L[3,2] <- h2(u, data, estim, admixMod, h) / 2
  L[3,3] <- 1 / estim$estimated_mixing_weights
  return(L)
}


##### All the following functions are defined p.28 of the paper (see p.29 for empirical versions).

## Function h1, see formula (3.9) p.13
h1 <- function(u, estim, admixMod)
{
  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("p", admixMod$comp.dist$known)
  comp_h1 <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  assign(x = names(comp_h1)[1], value = comp_h1[[1]])
  expr1 <- paste(names(comp_h1)[1],"(q=(estim$estimated_locations+u),",
                 paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  expr2 <- paste(names(comp_h1)[1],"(q=(estim$estimated_locations-u),",
                 paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  res <- (eval(parse(text = expr1)) + eval(parse(text = expr2)) - 1) / estim$estimated_mixing_weights
  return(res)
}

## Fonction h2, see formula (3.8) p.13 (see also the connexion with formula end p.6)
h2 <- function(u, data, estim, admixMod, h)
{
  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("d", admixMod$comp.dist$known)
  comp_h2 <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  assign(x = names(comp_h2)[1], value = comp_h2[[1]])
  expr <- paste(names(comp_h2)[1],"(x=(estim$estimated_locations+u),",
                paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")

  g  <- mean( kernel_density(u + estim$estimated_locations - data, h) )
  fo <- eval(parse(text = expr))
  f  <- (g - (1-estim$estimated_mixing_weights) * fo) / estim$estimated_mixing_weights			# cf (2.2) p.5

  res <- 2 * f * (f >= 0)         # end p.28
  return(res)
}

## Fonction h3:
h3 <- function(u, data, estim, admixMod)
{
  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("p", admixMod$comp.dist$known)
  comp_h3 <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  assign(x = names(comp_h3)[1], value = comp_h3[[1]])
  expr <- paste(names(comp_h3)[1],"(q=(estim$estimated_locations+u),",
                paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  G_ecdf <- stats::ecdf(data)

  res <- (eval(parse(text = expr)) - G_ecdf(u + estim$estimated_locations)) / estim$estimated_mixing_weights^2
  return(res)
}

## Function ka of the paper:
ka <- function(u, v, data, estim)
{
  return( l_fun(u, v, data, estim) + l_fun(-u, v, data, estim) )
}

## Old implementation of fonction k in the article:
ka_old <- function(u, v, data, estim)
{
  return( l_fun_old(u, v, data, estim) + l_fun_old(-u, v, data, estim) )
}

## Function l of the article:
l_fun <- function(u, v, data, estim)
{
  term1 <- .Call('_admix_Donsker_correl_cpp', PACKAGE = 'admix', estim$estimated_locations + u, estim$estimated_locations + v, data)
  term2 <- .Call('_admix_Donsker_correl_cpp', PACKAGE = 'admix', estim$estimated_locations + u, estim$estimated_locations - v, data)
  return(term1 + term2)
}

## Old implementation of fonction l in the article:
l_fun_old <- function(u, v, data, estim)
{
  return( Donsker_correl_old(estim$estimated_locations + u, estim$estimated_locations + v, data) +
            Donsker_correl_old(estim$estimated_locations + u, estim$estimated_locations - v, data) )
}

## Donsker correlation:
Donsker_correl_old <- function(u, v, obs.data)
{
  L.CDF <- stats::ecdf(obs.data)
  return( L.CDF(min(u,v)) * (1 - L.CDF(max(u,v))) )
}


#' Variance of the variance estimator for the unknown component in an admixture model
#'
#' Maximum Likelihood estimation of the variance of the variance parameter in Bordes & Vandekerkhove (2010) setting,
#' i.e. considering the admixture model with probability density function (pdf) l:
#'          l(x) = p*f(x-mu) + (1-p)*g,
#' where g is the known component of the two-component mixture, p is the mixture proportion, f is the unknown component with
#' symmetric density, and mu is the location shift parameter. The estimation of the variance of the variance related to the density
#' f is made by maximum likelihood optimization through the information matrix, with the assumption that the unknown f is gaussian.
#'
#' @param data The observed sample under study.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#' @param hat_w Estimate of the unknown component weight.
#' @param hat_loc Estimate of the location shift parameter.
#' @param hat_var Estimate of the variance of the symmetric density f, obtained by plugging-in the previous estimates. See 'Details'
#'                below for further information.
#'
#' @details Plug-in strategy is defined in Pommeret, D. and Vandekerkhove, P. (2019). The variance of the estimator variance of the unknown density f
#'          is needed in a testing perspective, since included in the variance of the test statistic. Other details about the information
#'          matrix can be found in Bordes, L. and Vandekerkhove, P. (2010).
#'
#' @references
#' \insertRef{BordesVandekerkhove2010}{admix}
#' \insertRef{PommeretVandekerkhove2019}{admix}
#'
#' @return The variance of the estimator of the variance of the unknown component density f.
#'
#' @examples
#' \donttest{
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 400, weight = 0.9,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 4, "sd" = 1),
#'                                         list("mean" = 7, "sd" = 0.5)))
#' data1 <- getmixtData(mixt1)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#'
#' estim <- estim_BVdk(data = data1, admixMod = admixMod1, method = "L-BFGS-B")
#' ## Estimation of the second-order moment of the known component distribution:
#' m2_knownComp <- mean(rnorm(n = 1000000, mean = 7, sd = 0.5)^2)
#' hat_s2 <- (1/estim$estimated_mixing_weights) * (mean(data1^2) -
#'           ((1-estim$estimated_mixing_weights)*m2_knownComp)) - estim$estimated_locations^2
#' ## Estimated variance of variance estimator related to the unknown symmetric component density:
#' BVdk_ML_varCov_estimators(data=data1, admixMod=admixMod1, hat_w=estim$estimated_mixing_weights,
#'                           hat_loc = estim$estimated_locations, hat_var = hat_s2)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

BVdk_ML_varCov_estimators <- function(data, admixMod, hat_w, hat_loc, hat_var)
{
  ## Same with density functions :
  comp.dist.dens <- paste0("d", admixMod$comp.dist$known)
  comp_dens <- sapply(X = comp.dist.dens, FUN = get, mode = "function")
  assign(x = names(comp_dens)[1], value = comp_dens[[1]])
  ## Creates the right expression depending on the component distributions:
  expr1 <- paste(names(comp_dens)[1],"(y, hat_loc, sqrt(hat_var))", sep = "")
  expr2 <- paste(names(comp_dens)[1],"(y,", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")

  f11 <- function(y) {
    f11.val <- (eval(parse(text=expr2)) - eval(parse(text=expr1)))^2 / (hat_w * eval(parse(text=expr2)) + (1-hat_w) * eval(parse(text=expr1)))
    return(f11.val)
  }
  f12 <- function(y) {
    f12.val <- (eval(parse(text=expr2)) - eval(parse(text=expr1))) * ( (1-hat_w) * (y - hat_loc) * eval(parse(text=expr1)) / hat_var ) /
      ( hat_w * eval(parse(text=expr2)) + (1-hat_w) * eval(parse(text=expr1)) )
    return(f12.val)
  }
  f22 <- function(y) {
    f22.val <- ( (1-hat_w) * (y - hat_loc) * eval(parse(text=expr1)) / hat_var )^2 /
      ( hat_w * eval(parse(text=expr2)) + (1-hat_w) * eval(parse(text=expr1)) )
    return(f22.val)
  }
  f13 <- function(y) {
    f13.val <- (eval(parse(text=expr2)) - eval(parse(text=expr1))) * ((1-hat_w) * eval(parse(text=expr1)) * ((y - hat_loc)^2 / (2*hat_var^2) - 1/(2*hat_var))) /
      (hat_w * eval(parse(text=expr2)) + (1-hat_w) * eval(parse(text=expr1)))
    return(f13.val)
  }
  f23 <- function(y) {
    f23.val <- ( (1-hat_w) * (y - hat_loc) * eval(parse(text=expr1)) / hat_var ) * ( (1-hat_w) * eval(parse(text=expr1)) * ((y - hat_loc)^2 / (2*hat_var^2) - 1/(2*hat_var)) ) /
      (hat_w * eval(parse(text=expr2)) + (1 - hat_w) * eval(parse(text=expr1)))
    return(f23.val)
  }
  f33 <- function(y) {
    f33.val <- ((1-hat_w) * eval(parse(text=expr1)) * ((y - hat_loc)^2 / (2*hat_var^2) - 1/(2*hat_var)))^2 /
      (hat_w * eval(parse(text=expr2)) + (1-hat_w) * eval(parse(text=expr1)))
    return(f33.val)
  }

  matvar <- matrix(rep(0,9), ncol = 3)
  ## Integrate such function the get the expectation of the Hessian matrix:
  ## (between lower bound=-20 and upper bound=20 => largely enough given that we look at gaussian distribution around 0)
  matvar[1,1] <- stats::integrate(f11, lower = (min(data)-0.1*abs(min(data))), upper = (max(data)+0.1*abs(min(data))))$value
  matvar[1,2] <- matvar[2,1] <- stats::integrate(f12, lower = (min(data)-0.1*abs(min(data))), upper = (max(data)+0.1*abs(min(data))))$value
  matvar[3,1] <- matvar[1,3] <- stats::integrate(f13, lower = (min(data)-0.1*abs(min(data))), upper = (max(data)+0.1*abs(min(data))))$value
  matvar[2,2] <- stats::integrate(f22, lower = (min(data)-0.1*abs(min(data))), upper = (max(data)+0.1*abs(min(data))))$value
  matvar[3,2] <- matvar[2,3] <- stats::integrate(f23, lower = (min(data)-0.1*abs(min(data))), upper = (max(data)+0.1*abs(min(data))))$value
  matvar[3,3] <- stats::integrate(f33, lower = (min(data)-0.1*abs(min(data))), upper = max(data))$value

  ## Take care: has to divide by 'n' (sample size) to get the estimated hessian
  variances <- (1/length(data)) * ( (-matvar[1,1] * matvar[1,2] + matvar[1,1] * matvar[2,2]) /
                                      (-(matvar[1,3])^2 * matvar[2,2] + matvar[1,1] * matvar[1,3] * matvar[2,3] + matvar[1,2] * matvar[1,3] * matvar[2,3] - matvar[1,1] * matvar[2,3]^2 -
                                         matvar[1,1] * matvar[1,2] * matvar[3,3] + matvar[1,1] * matvar[2,2] * matvar[3,3]) )
  if (variances < 0) warning("Plug-in strategy has led to a negative estimated variance with ML estimation.")

  return(variances)
}
