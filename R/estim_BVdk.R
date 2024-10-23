#' Estimation of the admixture parameters by Bordes & Vandekerkhove (2010)
#'
#' Estimates parameters in an admixture model where the unknown component is assumed to have a symmetric density.
#' More precisely, estimates the two parameters (mixture weight and location shift) in the admixture model with pdf:
#'          l(x) = p*f(x-mu) + (1-p)*g(x), x in R,
#' where g is the known component, p is the proportion and f is the unknown component with symmetric density.
#' The localization shift parameter is denoted mu, and the component weight p.
#' See 'Details' below for further information.
#'
#' @param data The observed sample under study.
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
#'
#' ## Perform the estimation of parameters in real-life:
#' estim_BVdk(data = data1, admixMod = admixMod, method = 'L-BFGS-B')
#'
#' ## Second example:
#' mixt2 <- twoComp_mixt(n = 200, weight = 0.65,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean" = -1, "sd" = 0.5),
#'                                         list("rate" = 1)))
#' data2 <- getmixtData(mixt2)
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                         knownComp_param = mixt2$comp.param[[2]])
#'
#' ## Perform the estimation of parameters in real-life:
#' estim_BVdk(data = data2, admixMod = admixMod2, method = 'L-BFGS-B')
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

estim_BVdk <- function(data, admixMod, method = c("L-BFGS-B","Nelder-Mead"))
{
  ## Extract useful information about known component distribution:
  comp.dist.sim <- paste0("r", admixMod$comp.dist$known)
  comp.sim <- sapply(X = comp.dist.sim, FUN = get, mode = "function")
  assign(x = names(comp.sim)[1], value = comp.sim[[1]])
  expr.sim <- paste(names(comp.sim)[1],"(n=100000,", paste(names(admixMod$comp.param$known),
                    "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  ## Initialization of the parameters: localization parameter is initialized depending on whether the global mean
  ## of the sample is lower than the mean of the known component (or not).
  if (mean(data) > mean(eval(parse(text = expr.sim)))) {
    init.param <- c(0.5, 0.99 * max(data))
  } else {
    init.param <- c(0.5, (min(data) + 0.01 * abs(min(data))))
  }

  ## Select the bandwith :
  bandw <- stats::density(data)$bw

  if (method == "Nelder-Mead") {
    sol <- stats::optim(par = init.param, fn = BVdk_contrast, gr = BVdk_contrast_gradient, data = data,
                        admixMod = admixMod, h = bandw, method = "Nelder-Mead")
  } else {
    sol <- stats::optim(par = init.param, fn = BVdk_contrast, gr = BVdk_contrast_gradient, data = data, admixMod = admixMod,
                        h = bandw, method = "L-BFGS-B", lower = c(0.001,min(data)), upper = c(0.999,max(data)))
  }

  obj <- list(
    n_populations = 1,
    population_sizes = length(data),
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
#'                       comp.param = list(c("mean" = 3, "sd" = 0.5),
#'                                         c("mean" = 0, "sd" = 1)))
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
#'                       comp.param = list(c("mean" = 3, "sd" = 0.5),
#'                                         c("mean" = 0, "sd" = 1)))
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
