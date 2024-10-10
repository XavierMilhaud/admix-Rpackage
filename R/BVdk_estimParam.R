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
#' @param method The method used throughout the optimization process, either 'L-BFGS-B' or 'Nelder-Mead' (see ?optim).
#' @param comp.dist A list with two elements corresponding to component distributions (specified with R native names for these distributions) involved
#'                  in the admixture model. Unknown elements must be specified as 'NULL' objects, e.g. when 'f' is unknown: list(f=NULL, g='norm').
#' @param comp.param A list with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Unknown elements must be specified as 'NULL' objects, e.g. if 'f' is unknown: list(f=NULL, g=list(mean=0,sd=1)).
#'
#' @references
#' \insertRef{BordesVandekerkhove2010}{admix}
#'
#' @return A numeric vector with the two estimated parameters (proportion first, and then location shift).
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f = 'norm', g = 'norm')
#' list.param <- list(f = list(mean = -2, sd = 0.5),
#'                    g = list(mean = 0, sd = 1))
#' data1 <- rsimmix(n = 200, unknownComp_weight = 0.4, list.comp, list.param)[['mixt.data']]
#' ## Perform the estimation of parameters in real-life:
#' list.comp <- list(f = NULL, g = 'norm')
#' list.param <- list(f = NULL, g = list(mean = 0, sd = 1))
#' BVdk_estimParam(data1, method = 'L-BFGS-B', list.comp, list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

BVdk_estimParam <- function(data, method = c("L-BFGS-B","Nelder-Mead"), comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")

  ## Select the bandwith :
  bandw <- stats::density(data)$bw

  ## Initialization of the parameters: localization parameter is initialized depending on whether the global mean
  ## of the sample is lower than the mean of the known component (or not).
  comp.dist.sim <- paste0("r", comp.dist[[2]])
#  comp.sim <- sapply(X = comp.dist.sim, FUN = get, pos = "package:stats", mode = "function")
  comp.sim <- sapply(X = comp.dist.sim, FUN = get, mode = "function")
  assign(x = names(comp.sim)[1], value = comp.sim[[1]])
  expr.sim <- paste(names(comp.sim)[1],"(n=100000,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  if (mean(data) > mean(eval(parse(text = expr.sim)))) {
    init.param <- c(0.5, 0.99 * max(data))
  } else {
    init.param <- c(0.5, (min(data) + 0.01 * abs(min(data))))
  }

  if (method == "Nelder-Mead") {
    sol <- stats::optim(par = init.param, fn = BVdk_contrast, gr = BVdk_contrast_gradient, data = data, h = bandw, comp.dist = comp.dist,
                        comp.param = comp.param, method = "Nelder-Mead")

  } else {
    sol <- stats::optim(par = init.param, fn = BVdk_contrast, gr = BVdk_contrast_gradient, data = data, h = bandw, comp.dist = comp.dist,
                        comp.param = comp.param, method = "L-BFGS-B", lower = c(0.001,min(data)), upper = c(0.999,max(data)))
  }

  return(sol$par)
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
#' @param data Numeric vector of observations following the admixture model given by the pdf l.
#' @param h Width of the window used in the kernel estimations.
#' @param comp.dist A list with two elements corresponding to component distributions (specified with R native names for these distributions) involved
#'                  in the admixture model. Unknown elements must be specified as 'NULL' objects, e.g. when 'f' is unknown: list(f=NULL, g='norm').
#' @param comp.param A list with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Unknown elements must be specified as 'NULL' objects, e.g. if 'f' is unknown: list(f=NULL, g=list(mean=0,sd=1)).
#'
#' @references
#' \insertRef{BordesVandekerkhove2010}{admix}
#'
#' @return The value of the contrast.
#'
#' @examples
#' ## Simulate data:
#' comp.dist <- list(f = 'norm', g = 'norm')
#' comp.param <- list(f = list(mean = 3, sd = 0.5), g = list(mean = 0, sd = 1))
#' data1 <- rsimmix(n = 1000, unknownComp_weight = 0.6, comp.dist, comp.param)[['mixt.data']]
#' ## Compute the contrast value for some given parameter vector in real-life framework:
#' comp.dist <- list(f = NULL, g = 'norm')
#' comp.param <- list(f = NULL, g = list(mean = 0, sd = 1))
#' BVdk_contrast(c(0.3,2), data1, density(data1)$bw, comp.dist, comp.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

BVdk_contrast <- function(param, data, h, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")

  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("p", comp.dist[[2]])
  #  comp_BVdk <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_BVdk <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  #for (i in 1:length(comp_BVdk)) assign(x = names(comp_BVdk)[i], value = comp_BVdk[[i]])
  assign(x = names(comp_BVdk)[1], value = comp_BVdk[[1]])

  ## Creates the expression allowing further to generate the right data:
  expr1 <- paste(names(comp_BVdk)[1],"(data[i] + mu,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  expr2 <- paste(names(comp_BVdk)[1],"(mu - data[i],", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")

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
#' @param h The window width used in the kernel estimations.
#' @param comp.dist A list with two elements corresponding to component distributions (specified with R native names for these distributions) involved
#'                  in the admixture model. Unknown elements must be specified as 'NULL' objects, e.g. when 'f' is unknown: list(f=NULL, g='norm').
#' @param comp.param A list with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Unknown elements must be specified as 'NULL' objects, e.g. if 'f' is unknown: list(f=NULL, g=list(mean=0,sd=1)).
#'
#' @references
#' \insertRef{BordesVandekerkhove2010}{admix}
#'
#' @return A numeric vector composed of the two partial derivatives w.r.t. the two parameters on which to optimize the contrast.
#'
#' @examples
#' ## Simulate data:
#' comp.dist <- list(f = 'norm', g = 'norm')
#' comp.param <- list(f = list(mean = 3, sd = 0.5), g = list(mean = 0, sd = 1))
#' data1 <- rsimmix(n = 1000, unknownComp_weight = 0.6, comp.dist, comp.param)[['mixt.data']]
#' ## Compute the contrast gradient for some given parameter vector in real-life framework:
#' comp.dist <- list(f = NULL, g = 'norm')
#' comp.param <- list(f = NULL, g = list(mean = 0, sd = 1))
#' BVdk_contrast_gradient(c(0.3,2), data1, density(data1)$bw, comp.dist, comp.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

BVdk_contrast_gradient <- function(param, data, h, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")

  ## Extracts the information on component distributions and stores in expressions:
  comp.dist.cdf <- paste0("p", comp.dist[[2]])
  #  comp_cdf <- sapply(X = comp.dist.cdf, FUN = get, pos = "package:stats", mode = "function")
  comp_cdf <- sapply(X = comp.dist.cdf, FUN = get, mode = "function")
  assign(x = names(comp_cdf)[1], value = comp_cdf[[1]])
  ## Creates the expression allowing further to generate the right data:
  expr1.cdf <- paste(names(comp_cdf)[1],"(data[i] + mu,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  expr2.cdf <- paste(names(comp_cdf)[1],"(mu - data[i],", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")

  ## Same with density functions :
  comp.dist.dens <- paste0("d", comp.dist[[2]])
  #  comp_dens <- sapply(X = comp.dist.dens, FUN = get, pos = "package:stats", mode = "function")
  comp_dens <- sapply(X = comp.dist.dens, FUN = get, mode = "function")
  assign(x = names(comp_dens)[1], value = comp_dens[[1]])
  ## Creates the expression allowing further to generate the right data:
  expr1.dens <- paste(names(comp_dens)[1],"(data[i] + mu,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  expr2.dens <- paste(names(comp_dens)[1],"(mu - data[i],", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")

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
