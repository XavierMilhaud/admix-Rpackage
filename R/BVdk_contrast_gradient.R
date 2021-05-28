#' Gradient of the contrast as defined in Bordes & Vandekerkhove (2010)
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
#' @details The contrast is defined in Bordes, L. and Vandekerkhove, P. (2010); Semiparametric two-component mixture model
#'          when a component is known: an asymptotically normal estimator; Math. Meth. Stat.; 19, pp. 22--41.
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

BVdk_contrast_gradient <- function(param, data, h, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")

  ## Extracts the information on component distributions and stores in expressions:
  comp.dist.cdf <- paste0("p", comp.dist[[2]])
  comp_cdf <- sapply(X = comp.dist.cdf, FUN = get, pos = "package:stats", mode = "function")
  assign(x = names(comp_cdf)[1], value = comp_cdf[[1]])
  ## Creates the expression allowing further to generate the right data:
  expr1.cdf <- paste(names(comp_cdf)[1],"(data[i] + mu,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  expr2.cdf <- paste(names(comp_cdf)[1],"(mu - data[i],", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")

  ## Same with density functions :
  comp.dist.dens <- paste0("d", comp.dist[[2]])
  comp_dens <- sapply(X = comp.dist.dens, FUN = get, pos = "package:stats", mode = "function")
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

