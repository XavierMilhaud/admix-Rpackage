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
#' @details The contrast is defined in Bordes, L. and Vandekerkhove, P. (2010); Semiparametric two-component mixture model
#'          when a component is known: an asymptotically normal estimator; Math. Meth. Stat.; 19, pp. 22--41.
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

BVdk_contrast <- function(param, data, h, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")

  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("p", comp.dist[[2]])
  comp_BVdk <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
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
