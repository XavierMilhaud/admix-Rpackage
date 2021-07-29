#' Maximum Likelihood estimation of the variance of the unknown density variance estimator in an admixture model
#'
#' Parametric estimation of the variance of the variance parameter in Bordes & Vandekerkhove (2010) setting, i.e. considering the
#' admixture model with probability density function (pdf) l:
#'          l(x) = p*f(x-mu) + (1-p)*g,
#' where g is the known component of the two-component mixture, p is the mixture proportion, f is the unknown component with
#' symmetric density, and mu is the location shift parameter. The estimation of the variance of the variance related to the density
#' f is made by maximum likelihood optimization through the information matrix, with the assumption that the unknown f is gaussian.
#'
#' @param data The observed sample under study.
#' @param hat_w Estimate of the unknown component weight.
#' @param hat_loc Estimate of the location shift parameter.
#' @param hat_var Estimate of the variance of the symmetric density f, obtained by plugging-in the previous estimates. See 'Details'
#'                below for further information.
#' @param comp.dist A list with two elements corresponding to component distributions (specified with R native names for these distributions) involved
#'                  in the admixture model. Unknown elements must be specified as 'NULL' objects, e.g. when 'f' is unknown: list(f=NULL, g='norm').
#' @param comp.param A list with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Unknown elements must be specified as 'NULL' objects, e.g. if 'f' is unknown: list(f=NULL, g=list(mean=0,sd=1)).
#'
#' @details Plug-in strategy is defined in Pommeret, D. and Vandekerkhove, P. (2019); Semiparametric density testing in the contamination model;
#'          Electronic Journal of Statistics, 13, pp. 4743--4793. The variance of the estimator variance of the unknown density f
#'          is needed in a testing perspective, since included in the variance of the test statistic. Other details about the information
#'          matrix can be found in Bordes, L. and Vandekerkhove, P. (2010); Semiparametric two-component mixture model when a component is known:
#'          an asymptotically normal estimator; Math. Meth. Stat.; 19, pp. 22--41.
#'
#' @return The variance of the estimator of the variance of the unknown component density f.
#'
#' @examples
#' \donttest{
#' ## Simulate data:
#' list.comp <- list(f = "norm", g = "norm")
#' list.param <- list(f = c(mean = 4, sd = 1), g = c(mean = 7, sd = 0.5))
#' sim.data <- rsimmix(n = 400, unknownComp_weight = 0.9, list.comp, list.param)$mixt.data
#' ## Estimate mixture weight and location shift parameters in real-life:
#' list.comp <- list(f = NULL, g = "norm")
#' list.param <- list(f = NULL, g = c(mean = 7, sd = 0.5))
#' estim <- BVdk_estimParam(data = sim.data, method = "L-BFGS-B",
#'                          comp.dist = list.comp, comp.param = list.param)
#' ## Estimation of the second-order moment of the known component distribution:
#' m2_knownComp <- mean(rnorm(n = 1000000, mean = 7, sd = 0.5)^2)
#' hat_s2 <- (1/estim[1]) * (mean(sim.data^2) - ((1-estim[1])*m2_knownComp)) - estim[2]^2
#' ## Estimated variance of variance estimator related to the unknown symmetric component density:
#' BVdk_ML_varCov_estimators(data = sim.data, hat_w = estim[1], hat_loc = estim[2],
#'                           hat_var = hat_s2, comp.dist = list.comp, comp.param = list.param)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

BVdk_ML_varCov_estimators <- function(data, hat_w, hat_loc, hat_var, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")
  if (is.null(comp.dist[[1]])) comp.dist[[1]] <- 'norm'

  ## Extracts the information on component distributions and stores in expressions:
  comp.dist.dens <- paste0("d", comp.dist)
#  comp_dens <- sapply(X = comp.dist.dens, FUN = get, pos = "package:stats", mode = "function")
  comp_dens <- sapply(X = comp.dist.dens, FUN = get, mode = "function")
  for (i in 1:length(comp_dens)) assign(x = names(comp_dens)[i], value = comp_dens[[i]])
  ## Creates the expression allowing further to compute the hessian (the first component is the unknown one 'f', whereas the second one is 'g') :
  expr1 <- paste(names(comp_dens)[1],"(y, hat_loc, sqrt(hat_var))", sep = "")
  expr2 <- paste(names(comp_dens)[2],"(y,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")

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

