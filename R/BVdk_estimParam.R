#' Estimation of the parameters in a two-component admixture model with symmetric unknown density
#'
#' Estimation of the two parameters (mixture weight as well as location shift) in the admixture model with pdf:
#'          l(x) = p*f(x-mu) + (1-p)*g(x), x in R,
#' where g is the known component, p is the proportion and f is the unknown component with symmetric density.
#' The localization shift parameter is thus denoted mu, and the component weight p.
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
#' @details Parameters are estimated by minimization of the contrast function, where the contrast is defined in
#'          Bordes, L. and Vandekerkhove, P. (2010); Semiparametric two-component mixture model when a component is known:
#'          an asymptotically normal estimator; Math. Meth. Stat.; 19, pp. 22--41.
#'
#' @return A numeric vector with the two estimated parameters (proportion first, and then location shift).
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f = 'norm', g = 'norm')
#' list.param <- list(f = list(mean = 3, sd = 0.5),
#'                    g = list(mean = 0, sd = 1))
#' data1 <- rsimmix(n = 150, unknownComp_weight = 0.9, list.comp, list.param)[['mixt.data']]
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
