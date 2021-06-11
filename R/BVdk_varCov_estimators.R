#' Estimation of the variance of the estimators in admixture models with symmetric unknown density
#'
#' Semiparametric estimation of the variance of the estimators, i.e. the mixture weight p and the location shift parameter mu
#' considering the admixture model with probability density function l:
#'          l(x) = p*f(x-mu) + (1-p)*g(x), x in R,
#' where g is the known component of the two-component mixture, p is the unknown proportion, f is the unknown component density and
#' mu is the location shift. See 'Details' below for more information.
#'
#' @param data The observed sample under study.
#' @param loc The estimated location shift parameter, related to the unknown symmetric density.
#' @param p The estimated unknown component weight.
#' @param comp.dist A list with two elements corresponding to component distributions (specified with R native names for these distributions) involved
#'                  in the admixture model. Unknown elements must be specified as 'NULL' objects, e.g. when 'f' is unknown: list(f=NULL, g='norm').
#' @param comp.param A list with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Unknown elements must be specified as 'NULL' objects, e.g. if 'f' is unknown: list(f=NULL, g=list(mean=0,sd=1)).
#'
#' @details See formulas pp.28--30 in Appendix of Bordes, L. and Vandekerkhove, P. (2010); Semiparametric two-component mixture
#'          model when a component is known: an asymptotically normal estimator; Math. Meth. Stat.; 19, pp. 22--41.
#'
#' @return A list containing 1) the variance-covariance matrix of the estimators (assessed at the specific time points 'u' and 'v'
#'         such that u=v=mean(data)); 2) the variance of the mixture weight estimator; 3) the variance
#'         of the location shift estimator; 4) the variance of the unknown component cumulative distribution
#'         function at points 'u' and 'v' (useless for most of applications, explaining why 'u' and 'v'
#'         are set equal to mean(data) by default, with no corresponding arguments here).
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f = 'norm', g = 'norm')
#' list.param <- list(f = c(mean = 4, sd = 1), g = c(mean = 7, sd = 0.5))
#' sim.data <- rsimmix(n=140, unknownComp_weight=0.9, comp.dist=list.comp, comp.param=list.param)
#' ## Estimate the location shift and mixture weight parameters in real-life setting:
#' list.comp <- list(f = NULL, g = 'norm')
#' list.param <- list(f = NULL, g = c(mean = 7, sd = 0.5))
#' estimators <- BVdk_estimParam(data = sim.data[['mixt.data']], method = "L-BFGS-B",
#'                               comp.dist = list.comp, comp.param = list.param)
#' ## Estimate the variance of the two estimators (first mixture weight, then location shift):
#' BVdk_varCov_estimators(data = sim.data[['mixt.data']], loc = estimators[2], p = estimators[1],
#'                        comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

BVdk_varCov_estimators <- function(data, loc, p, comp.dist, comp.param)
{
  n <- length(data)
  bandw <- stats::density(data)$bw
  u <- v <- round(mean(data))

  inv_J <- solve(BVdk_mat_J(data = data, loc = loc, p = p, h = bandw, comp.dist = comp.dist, comp.param = comp.param))
  mat_L.u <- BVdk_mat_L(u = u, data = data, loc = loc, p = p, h = bandw, comp.dist, comp.param)
  mat_L.v <- BVdk_mat_L(u = v, data = data, loc = loc, p = p, h = bandw, comp.dist, comp.param)
  sigma.mat.uv <- BVdk_mat_Sigma(u = u, v = v, data = data, loc = loc, p = p, h = bandw, comp.dist, comp.param)

  varCov.matrix <- (1/n) * (mat_L.u %*% inv_J %*% sigma.mat.uv %*% inv_J %*% t(mat_L.v))
  var_hat_p <- varCov.matrix[1,1]
  var_hat_loc <- varCov.matrix[2,2]
  var_hat_F <- varCov.matrix[3,3]

  return( list(varCov.mat = varCov.matrix, var_pEstim = var_hat_p, var_muEstim = var_hat_loc, var_FEstim = var_hat_F) )
}


## Estimation of matrix Sigma(u,v) from empirical versions (cf p.28 et p.29), with arguments:
## - u the time point at which the first (related to the first parameter) underlying empirical process is looked through.
## - v the time point at which the second (related to the second parameter) underlying empirical process is looked through.
BVdk_mat_Sigma <- function(u, v, data, loc, p, h, comp.dist, comp.param)
{
  n <- length(data)
  aux1 <- aux2 <- rep(0, n)
  aux1 <- sapply(X = data, FUN = h1, loc, p, comp.dist, comp.param)
  aux2 <- sapply(X = data, FUN = h2, data, loc, p, h, comp.dist, comp.param)

  sorted_data <- sort(data)
  ul <- sapply(X = data, FUN = l_fun, u, sorted_data, loc)
  vl <- sapply(X = data, FUN = l_fun, v, sorted_data, loc)

  aux12 <- terms_sigma11 <- terms_sigma22 <- terms_sigma12 <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      aux12[i,j] <- ka(data[i], data[j], sorted_data, loc)
      terms_sigma11[i,j] <- aux1[i] * aux1[j] * aux12[i,j]
      terms_sigma22[i,j] <- aux2[i] * aux2[j] * aux12[i,j]
      terms_sigma12[i,j] <- (aux1[i] * aux2[j] + aux1[j] * aux2[i]) * aux12[i,j]
    }
  }

  sigma_mat <- matrix(NA, nrow = 3, ncol = 3)
  sigma_mat[3,3] <- .Call('_admix_Donsker_correl_cpp', PACKAGE = 'admix', u + loc, v + loc, sorted_data)
  sigma_mat[1,3] <- (2 / (n * p)) * sum(aux1 * vl)
  sigma_mat[3,1] <- (2 / (n * p)) * sum(aux1 * ul)
  sigma_mat[2,3] <- (2 / (n * p)) * sum(aux2 * vl)
  sigma_mat[3,2] <- (2 / (n * p)) * sum(aux2 * ul)
  sigma_mat[1,1] <- (8 / (n * (n-1) * p^2)) * sum(terms_sigma11, na.rm = TRUE)
  sigma_mat[2,2] <- (8 / (n * (n-1) * p^2)) * sum(terms_sigma22, na.rm = TRUE)
  sigma_mat[1,2] <- sigma_mat[2,1] <- (4 / (n * (n-1) * p^2)) * sum(terms_sigma12, na.rm = TRUE)

  return(sigma_mat)
}


## Define the matrix J that is useful for the consistency theorem of estimators (cf p.29-30) .
BVdk_mat_J <- function(data, loc, p, h, comp.dist, comp.param)
{
  n <- length(data)
  aux1 <- aux2 <- rep(0, n)
  aux1 <- sapply(X = data, FUN = h1, loc, p, comp.dist, comp.param)
  aux2 <- sapply(X = data, FUN = h2, data, loc, p, h, comp.dist, comp.param)
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
BVdk_mat_L <- function(u, data, loc, p, h, comp.dist, comp.param)
{
  L <- matrix(0, nrow = 3, ncol = 3)
  L[1,1] <- L[2,2] <- 1
  L[3,1] <- h3(u, data, loc, p, comp.dist, comp.param)
  L[3,2] <- h2(u, data, loc, p, h, comp.dist, comp.param) / 2
  L[3,3] <- 1 / p
  return(L)
}


##### All the following functions are defined p.28 of the paper (see p.29 for empirical versions).

## Function h1, see formula (3.9) p.13. Arguments are :
##  - 'hat_loc' : the estimator of the localization parameter 'mu',
##  - 'hat_p' : the estimator of the component weight, related to the unknown mixture component.
h1 <- function(u, loc, p, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")
  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("p", comp.dist[[2]])
  comp_h1 <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  assign(x = names(comp_h1)[1], value = comp_h1[[1]])
  expr1 <- paste(names(comp_h1)[1],"(q=(loc+u),", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  expr2 <- paste(names(comp_h1)[1],"(q=(loc-u),", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")

  res <- (eval(parse(text = expr1)) + eval(parse(text = expr2)) - 1) / p
  return(res)
}

## Fonction h2, see formula (3.8) p.13 (see also the connexion with formula end p.6)
h2 <- function(u, data, loc, p, h, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")
  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("d", comp.dist[[2]])
  comp_h2 <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  assign(x = names(comp_h2)[1], value = comp_h2[[1]])
  expr <- paste(names(comp_h2)[1],"(x=(loc+u),", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")

  g  <- mean( kernel_density(u + loc - data, h) )
  fo <- eval(parse(text = expr))
  f  <- (g - (1-p) * fo) / p			# cf (2.2) p.5

  res <- 2 * f * (f >= 0)         # end p.28
  return(res)
}

## Fonction h3:
h3 <- function(u, data, loc, p, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")
  ## Extracts the information on component distributions and stores in expressions:
  exp.comp.dist <- paste0("p", comp.dist[[2]])
  comp_h3 <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  assign(x = names(comp_h3)[1], value = comp_h3[[1]])
  expr <- paste(names(comp_h3)[1],"(q=(loc+u),", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  G_ecdf <- stats::ecdf(data)

  res <- (eval(parse(text = expr)) - G_ecdf(u + loc)) / p^2
  return(res)
}

## Function ka of the paper:
ka <- function(u, v, data, loc)
{
  return( l_fun(u, v, data, loc) + l_fun(-u, v, data, loc) )
}

## Old implementation of fonction k in the article:
ka_old <- function(u, v, data, loc)
{
  return( l_fun_old(u, v, data, loc) + l_fun_old(-u, v, data, loc) )
}

## Function l of the article:
l_fun <- function(u, v, data, loc)
{
  term1 <- .Call('_admix_Donsker_correl_cpp', PACKAGE = 'admix', loc + u, loc + v, data)
  term2 <- .Call('_admix_Donsker_correl_cpp', PACKAGE = 'admix', loc + u, loc - v, data)
  return(term1 + term2)
}

## Old implementation of fonction l in the article:
l_fun_old <- function(u, v, data, loc)
{
  return( Donsker_correl_old(loc + u, loc + v, data) + Donsker_correl_old(loc + u, loc - v, data) )
}

## Donsker correlation:
Donsker_correl_old <- function(u, v, obs.data)
{
  L.CDF <- stats::ecdf(obs.data)
  return( L.CDF(min(u,v)) * (1 - L.CDF(max(u,v))) )
}
