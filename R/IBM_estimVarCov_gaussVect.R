#' Nonparametric estimation of the variance-covariance matrix of the gaussian vector in IBM approach
#'
#' Estimate the variance-covariance matrix of the gaussian vector at point 'z', considering the use of Inversion - Best
#' Matching (IBM) method to estimate the model parameters in two-sample admixture models.
#' Recall that the two admixture models have respective probability density functions (pdf) l1 and l2, such that:
#'   l1 = p1*f1 + (1-p1)*g1 and l2 = p2*f2 + (1-p2)*g2, where g1 and g2 are the known component densities.
#' Further information for the IBM approach are given in 'Details' below.
#'
#' @param x Time point at which the first (related to the first parameter) underlying empirical process is looked through.
#' @param y Time point at which the second (related to the second parameter) underlying empirical process is looked through.
#' @param estim.obj Object obtained from the estimation of the component weights related to the proportions of the
#'                  unknown component in each of the two admixture models.
#' @param fixed.p1 Arbitrary value chosen by the user for the component weight related to the unknown component of the first
#'                  admixture model. Only useful for optimization when the known components of the two models are identical
#'                  (G1=G2, leading to unidimensional optimization).
#' @param known.p (optional, NULL by default) Numeric vector with two elements, the known (true) mixture weights.
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param min_size (optional, NULL by default) in the k-sample case, useful to provide the minimal size among all samples
#'                 (needed to take into account the correction factor in variance-covariance assessment). Otherwise, useless.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#'
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03201760
#'
#' @return The estimated variance-covariance matrix of the gaussian vector Z = (hat(p1),(hat(p2),Dn(z)), at location '(x,y)'.
#'
#' @examples
#' \donttest{
#' ######## Analysis by simulated data:
#' ## Simulate Gamma - Exponential admixtures :
#' list.comp <- list(f1 = "gamma", g1 = "exp",
#'                   f2 = "gamma", g2 = "exp")
#' list.param <- list(f1 = list(shape = 2, scale = 3), g1 = list(rate = 1/3),
#'                    f2 = list(shape = 2, scale = 3), g2 = list(rate = 1/5))
#' X.sim <- rsimmix(n=400, unknownComp_weight=0.8, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' Y.sim <- rsimmix(n=350, unknownComp_weight=0.9, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' ## Real-life setting:
#' list.comp <- list(f1 = NULL, g1 = "exp",
#'                   f2 = NULL, g2 = "exp")
#' list.param <- list(f1 = NULL, g1 = list(rate = 1/3),
#'                    f2 = NULL, g2 = list(rate = 1/5))
#' ## Estimate the unknown component weights in the two admixture models:
#' estim <- IBM_estimProp(sample1 =X.sim, sample2 =Y.sim, known.prop = NULL, comp.dist = list.comp,
#'                        comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' IBM_estimVarCov_gaussVect(x = mean(X.sim), y = mean(Y.sim), estim.obj = estim,
#'                           fixed.p1 = estim[["p.X.fixed"]], known.p = NULL, sample1=X.sim,
#'                           sample2 = Y.sim, min_size = NULL,
#'                           comp.dist = list.comp, comp.param = list.param)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_estimVarCov_gaussVect <- function(x, y, estim.obj, fixed.p1 = NULL, known.p = NULL, sample1, sample2, min_size = NULL,
                                     comp.dist = NULL, comp.param = NULL)
{
  if (!is.null(known.p)) {
    estimators <- estim.obj[["theo.prop.estim"]]
  } else {
    estimators <- estim.obj[["prop.estim"]]
  }
  ## Location at which the variance-covariance matrix is evaluated:
  varCov.gaussianVect <- IBM_mat_Sigma(x = x, y = y, par = estimators, fixed.p1 = fixed.p1, known.p = known.p, sample1 = sample1,
                                       sample2 = sample2, G = estim.obj[["integ.supp"]], comp.dist = comp.dist, comp.param = comp.param)
  ## Adjusts by the normalization factor to get the distribution of the gaussian vector Z=(hat(p1), hat(p2), Dn(z)):
  normalization.factor <- IBM_normalization_term(z = x, estim.obj = estim.obj, fixed.p1 = fixed.p1, known.p = known.p, sample1 = sample1,
                                                 sample2 = sample2, min_size = min_size, comp.dist = comp.dist, comp.param = comp.param)
  varCov.transfo_gaussVect <- normalization.factor %*% varCov.gaussianVect %*% t(normalization.factor)

  return(varCov.transfo_gaussVect)
}


## Normalization of the variance-covariance by matrix 'M' in the paper
IBM_normalization_term <- function(z, estim.obj, fixed.p1 = NULL, known.p = NULL, sample1, sample2, min_size = NULL, comp.dist = NULL, comp.param = NULL)
{
  n1 <- length(sample1)
  n2 <- length(sample2)
  if (is.null(min_size)) {
    xi <- 1 / sqrt(max(n1,n2)/min(n1,n2))
  } else {
    xi1 <- 1 / sqrt(n1 / min_size)
    xi2 <- 1 / sqrt(n2 / min_size)
  }

  if (!is.null(known.p)) {
    estimators <- estim.obj[["theo.prop.estim"]]
  } else {
    estimators <- estim.obj[["prop.estim"]]
  }
  ## Adjusts by the normalization factor to get the distribution of the gaussian vector Z=(hat(p1), hat(p2), Dn(z)):
  L <- IBM_mat_L(z = z, par = estimators, fixed.p1 = fixed.p1, known.p = known.p, sample1 = sample1, sample2 = sample2,
             min_size = min_size, comp.dist = comp.dist, comp.param = comp.param)
  inv_J <- solve(IBM_mat_J(par = estimators, fixed.p1 = fixed.p1, known.p = known.p, sample1 = sample1, sample2 = sample2,
                       G = estim.obj[["integ.supp"]], comp.dist = comp.dist, comp.param = comp.param))

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)
  if (G1equalG2) {
    if (is.null(min_size)) {
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
      C <- matrix(c(-xi1,0,-xi2,0,
                    0,1,0,0,
                    0,0,0,1), nrow = 3, ncol = 4, byrow = T)
    }

  } else {
    if (is.null(min_size)) {
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
    } else {
      C <- matrix(c(-xi1,0,0,0,-xi2,0,
                    0,-xi1,0,-xi2,0,0,
                    0,0,1,0,0,0,
                    0,0,0,0,0,1), nrow = 4, ncol = 6, byrow = T)
    }
  }

  return(L %*% inv_J %*% C)
}


## Matrice J de l'article (normalisation par la variance):
IBM_mat_J <- function(par, fixed.p1 = NULL, known.p = NULL, sample1, sample2, G, comp.dist = NULL, comp.param = NULL)
{
  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  if (G1equalG2) {
    stopifnot(!is.null(fixed.p1))
    J <- diag(1, nrow = 3, ncol = 3)
    J[1,1] <- IBM_hessian_contrast(par = par, fixed.p1 = fixed.p1, known.p = known.p, sample1 = sample1, sample2 = sample2, G = G,
                                   comp.dist = comp.dist, comp.param = comp.param)
  } else {
    J <- diag(1, nrow = 4, ncol = 4)
    J[1:2,1:2] <- IBM_hessian_contrast(par = par, fixed.p1 = fixed.p1, known.p = known.p, sample1 = sample1, sample2 = sample2, G = G,
                                       comp.dist = comp.dist, comp.param = comp.param)
  }
  return(J)
}

## Matrice L de l'article (matrice muette pour les estimateurs de p1 et p2):
IBM_mat_L <- function(z, par, fixed.p1 = NULL, known.p = NULL, sample1, sample2, min_size = NULL, comp.dist, comp.param)
{
  n1 <- length(sample1)
  n2 <- length(sample2)
  if (is.null(min_size)) {
    xi <- 1 / sqrt(max(n1,n2)/min(n1,n2))
  } else {
    xi1 <- 1 / sqrt(n1 / min_size)
    xi2 <- 1 / sqrt(n2 / min_size)
  }

  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
    stop("Known components must be specified in the two admixture models.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
  }
  ## Extracts the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
#  comp_L <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_L <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_L)) assign(x = names(comp_L)[i], value = comp_L[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_L)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                  paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_L)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Differentiates when the 4 components were specified (simulations) and real-life cases (only known components are given).
  ## First, real-life case:
  if (all(unlist(sapply(comp.param, is.na)[c(1,3)]))) {
    if (any(exp.comp.dist == "stepfun")) {
      G1.fun <- eval(parse(text = expr[2]))
      G2.fun <- eval(parse(text = expr[4]))
      G1 <- function(z) G1.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      G1 <- function(z) { eval(parse(text = expr[2])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  } else {
    ## case of simulated data:
    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      G1.fun <- eval(parse(text = expr[2]))
      F2.fun <- eval(parse(text = expr[3]))
      G2.fun <- eval(parse(text = expr[4]))
      F1 <- function(z) F1.fun(z)
      G1 <- function(z) G1.fun(z)
      F2 <- function(z) F2.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
      F2 <- function(z) { eval(parse(text = expr[3])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  }

  ## Cumulative distribution functions (empirical or exact) at point 'z':
  if ( !is.null(known.p) ) {
    L1 <- function(z) { known.p[1] * F1(z) + (1-known.p[1]) * G1(z) }
    L2 <- function(z) { known.p[2] * F2(z) + (1-known.p[2]) * G2(z) }
  } else {
    L1 <- stats::ecdf(sample1)
    L2 <- stats::ecdf(sample2)
  }

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  if (G1equalG2) {
    stopifnot(!is.null(fixed.p1))
    L <- matrix(0, nrow = 2, ncol = 3)
    if (is.null(min_size)) {
      if (n1 <= n2) {
        L[1,1] <- 1
        L[2,1] <- (1/par^2) * (L2(z) - G2(z))
        L[2,2] <- 1 / fixed.p1
        L[2,3] <- -xi / par
      } else {
        L[1,1] <- 1
        L[2,1] <- (1/par^2) * (L2(z) - G2(z))
        L[2,2] <- xi / fixed.p1
        L[2,3] <- -1 / par
      }
    } else {
      L[1,1] <- 1
      L[2,1] <- (1/par^2) * (L2(z) - G2(z))
      L[2,2] <- xi1 / fixed.p1
      L[2,3] <- -xi2 / par
    }
  } else {
    L <- matrix(0, nrow = 3, ncol = 4)
    if (is.null(min_size)) {
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
    } else {
      L[1,1] <- L[2,2] <- 1
      L[3,1] <- -(1/par[1]^2) * (L1(z) - G1(z))
      L[3,2] <- (1/par[2]^2) * (L2(z) - G2(z))
      L[3,3] <- xi1 / par[1]
      L[3,4] <- -xi2 / par[2]
    }
  }

  return(L)
}


## Variance-covariance matrix of the Gaussian random vector. Arguments are defined as follows:
## This function returns the non-normalized variance-covariance matrix of the gaussian random vector.
IBM_mat_Sigma <- function(x, y, par, fixed.p1 = NULL, known.p = NULL, sample1, sample2, G, comp.dist = NULL, comp.param = NULL)
{
  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  if (G1equalG2) {
    Sigma <- matrix(0, nrow = 4, ncol = 4)
    ## Upper block of the matrix :
    Sigma[1:2,1:2] <- IBM_Sigma1(x = x, y = y, par = par, fixed.p1 = fixed.p1, known.p = known.p, sample1 = sample1, sample2 = sample2, G = G,
                                 comp.dist = comp.dist, comp.param = comp.param)
    ## Lower block of the matrix :
    Sigma[3:4,3:4] <- IBM_Sigma2(x = x, y = y, par = par, fixed.p1 = fixed.p1, known.p = known.p, sample1 = sample1, sample2 = sample2, G = G,
                                 comp.dist = comp.dist, comp.param = comp.param)
  } else {
    Sigma <- matrix(0, nrow = 6, ncol = 6)
    ## Upper block of the matrix :
    Sigma[1:3,1:3] <- IBM_Sigma1(x = x, y = y, par = par, fixed.p1 = NULL, known.p = known.p, sample1 = sample1, sample2 = sample2, G = G,
                                 comp.dist = comp.dist, comp.param = comp.param)
    ## Lower block of the matrix :
    Sigma[4:6,4:6] <- IBM_Sigma2(x = x, y = y, par = par, fixed.p1 = NULL, known.p = known.p, sample1 = sample1, sample2 = sample2, G = G,
                                 comp.dist = comp.dist, comp.param = comp.param)
  }
  return(Sigma)
}


## Upper block of the variance-covariance matrix:
IBM_Sigma1 <- function(x, y, par, fixed.p1 = NULL, known.p = NULL, sample1, sample2, G, comp.dist = NULL, comp.param = NULL)
{
  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
    stop("Known components must be specified in the two admixture models.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
  }
  ## Extracts the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
#  comp_sigma1 <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_sigma1 <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_sigma1)) assign(x = names(comp_sigma1)[i], value = comp_sigma1[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_sigma1)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                                     paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_sigma1)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Differentiates when the 4 components were specified (simulations) and real-life cases (only known components are given).
  ## First, real-life case:
  if (all(unlist(sapply(comp.param, is.na)[c(1,3)]))) {
    if (any(exp.comp.dist == "stepfun")) {
      G1.fun <- eval(parse(text = expr[2]))
      G2.fun <- eval(parse(text = expr[4]))
      G1 <- function(z) G1.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      G1 <- function(z) { eval(parse(text = expr[2])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  } else {
    ## case of simulated data:
    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      G1.fun <- eval(parse(text = expr[2]))
      F2.fun <- eval(parse(text = expr[3]))
      G2.fun <- eval(parse(text = expr[4]))
      F1 <- function(z) F1.fun(z)
      G1 <- function(z) G1.fun(z)
      F2 <- function(z) F2.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
      F2 <- function(z) { eval(parse(text = expr[3])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  }

  ## Empirical cumulative distribution functions:
  if (!is.null(known.p)) {
    L1 <- function(z) { known.p[1] * F1(z) + (1-known.p[1]) * G1(z) }
    L2 <- function(z) { known.p[2] * F2(z) + (1-known.p[2]) * G2(z) }
  } else {
    L1 <- stats::ecdf(sample1)
    L2 <- stats::ecdf(sample2)
  }

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  ## Loi de G et support d'integration :
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  ## Numerical computation of each term included in the variance-covariance matrix of the gaussian vector:
  if (G1equalG2) {
    stopifnot(length(par) == 1)

    psi1 <- function(z) 2*( ((2-par)/par^3) * G1(z) - (2/par^3)*L2(z) + (1/(par^2*fixed.p1))*L1(z) - ((1-fixed.p1)/(par^2*fixed.p1)) * G1(z) )
    psi2 <- function(z) 2*( (1/(par^2*fixed.p1)) * (L2(z) - G1(z)) )
    integrand.sigma1.11_g1eqg2 <- function(x,y) {
      if (support == "continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi2(x) * psi2(y) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma1.12_g1eqg2 <- function(x) {
      if (support == "continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
      }
      psi2(x) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.x
    }
    integrand.sigma1.21_g1eqg2 <- function(y) {
      if (support == "continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi2(y) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.y
    }

    sigma <- matrix(NA, nrow = 2, ncol = 2)
      if (support == "continuous") {
        sigma[1,1] <- pracma::integral2(Vectorize(integrand.sigma1.11_g1eqg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
        sigma[1,2] <- stats::integrate(Vectorize(integrand.sigma1.12_g1eqg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
        sigma[2,1] <- stats::integrate(Vectorize(integrand.sigma1.21_g1eqg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
        sigma[2,2] <- estimVarCov_empProcess(x, y, sample1, known.p[1], list(comp.dist$f1,comp.dist$g1), list(comp.param$f1,comp.param$g1))
      } else {
        sigma[1,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma1.11_g1eqg2, x, supp.integration)) )
        sigma[1,2] <- sum( unlist(sapply(supp.integration, integrand.sigma1.12_g1eqg2)) )
        sigma[2,1] <- sum( unlist(sapply(supp.integration, integrand.sigma1.21_g1eqg2)) )
        sigma[2,2] <- estimVarCov_empProcess(x, y, sample1, known.p[1], list(comp.dist$f1,comp.dist$g1), list(comp.param$f1,comp.param$g1))
      }

  } else {

    psi1.1 <- function(z) 2*( ((2-par[1])/par[1]^3) * G1(z) - (2/par[1]^3)*L1(z) + (1/(par[1]^2*par[2]))*L2(z) - ((1-par[2])/(par[1]^2*par[2])) * G2(z) )
    psi1.2 <- function(z) 2*( (1/(par[1]^2*par[2])) * (L1(z) - G1(z)) )
    psi2.1 <- function(z) 2*( ((2-par[2])/par[2]^3) * G2(z) - (2/par[2]^3)*L2(z) + (1/(par[2]^2*par[1]))*L1(z) - ((1-par[1])/(par[2]^2*par[1])) * G1(z) )
    psi2.2 <- function(z) 2*( (1/(par[2]^2*par[1])) * (L2(z) - G2(z)) )
    integrand.sigma1.11_g1diffg2 <- function(x,y) {
      if (support == "continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi1.1(x) * psi1.1(y) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma1.12_g1diffg2 <- function(x,y) {
      if (support == "continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi1.1(x) * psi2.2(y) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma1.13_g1diffg2 <- function(x) {
      if (support == "continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
      }
      psi1.1(x) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.x
    }
    integrand.sigma1.22_g1diffg2 <- function(x,y) {
      if (support == "continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi2.2(x) * psi2.2(y) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma1.23_g1diffg2 <- function(x) {
      if (support == "continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
      }
      psi2.2(x) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.x
    }
    integrand.sigma1.31_g1diffg2 <- function(y) {
      if (support == "continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi1.1(y) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.y
    }
    integrand.sigma1.32_g1diffg2 <- function(y) {
      if (support == "continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi2.2(y) * estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1)) * densite.G.dataPoint.y
    }

    sigma <- matrix(NA, nrow = 3, ncol = 3)
      if (support == "continuous") {
        sigma[1,1] <- pracma::integral2(Vectorize(integrand.sigma1.11_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
        sigma[1,2] <- sigma[2,1] <- pracma::integral2(Vectorize(integrand.sigma1.12_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
        sigma[1,3] <- stats::integrate(Vectorize(integrand.sigma1.13_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
        sigma[2,2] <- pracma::integral2(Vectorize(integrand.sigma1.22_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
        sigma[2,3] <- stats::integrate(Vectorize(integrand.sigma1.23_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
        sigma[3,1] <- stats::integrate(Vectorize(integrand.sigma1.31_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
        sigma[3,2] <- stats::integrate(Vectorize(integrand.sigma1.32_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
        sigma[3,3] <- estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1))
      } else {
        sigma[1,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma1.11_g1diffg2, x, supp.integration)) )
        sigma[1,2] <- sigma[2,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma1.12_g1diffg2, x, supp.integration)) )
        sigma[1,3] <- sum( unlist(sapply(supp.integration, integrand.sigma1.13_g1diffg2)) )
        sigma[2,2] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma1.22_g1diffg2, x, supp.integration)) )
        sigma[2,3] <- sum( unlist(sapply(supp.integration, integrand.sigma1.23_g1diffg2)) )
        sigma[3,1] <- sum( unlist(sapply(supp.integration, integrand.sigma1.31_g1diffg2)) )
        sigma[3,2] <- sum( unlist(sapply(supp.integration, integrand.sigma1.32_g1diffg2)) )
        sigma[3,3] <- estimVarCov_empProcess(x,y,sample1,known.p[1],list(comp.dist$f1,comp.dist$g1),list(comp.param$f1,comp.param$g1))
      }
  }

  return(sigma)
}


## Lower block:
IBM_Sigma2 <- function(x, y, par, fixed.p1 = NULL, known.p = NULL, sample1, sample2, G, comp.dist = NULL, comp.param = NULL)
{
  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
    stop("Known components must be specified in the two admixture models.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
  }
  ## Extracts the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
#  comp_sigma2 <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_sigma2 <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_sigma2)) assign(x = names(comp_sigma2)[i], value = comp_sigma2[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_sigma2)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                      paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_sigma2)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Differentiates when the 4 components were specified (simulations) and real-life cases (only known components are given).
  ## First, real-life case:
  if (all(unlist(sapply(comp.param, is.na)[c(1,3)]))) {
    if (any(exp.comp.dist == "stepfun")) {
      G1.fun <- eval(parse(text = expr[2]))
      G2.fun <- eval(parse(text = expr[4]))
      G1 <- function(z) G1.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      G1 <- function(z) { eval(parse(text = expr[2])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  } else {
    ## case of simulated data:
    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      G1.fun <- eval(parse(text = expr[2]))
      F2.fun <- eval(parse(text = expr[3]))
      G2.fun <- eval(parse(text = expr[4]))
      F1 <- function(z) F1.fun(z)
      G1 <- function(z) G1.fun(z)
      F2 <- function(z) F2.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
      F2 <- function(z) { eval(parse(text = expr[3])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  }

  ## Cumulative distribution function (either empirical from the two observed samples, or theoretical):
  if ( !is.null(known.p) ) {
    L1 <- function(z) { known.p[1] * F1(z) + (1-known.p[1]) * G1(z) }
    L2 <- function(z) { known.p[2] * F2(z) + (1-known.p[2]) * G2(z) }
  } else {
    L1 <- stats::ecdf(sample1)
    L2 <- stats::ecdf(sample2)
  }

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  ## Loi de G et support d'integration:
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  ## Numerical computation of each term included in the variance-covariance matrix of the gaussian vector:
  if (G1equalG2) {

    stopifnot(length(par) == 1)
    psi1 <- function(z) 2*( ((2-par)/par^3) * G1(z) - (2/par^3)*L2(z) + (1/(par^2*fixed.p1))*L1(z) - ((1-fixed.p1)/(par^2*fixed.p1)) * G1(z) )
    psi2 <- function(z) 2*( (1/(par^2*fixed.p1)) * (L2(z) - G1(z)) )

    integrand.sigma2.11_g1eqg2 <- function(x,y) {
      if (support == "continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi1(x) * psi1(y) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma2.12_g1eqg2 <- function(x) {
      if (support == "continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
      }
      psi1(x) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.x
    }
    integrand.sigma2.21_g1eqg2 <- function(y) {
      if (support == "continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi1(y) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.y
    }

    sigma <- matrix(NA, nrow = 2, ncol = 2)
    if (support == "continuous") {
      sigma[1,1] <- pracma::integral2(Vectorize(integrand.sigma2.11_g1eqg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,2] <- stats::integrate(Vectorize(integrand.sigma2.12_g1eqg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,1] <- stats::integrate(Vectorize(integrand.sigma2.21_g1eqg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,2] <- estimVarCov_empProcess(x, y, sample2, known.p[2], list(comp.dist$f2,comp.dist$g2), list(comp.param$f2,comp.param$g2))
    } else {
      sigma[1,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma2.11_g1eqg2, x, supp.integration)) )
      sigma[1,2] <- sum( unlist(sapply(supp.integration, integrand.sigma2.12_g1eqg2)) )
      sigma[2,1] <- sum( unlist(sapply(supp.integration, integrand.sigma2.21_g1eqg2)) )
      sigma[2,2] <- estimVarCov_empProcess(x, y, sample2, known.p[2], list(comp.dist$f2,comp.dist$g2), list(comp.param$f2,comp.param$g2))
    }

  } else {

    psi1.1 <- function(z) 2*( ((2-par[1])/par[1]^3) * G1(z) - (2/par[1]^3)*L1(z) + (1/(par[1]^2*par[2]))*L2(z) - ((1-par[2])/(par[1]^2*par[2])) * G2(z) )
    psi1.2 <- function(z) 2*( (1/(par[1]^2*par[2])) * (L1(z) - G1(z)) )
    psi2.1 <- function(z) 2*( ((2-par[2])/par[2]^3) * G2(z) - (2/par[2]^3)*L2(z) + (1/(par[2]^2*par[1]))*L1(z) - ((1-par[1])/(par[2]^2*par[1])) * G1(z) )
    psi2.2 <- function(z) 2*( (1/(par[2]^2*par[1])) * (L2(z) - G2(z)) )

    integrand.sigma2.11_g1diffg2 <- function(x,y) {
      if (support == "continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi2.1(x) * psi2.1(y) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma2.12_g1diffg2 <- function(x,y) {
      if (support == "continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi2.1(x) * psi1.2(y) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma2.13_g1diffg2 <- function(x) {
      if (support == "continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
      }
      psi2.1(x) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.x
    }
    integrand.sigma2.22_g1diffg2 <- function(x,y) {
      if (support == "continuous") {
        densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
        densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi1.2(x) * psi1.2(y) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.x * densite.G.dataPoint.y
    }
    integrand.sigma2.23_g1diffg2 <- function(x) {
      if (support == "continuous") { densite.G.dataPoint.x <- stats::approx(densite.G$x, densite.G$y, xout = x)$y
      } else {
        densite.G.dataPoint.x <- 1 / length(table(c(sample1,sample2)))
      }
      psi1.2(x) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.x
    }
    integrand.sigma2.31_g1diffg2 <- function(y) {
      if (support == "continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi2.1(y) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.y
    }
    integrand.sigma2.32_g1diffg2 <- function(y) {
      if (support == "continuous") { densite.G.dataPoint.y <- stats::approx(densite.G$x, densite.G$y, xout = y)$y
      } else {
        densite.G.dataPoint.y <- 1 / length(table(c(sample1,sample2)))
      }
      psi1.2(y) * estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2)) * densite.G.dataPoint.y
    }

    sigma <- matrix(NA, nrow = 3, ncol = 3)
    if (support == "continuous") {
      sigma[1,1] <- pracma::integral2(Vectorize(integrand.sigma2.11_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,2] <- sigma[2,1] <- pracma::integral2(Vectorize(integrand.sigma2.12_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[1,3] <- stats::integrate(Vectorize(integrand.sigma2.13_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[2,2] <- pracma::integral2(Vectorize(integrand.sigma2.22_g1diffg2), xmin=supp.integration[1], xmax=supp.integration[2], ymin=supp.integration[1], ymax=supp.integration[2])$Q
      sigma[2,3] <- stats::integrate(Vectorize(integrand.sigma2.23_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,1] <- stats::integrate(Vectorize(integrand.sigma2.31_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,2] <- stats::integrate(Vectorize(integrand.sigma2.32_g1diffg2), lower=supp.integration[1], upper=supp.integration[2], subdivisions = 10000L, rel.tol = 1e-3)$value
      sigma[3,3] <- estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2))
    } else {
      sigma[1,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma2.11_g1diffg2, x, supp.integration)) )
      sigma[1,2] <- sigma[2,1] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma2.12_g1diffg2, x, supp.integration)) )
      sigma[1,3] <- sum( unlist(sapply(supp.integration, integrand.sigma2.13_g1diffg2)) )
      sigma[2,2] <- sum( sapply(supp.integration, function(x) mapply(integrand.sigma2.22_g1diffg2, x, supp.integration)) )
      sigma[2,3] <- sum( unlist(sapply(supp.integration, integrand.sigma2.23_g1diffg2)) )
      sigma[3,1] <- sum( unlist(sapply(supp.integration, integrand.sigma2.31_g1diffg2)) )
      sigma[3,2] <- sum( unlist(sapply(supp.integration, integrand.sigma2.32_g1diffg2)) )
      sigma[3,3] <- estimVarCov_empProcess(x,y,sample2,known.p[2],list(comp.dist$f2,comp.dist$g2),list(comp.param$f2,comp.param$g2))
    }
  }

  return(sigma)
}

