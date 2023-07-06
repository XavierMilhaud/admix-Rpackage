#' Hessian matrix of the contrast function in the Inversion - Best Matching (IBM) method
#'
#' Compute the hessian matrix of the contrast as defined in the IBM approach, at point (p1,p2). Here, based on
#' two samples following admixture models, where we recall that admixture models have probability distribution function
#' (pdf) given by l where l = p*f + (1-p)*g, where g represents the only known quantity and l is the pdf of the observed sample.
#' See 'Details' below for further information about the definition of the contrast.
#'
#' @param par Numeric vector with two elements (corresponding to the two unknown component weights) at which the hessian is computed.
#' @param fixed.p1 (optional, NULL by default) Arbitrary value chosen by the user for the component weight related to the unknown
#'                  component of the first admixture model. Only useful for optimization when the known components of the two
#'                  models are identical (G1=G2, leading to unidimensional optimization).
#' @param known.p (optional, NULL by default) Numeric vector with two elements, the known (true) mixture weights.
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param G Distribution on which to integrate when calculating the contrast.
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
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.science/hal-03201760
#'
#' @return the hessian matrix of the contrast.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Define the distribution over which to integrate:
#' fit.all <- stats::density(x = c(sample1[['mixt.data']],sample2[['mixt.data']]))
#' G <- stats::rnorm(n = 1000, mean = sample(c(sample1[['mixt.data']], sample2[['mixt.data']]),
#'                                           size = 1000, replace = TRUE), sd = fit.all$bw)
#' ## Evaluate the hessian matrix at point (p1,p2) = (0.3,0.6):
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' IBM_hessian_contrast(par = c(0.3,0.6), fixed.p1 = NULL, known.p = NULL,
#'                      sample1 = sample1[['mixt.data']],  sample2 = sample2[['mixt.data']], G = G,
#'                      comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_hessian_contrast <- function(par, fixed.p1 = NULL, known.p = NULL, sample1, sample2, G, comp.dist = NULL, comp.param = NULL)
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
#  comp_hessian <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_hessian <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_hessian)) assign(x = names(comp_hessian)[i], value = comp_hessian[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_hessian)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_hessian)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
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

  ## Calcul des fonctions de repartition (empiriques provenant de l'echantillon X, ou exactes) au point 'z':
  if ( !is.null(known.p) ) {
    L1 <- function(z) { known.p[1] * F1(z) + (1-known.p[1]) * G1(z) }
    L2 <- function(z) { known.p[2] * F2(z) + (1-known.p[2]) * G2(z) }
  } else {
    L1 <- stats::ecdf(sample1)  # hat(L1)
    L2 <- stats::ecdf(sample2)  # hat(L2)
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

  ## Differentiates cases where G1 = G2 and G1 != G2 :
  if (G1equalG2) {

    stopifnot(length(par) == 1)
    integrand <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/fixed.p1) * (L1(z) - (1-fixed.p1) * G1(z))
      F.Y.dataPoint <- (1/par) * (L2(z) - (1-par) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.dataPoint <- (1/par^2) * (L2(z) - G2(z))
      hessianD.dataPoint <- -(2/par^3) * (L2(z) - G2(z))
      (hessianD.dataPoint * D + gradD.dataPoint^2) * densite.G.dataPoint
    }

    if (support == "continuous") {
      err.tol <- 1e-03
      hessian.contr <- try(suppressWarnings(
        2 * stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L,
                             rel.tol = err.tol)$value), silent = TRUE)
      while ( (class(hessian.contr) == "character") & (err.tol < 1) ) {
        err.tol <- err.tol * 10
        hessian.contr <- try(suppressWarnings(
          2 * stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L,
                               rel.tol = err.tol)$value), silent = TRUE)
      }
    } else {
      hessian.contr <- 2 * sum( unlist(sapply(supp.integration, integrand, par)) )
    }

  } else {

    integrand.1.1 <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        #densite.G.dataPoint <- (table(c(sample1,sample2)) / sum(table(c(sample1,sample2))))[which(names(table(c(sample1,sample2))) == z)]
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
      F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.pX.dataPoint <- (1/par[1]^2) * (G1(z) - L1(z))
      hessianD.dataPoint <- -(2/par[1]^3) * (G1(z) - L1(z))
      (gradD.pX.dataPoint^2 + D * hessianD.dataPoint) * densite.G.dataPoint
    }
    integrand.1.2 <- integrand.2.1 <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      gradD.pX.dataPoint <- (1/par[1]^2) * (G1(z) - L1(z))
      gradD.pY.dataPoint <- -(1/par[2]^2) * (G2(z) - L2(z))
      (gradD.pX.dataPoint * gradD.pY.dataPoint) * densite.G.dataPoint
    }
    integrand.2.2 <- function(z, par) {
      if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
      } else {
        densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
      }
      F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
      F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
      D <- (F.X.dataPoint - F.Y.dataPoint)
      gradD.pY.dataPoint <- -(1/par[2]^2) * (G2(z) - L2(z))
      hessianD.dataPoint <- (2/par[1]^3) * (G2(z) - L2(z))
      (gradD.pY.dataPoint^2 + D * hessianD.dataPoint) * densite.G.dataPoint
    }

    hessian.contr <- matrix(NA, nrow = 2, ncol = 2)
    if (support == "continuous") {

      err.tol <- 1e-03
      term1 <- try(2 * stats::integrate(integrand.1.1, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term1) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term1 <- try(2 * stats::integrate(integrand.1.1, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      }
      hessian.contr[1,1] <- term1

      err.tol <- 1e-03
      term2 <- try(2 * stats::integrate(integrand.1.2, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term2) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term2 <- try(2 * stats::integrate(integrand.1.2, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      }
      hessian.contr[1,2] <- term2
      hessian.contr[2,1] <- hessian.contr[1,2]

      err.tol <- 1e-03
      term3 <- try(2 * stats::integrate(integrand.2.2, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
      while ((class(term3) == "try-error") & (err.tol < 1)) {
        err.tol <- err.tol * 10
        term3 <- try(2 * stats::integrate(integrand.2.2, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L, rel.tol = err.tol)$value, silent = TRUE)
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
