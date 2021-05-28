#' Empirical computation of the contrast in the Inversion - Best Matching (IBM) method
#'
#' Defines the empirical version of the contrast in the IBM method, to be minimized in the optimization process. For further details
#' about the contrast definition, see 'Details' below.
#'
#' @param par Numeric vector with two elements, corresponding to the two parameter values at which to compute the contrast. In practice
#'            the component weights for the two admixture models.
#' @param fixed.p.X Arbitrary value chosen by the user for the component weight related to the unknown component of the first
#'                  admixture model. Only useful for optimization when the known components of the two models are identical
#'                  (G1=G2, leading to unidimensional optimization).
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
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03201760
#'
#' @return The empirical contrast value evaluated at parameter values.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                comp.param = list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                comp.param = list(list.param$f2,list.param$g2))
#' ## Create the distribution on which the contrast will be integrated:
#' G <- stats::rnorm(n = 1000, mean = sample(c(sample1[['mixt.data']], sample2[['mixt.data']]),
#'                                           size = 1000, replace = TRUE),
#'                   sd = density(c(sample1[['mixt.data']], sample2[['mixt.data']]))$bw)
#' ## Compute the empirical contrast at parameters (p1,p2) = (0.2,0.7) in a real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' IBM_empirical_contrast(par = c(0.2,0.7), fixed.p.X = NULL, sample1 = sample1[['mixt.data']],
#'            sample2= sample2[['mixt.data']], G=G, comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_empirical_contrast <- function(par, fixed.p.X = NULL, sample1, sample2, G, comp.dist, comp.param)
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

  ## Extract the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
  comp_emp <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_emp)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                      paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_emp)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
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

  ## Empirical cumulative distribution function from the two observed samples:
  L1 <- stats::ecdf(sample1)						# hat(L1)
  L2 <- stats::ecdf(sample2)						# hat(L2)

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)

  ##------- Computes the empirical contrast --------##
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    densite.G <- stats::density(G, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integration <- c(min(G), max(G))
  } else {
    supp.integration <- G
  }

  if (is.null(fixed.p.X)) {
      integrand <- function(z, par) {
        if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
        } else {
          densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
        }
        F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
        F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
        weighted.difference <- (F.X.dataPoint - F.Y.dataPoint)^2 * densite.G.dataPoint
        weighted.difference
      }
  } else {   # for one-dimensional optimization
    if (G1equalG2) stopifnot(!is.null(fixed.p.X))
      integrand <- function(z, par) {
        if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
        } else {
          densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
        }
        F.X.dataPoint <- (1/fixed.p.X) * (L1(z) - (1-fixed.p.X) * G1(z))
        F.Y.dataPoint <- (1/par) * (L2(z) - (1-par) * G2(z))
        weighted.difference <- (F.X.dataPoint - F.Y.dataPoint)^2 * densite.G.dataPoint
        weighted.difference
      }
  }

  if (support == "continuous") {
    res <- stats::integrate(integrand, lower = supp.integration[1], upper = supp.integration[2], par, subdivisions = 10000L,
                            rel.tol = 1e-03)$value
  } else {
    res <- sum( unlist(sapply(supp.integration, integrand, par)) )
  }

  return(res)
}
