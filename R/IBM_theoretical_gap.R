#' Difference between unknown cumulative distribution functions of admixture models at some given point
#'
#' Compute the gap between the unknown cumulative distribution functions of the two considered admixture models at some given point,
#' where each admixture model has probability distribution function (pdf) given by l where l = p*f + (1-p)*g.
#' Uses the inversion method to do so, i.e. f = (1/p) (l - (1-p)g), where g represents the known component of the admixture
#' model and p is the proportion of the unknown component. This difference must be integrated over some domain to compute the
#' global discrepancy, as introduced in the paper presenting the IBM approach (see 'Details' below).
#'
#' @param z Point at which the difference between the unknown component distributions of the two considered admixture models is computed.
#' @param par Numeric vector with two elements, corresponding to the two parameter values at which to compute the gap. In practice
#'            the component weights for the two admixture models.
#' @param known.p Numeric vector with two elements, the known (true) mixture weights.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. No unknown elements permitted.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1='rnorm', g1='norm', f2='rnorm', g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. No unknown elements permitted. For instance, 'comp.param' could be specified
#'                   as follows: : list(f1 = list(mean=2,sd=0.3), g1 = list(mean=0,sd=1), f2 = list(mean=2,sd=0.3), g2 = list(mean=3,sd=1.1)).
#'
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03201760
#'
#' @return The gap between F1 and F2 (unknown components of the two admixture models), evaluated at the specified point.
#'
#' @examples
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' IBM_theoretical_gap(z = 2.8, par = c(0.3,0.6), known.p = c(0.5,0.5),
#'                     comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_theoretical_gap <- function(z, par, known.p = c(0.5,0.5), comp.dist, comp.param)
{
  if (is.null(known.p)) stop("In 'IBM_theoretical_gap': argument 'known.p' must be specified.")
  if ( (length(comp.dist) != 4) | (length(comp.param) != 4) | is.null(known.p) )
  { stop("The theoretical version of the gap can only be computed when all mixture components / weights are specified") }

  ## Extracts the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
  comp_gap <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
#  comp_gap <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  for (i in 1:length(comp_gap)) assign(x = names(comp_gap)[i], value = comp_gap[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_gap)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                           paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_gap)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  if (any(exp.comp.dist == "stepfun")) {
    F1.fun <- eval(parse(text=expr[1]))
    G1.fun <- eval(parse(text=expr[2]))
    F2.fun <- eval(parse(text=expr[3]))
    G2.fun <- eval(parse(text=expr[4]))
    F1 <- function(z) F1.fun(z)
    G1 <- function(z) G1.fun(z)
    F2 <- function(z) F2.fun(z)
    G2 <- function(z) G2.fun(z)
  } else {
    F1 <- function(z) { eval(parse(text=expr[1])) }
    G1 <- function(z) { eval(parse(text=expr[2])) }
    F2 <- function(z) { eval(parse(text=expr[3])) }
    G2 <- function(z) { eval(parse(text=expr[4])) }
  }

  L1.theo <- function(z) { known.p[1] * F1(z) + (1-known.p[1]) * G1(z) }
  L2.theo <- function(z) { known.p[2] * F2(z) + (1-known.p[2]) * G2(z) }

  ## Calcul explicite de la difference:
  if (length(par) == 2) {
    F.X.dataPoint <- (1/par[1]) * (L1.theo(z) - (1-par[1]) * G1(z))
    F.Y.dataPoint <- (1/par[2]) * (L2.theo(z) - (1-par[2]) * G2(z))
  } else {
    F.X.dataPoint <- (1/known.p[1]) * (L1.theo(z) - (1-known.p[1]) * G1(z))
    F.Y.dataPoint <- (1/par) * (L2.theo(z) - (1-par) * G2(z))
  }
  difference <- F.X.dataPoint - F.Y.dataPoint

  return(difference)
}
