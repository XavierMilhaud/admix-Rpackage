#' Difference between the unknown empirical cumulative distribution functions in two admixture models
#'
#' Compute the 'gap' between two unknown cumulative distribution functions (ecdf) at some given point, in admixture models
#' with probability distribution function (pdf) given by l where l = p*f + (1-p)*g.
#' Uses the inversion method to do so, i.e. f = (1/p) (l - (1-p)*g), where g represents the known component of the admixture
#' model and p is the unknown proportion of the unknown component. Therefore, compute:
#'    D(z,L1,L2,p1,p2) = F1(z,L1,p1) - F2(z,L2,p2)
#' This measure should be integrated over some domain to compute the global discrepancy, see further information in 'Details' below.
#'
#' @param z the point at which the difference between both unknown (estimated) component distributions is computed.
#' @param par Numeric vector with two elements, corresponding to the weights of the unknown component for the two admixture models.
#' @param fixed.p1 (optional, NULL by default) Arbitrary value chosen by the user for the component weight related to the unknown
#'                  component of the first admixture model. Only useful for optimization when the known components of the two
#'                  models are identical (G1=G2, leading to unidimensional optimization).
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
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
#' @return the gap evaluated at the specified point between the unknown components of the two observed samples.
#'
#' @examples
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' IBM_gap(z = 2.8, par = c(0.3,0.6), fixed.p1 = NULL, sample1 = sample1[['mixt.data']],
#'         sample2 = sample2[['mixt.data']], comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_gap <- function(z, par, fixed.p1 = NULL, sample1, sample2, comp.dist, comp.param)
{
  if (is.null(fixed.p1)) stopifnot(length(par) == 2)
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
#  comp_IBM_gap <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_IBM_gap <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_IBM_gap)) assign(x = names(comp_IBM_gap)[i], value = comp_IBM_gap[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_IBM_gap)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                          paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_IBM_gap)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  if (any(exp.comp.dist == "stepfun")) {
    G1.fun <- eval(parse(text = expr[2]))
    G2.fun <- eval(parse(text = expr[4]))
    G1 <- function(z) G1.fun(z)
    G2 <- function(z) G2.fun(z)
  } else {
    G1 <- function(z) { eval(parse(text = expr[2])) }
    G2 <- function(z) { eval(parse(text = expr[4])) }
  }
  ## Calcul des fonctions de repartition empiriques utiles provenant des echantillons X et Y:
  L1 <- stats::ecdf(sample1)								    # L1
  L2 <- stats::ecdf(sample2)								    # L2

  ## Calcul explicite de la difference:
  if (is.null(fixed.p1)) {
    F.X.dataPoint <- (1/par[1]) * (L1(z) - (1-par[1]) * G1(z))
    F.Y.dataPoint <- (1/par[2]) * (L2(z) - (1-par[2]) * G2(z))
  } else {
    F.X.dataPoint <- (1/fixed.p1) * (L1(z) - (1-fixed.p1) * G1(z))
    F.Y.dataPoint <- (1/par) * (L2(z) - (1-par) * G2(z))
  }
  difference <- F.X.dataPoint - F.Y.dataPoint

  return(difference)
}

