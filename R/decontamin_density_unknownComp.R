#' Provide the decontaminated density of the unknown component in an admixture model
#'
#' Estimate the decontaminated density of the unknown component in the admixture model under study, after inversion of the admixture
#' cumulative distribution function. Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with g a known CDF and p and f unknown quantities.
#'
#' @param sample1 Observations of the first sample under study.
#' @param comp.dist A list with two elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the admixture model. The two elements refer to the unknown and known components of the admixture model,
#'                  If there are unknown elements, they must be specified as 'NULL' objects (e.g. 'comp.dist' could be set to list(f1=NULL, g1='norm')).
#' @param comp.param A list with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two elements refer to the parameters of unknown and known components of the admixture model. If there are unknown
#'                   elements, they must be specified as 'NULL' objects (e.g. 'comp.param' could be set to list(f1=NULL, g1=list(mean=0,sd=1))).
#' @param estim.p The estimated weight of the unknown component distribution, related to the proportion of the unknown component
#'                  in the admixture model studied.
#'
#' @details The decontaminated density is obtained by inverting the admixture density, given by l = p*f + (1-p)*g, to isolate the
#'          unknown component f after having estimated p.
#'
#' @return A list containing the decontaminated density of the admixture model (of class 'function'),
#'         and the support of the observations (either discrete or continuous).
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=8000, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=7000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' decontamin_density_unknownComp(sample1 = sample1[['mixt.data']],
#'                                comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                                estim.p = estimate$prop.estim[1])
#' ####### Discrete support:
#' list.comp <- list(f1 = 'pois', g1 = 'pois',
#'                   f2 = 'pois', g2 = 'pois')
#' list.param <- list(f1 = list(lambda = 3), g1 = list(lambda = 2),
#'                    f2 = list(lambda = 3), g2 = list(lambda = 4))
#' sample1 <- rsimmix(n=7000, unknownComp_weight=0.6, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=6000, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'pois',
#'                   f2 = NULL, g2 = 'pois')
#' list.param <- list(f1 = NULL, g1 = list(lambda = 2),
#'                    f2 = NULL, g2 = list(lambda = 4))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' decontamin_density_unknownComp(sample1 = sample1[['mixt.data']],
#'                                comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                                estim.p = estimate$prop.estim[1])
#' ####### Finite discrete support:
#' list.comp <- list(f1 = 'multinom', g1 = 'multinom',
#'                   f2 = 'multinom', g2 = 'multinom')
#' list.param <- list(f1 = list(size=1, prob=c(0.3,0.4,0.3)), g1 = list(size=1, prob=c(0.6,0.3,0.1)),
#'                    f2 = list(size=1, prob=c(0.3,0.4,0.3)), g2 = list(size=1, prob=c(0.2,0.6,0.2)))
#' sample1 <- rsimmix(n=12000, unknownComp_weight=0.6, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=10000, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'multinom',
#'                   f2 = NULL, g2 = 'multinom')
#' list.param <- list(f1 = NULL, g1 = list(size=1, prob=c(0.6,0.3,0.1)),
#'                    f2 = NULL, g2 = list(size=1, prob=c(0.2,0.6,0.2)))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' decontamin_density_unknownComp(sample1 = sample1[['mixt.data']],
#'                                comp.dist = list.comp[1:2], comp.param = list.param[1:2],
#'                                estim.p = estimate$prop.estim[1])
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

decontamin_density_unknownComp <- function(sample1, comp.dist, comp.param, estim.p)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) {
    stop("The known component must be specified in the admixture model under study.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1)]))) stop("Mixture distributions/parameters were not correctly specified")
  }

  ## Extract the information on component distributions:
  exp.comp.dist <- paste0("d", comp.dist)
  if (any(exp.comp.dist == "dmultinom")) { exp.comp.dist[which(exp.comp.dist == "dmultinom")] <- "approxfun" }
  comp_decontamin <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  for (i in 1:length(comp_decontamin)) assign(x = names(comp_decontamin)[i], value = comp_decontamin[[i]])

  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_decontamin)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = c(",
                                      paste(paste(comp.param[[i]][[2]], collapse = ","), ")", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_decontamin)[i],"(x,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "approxfun")] <- sapply(which(exp.comp.dist == "approxfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Evaluate density of known components, differentiates the case of the multinomial distribution from others:
  if (any(exp.comp.dist == "approxfun")) {
    g1.fun <- eval(parse(text = expr[2]))
    g1 <- function(x) g1.fun(x)
  } else {
    g1 <- function(x) { eval(parse(text = expr[2])) }
  }

  ##------- Defines the support --------##
  support <- detect_support_type(sample1)

  ## Retrieves the empirical cumulative distribution functions:
  if (support == "continuous") {
    l1_dens <- stats::density(sample1)
    l1_emp <- stats::approxfun(x = l1_dens$x, y = l1_dens$y)
  } else {
    l1_emp <- stats::approxfun(x = as.numeric(names(table(sample1))),
                               y = table(sample1) / sum(table(sample1)), method = "constant")
  }

  f1_decontamin <- function(x) (1/estim.p) * (l1_emp(x) - (1-estim.p) * g1(x))

  return( list(decontamin_f = f1_decontamin, supp = support) )
}
