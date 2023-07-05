#' Provide the decontaminated cumulative distribution function (CDF) of the unknown component in an admixture model
#'
#' Estimate the decontaminated CDF of the unknown component in the admixture model under study, after inversion of the admixture
#' cumulative distribution function. Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with g a known CDF and p and f unknown quantities.
#'
#' @param sample1 Observations of the sample under study.
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
#' @details The decontaminated CDF is obtained by inverting the admixture CDF, given by L = p*F + (1-p)*G, to isolate the
#'          unknown component F after having estimated p. This means that F = (1/hat(p)) * (hat(L)-(1-p)*G).
#'
#' @return The decontaminated CDF F of the admixture model, as an of class 'stepfun' (step function).
#'
#' @examples
#' ####### Continuous support:
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=3500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=3000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' ## Determine the decontaminated version of the unknown CDF by inversion:
#' decontaminated_cdf(sample1 = sample1[['mixt.data']], comp.dist = list.comp[1:2],
#'                     comp.param = list.param[1:2], estim.p = estimate$prop.estim[1])
#' ####### Countable discrete support:
#' list.comp <- list(f1 = 'pois', g1 = 'pois',
#'                   f2 = 'pois', g2 = 'pois')
#' list.param <- list(f1 = list(lambda = 3), g1 = list(lambda = 2),
#'                    f2 = list(lambda = 3), g2 = list(lambda = 4))
#' sample1 <- rsimmix(n=6000, unknownComp_weight=0.6, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=4500, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'pois',
#'                   f2 = NULL, g2 = 'pois')
#' list.param <- list(f1 = NULL, g1 = list(lambda = 2),
#'                    f2 = NULL, g2 = list(lambda = 4))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' decontaminated_cdf(sample1 = sample1[['mixt.data']], comp.dist = list.comp[1:2],
#'                    comp.param = list.param[1:2], estim.p = estimate$prop.estim[1])
#' ####### Finite discrete support:
#' list.comp <- list(f1 = 'multinom', g1 = 'multinom',
#'                   f2 = 'multinom', g2 = 'multinom')
#' list.param <- list(f1 = list(size=1, prob=c(0.3,0.4,0.3)), g1 = list(size=1, prob=c(0.6,0.3,0.1)),
#'                    f2 = list(size=1, prob=c(0.3,0.4,0.3)), g2 = list(size=1, prob=c(0.2,0.6,0.2)))
#' sample1 <- rsimmix(n=8000, unknownComp_weight=0.6, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=6000, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' list.comp <- list(f1 = NULL, g1 = 'multinom',
#'                   f2 = NULL, g2 = 'multinom')
#' list.param <- list(f1 = NULL, g1 = list(size=1, prob=c(0.6,0.3,0.1)),
#'                    f2 = NULL, g2 = list(size=1, prob=c(0.2,0.6,0.2)))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' decontaminated_cdf(sample1 = sample1[['mixt.data']], comp.dist = list.comp[1:2],
#'                     comp.param = list.param[1:2], estim.p = estimate$prop.estim[1])
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

decontaminated_cdf <- function(sample1, comp.dist, comp.param, estim.p)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) {
    stop("Known components must be specified in the admixture model.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1)]))) stop("Mixture distributions/parameters were not correctly specified")
  }

  ## Extract the information on component distributions:
  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
#  comp_decontamin <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_decontamin <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_decontamin)) assign(x = names(comp_decontamin)[i], value = comp_decontamin[[i]])

  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_decontamin)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                      paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_decontamin)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Differentiates the case of the multinomial distribution from other cases:
  if (any(exp.comp.dist == "stepfun")) {
    G1.fun <- eval(parse(text = expr[2]))
    G1 <- function(z) G1.fun(z)
  } else {
    G1 <- function(z) { eval(parse(text = expr[2])) }
  }
  ## Empirical cumulative distribution function from the observed sample:
  L1 <- stats::ecdf(sample1)						# hat(L1)
  #class(L1) ; plot.ecdf(L1)

  ##------- Defines the support --------##
  support <- detect_support_type(sample1)
  ## Range of observations:
  if (support == "continuous") {
    x_values <- seq(from = min(sample1), to = max(sample1), length.out = 3000)
  } else {
    x_values <- unique(sort(c(unique(sample1))))
  }

  F1_decontamin <- sapply(X = x_values, FUN = function(x) (1/estim.p) * (L1(x) - (1-estim.p) * G1(x)))
  F1_decontamin[which(is.na(F1_decontamin))] <- 0
  F1_decontamin <- pmax(F1_decontamin, 0)
  F1_decontamin <- pmin(F1_decontamin, 1)

  F1_fun <- stats::stepfun(x = x_values, y = c(0,F1_decontamin))
  #F1_fun <- Iso::pava(y = c(0,F1_decontamin), decreasing = FALSE, stepfun = TRUE)
  #plot.stepfun(x = F1_fun, xlim = c(min(x_values[knots(F1_fun)-1]), max(x_values[knots(F1_fun)-1])))
  return(F1_fun)
}
