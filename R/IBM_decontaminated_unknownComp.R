#' Provide the decontaminated density of the unknown component in an admixture model
#'
#' Estimate (and plot) the decontaminated density of the unknown component in the admixture models to compare, thanks to the Inversion
#' step related to the Inversion - Best Matching (IBM) method. See 'Details' for further information on this estimation technique.
#'
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
#' @param estim.obj an R object obtained from the estimation of the component weights related to the proportions of the unknown component
#'                  in each of the two admixture models studied.
#' @param add_plot a boolean (TRUE by default) specifying if one plots the decontaminated densities of the two admixture models,
#'                 for visual comparison purpose.
#'
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03201760.
#'          Still not implemented for admixture of multinomial distributions.
#'
#' @return A list containing two elements: 1) the decontaminated density f1 of the 1st admixture model, 2) the same for the 2nd admixture model.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the mixture weight in each of the sample in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
#'                           comp.param = list.param, with.correction = FALSE, n.integ = 1000)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' res <- IBM_decontaminated_unknownComp(sample1 = sample1[['mixt.data']],
#'                                       sample2 = sample2[['mixt.data']],
#'                                       comp.dist = list.comp, comp.param = list.param,
#'                                       estim.obj = estimate, add_plot = TRUE)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_decontaminated_unknownComp <- function(sample1, sample2, comp.dist, comp.param, estim.obj, add_plot = TRUE)
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
  exp.comp.dist <- paste0("d", comp.dist)
  if (any(exp.comp.dist == "dmultinom")) { exp.comp.dist[which(exp.comp.dist == "dmultinom")] <- "stepfun" }
  comp_decontamin <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  for (i in 1:length(comp_decontamin)) assign(x = names(comp_decontamin)[i], value = comp_decontamin[[i]])

  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_decontamin)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                          paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_decontamin)[i],"(x,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Evaluate the density of known components:
  if (any(exp.comp.dist == "stepfun")) {
    g1.fun <- eval(parse(text = expr[2]))
    g2.fun <- eval(parse(text = expr[4]))
    g1 <- function(x) g1.fun(x)
    g2 <- function(x) g2.fun(x)
  } else {
    g1 <- function(x) { eval(parse(text = expr[2])) }
    g2 <- function(x) { eval(parse(text = expr[4])) }
  }
  ## Retrieves the empirical cumulative distribution functions:
  l1_emp <- stats::approxfun(stats::density(sample1))
  l2_emp <- stats::approxfun(stats::density(sample2))

  ## Range of observations:
  x_values <- seq(from = min(sample1, sample2), to = max(sample1,sample2), length.out = 1000)

  if (length(estim.obj[["prop.estim"]]) == 2) {
    ## Case when G1 != G2
    f1_decontamin <- f2_decontamin <- numeric(length = 1000)
    f1_decontamin <- sapply(X = x_values, FUN = function(x) (1/estim.obj[["prop.estim"]][1]) * (l1_emp(x) - (1-estim.obj[["prop.estim"]][1]) * g1(x)))
    f1_decontamin[which(is.na(f1_decontamin))] <- 0
    f2_decontamin <- sapply(X = x_values, FUN = function(x) (1/estim.obj[["prop.estim"]][2]) * (l2_emp(x) - (1-estim.obj[["prop.estim"]][2]) * g2(x)))
    f2_decontamin[which(is.na(f2_decontamin))] <- 0
  } else {
    ## Case when G1 = G2
    f1_decontamin <- f2_decontamin <- numeric(length = 1000)
    f1_decontamin <- sapply(X = x_values, FUN = function(x) (1/estim.obj[["p.X.fixed"]]) * (l1_emp(x) - (1-estim.obj[["p.X.fixed"]]) * g1(x)))
    f1_decontamin[which(is.na(f1_decontamin))] <- 0
    f2_decontamin <- sapply(X = x_values, FUN = function(x) (1/estim.obj[["prop.estim"]]) * (l2_emp(x) - (1-estim.obj[["prop.estim"]]) * g2(x)))
    f2_decontamin[which(is.na(f2_decontamin))] <- 0
  }

  if (add_plot) {
    x_range <- c(min(x_values), max(x_values))
    y_range <- c(min(f1_decontamin,f2_decontamin), max(f1_decontamin,f2_decontamin))
    plot(x = x_values, y = f1_decontamin, xlim = c(x_range[1],x_range[2]), ylim = c(y_range[1],y_range[2]), type = "l")
    old_par <- graphics::par()$new
    graphics::par(new = TRUE)
    plot(x = x_values, y = f2_decontamin, xlim = c(x_range[1],x_range[2]), ylim = c(y_range[1],y_range[2]), type = "l", col = "blue")
    on.exit(graphics::par(old_par))
  }

  return( list(density_decontamin_F1 = f1_decontamin, density_decontamin_F2 = f2_decontamin) )
}
