#' Variance-covariance matrix of the empirical process in an admixture model
#'
#' Estimate the variance-covariance matrix of some given empirical process, based on the Donsker correlation.
#' Compute Donsker correlation between two time points (x,y) for some given empirical process with R code
#' (another implementation in C++ is also available to speed up this computation).
#'
#' @param x First time point considered for the computation of the correlation given the empirical process.
#' @param y Second time point considered for the computation of the correlation given the same empirical process.
#' @param obs.data Sample that permits to estimate the cumulative distribution function (cdf).
#' @param known.p NULL by default (only useful to compute the exact Donsker correlation). The component weight dedicated to
#'                the unknown mixture component if available (in case of simulation studies)
#' @param comp.dist NULL by default (only useful to compute the exact Donsker correlation). Otherwise, a list with two elements
#'                  corresponding to component distributions (specified with R native names for these distributions) involved
#'                  in the admixture model. All elements must be specified, for instance list(f='norm', g='norm').
#' @param comp.param NULL by default (only useful to compute the exact Donsker correlation). Otherwise, a list with two elements
#'                   corresponding to the parameters of the component distributions, each element being a list itself. The names
#'                   used in this list must correspond to the native R argument names for these distributions.
#'                   All elements must be specified, for instance list(f=NULL, g=list(mean=0,sd=1)).
#'
#' @return The estimated variance-covariance matrix.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm')
#' list.param <- list(f1 = list(mean = 12, sd = 0.4),
#'                    g1 = list(mean = 16, sd = 0.7))
#' obs.data <- rsimmix(n=2500, unknownComp_weight=0.5, comp.dist=list.comp, comp.param= list.param)
#' ## Compute the variance-covariance matrix of the corresponding empirical process:
#' t <- seq(from = min(obs.data$mixt.data), to = max(obs.data$mixt.data), length = 50)
#' S2 <- sapply(t, function(s1) {
#'                 sapply(t, function(s2) {
#'                      estimVarCov_empProcess(x = s1, y = s2, obs.data = obs.data$mixt.data) })
#'                 })
#' lattice::wireframe(S2)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

estimVarCov_empProcess <- function(x, y, obs.data, known.p = NULL, comp.dist = NULL, comp.param = NULL)
{
  if (is.null(obs.data)) {
    stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
    if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) | is.null(known.p)) {
      stop("All parameters of the admixture model must be specified to compute the exact Donsker correlation.")
    }
    ## Extracts the information on component distributions:
    exp.comp.dist <- paste0("p", comp.dist)
    if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
    comp.ro <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
    for (i in 1:length(comp.ro)) assign(x = names(comp.ro)[i], value = comp.ro[[i]])
    ## Creates the expression involved in future assessments of the CDF:
    make.expr.step <- function(i) paste(names(comp.ro)[i], "(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                  paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
    make.expr <- function(i) paste(names(comp.ro)[i], "(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
    expr <- vector(mode = "character", length = length(exp.comp.dist))
    expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
    expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
    expr <- unlist(expr)

    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      F1 <- function(z) F1.fun(z)
      G1.fun <- eval(parse(text = expr[2]))
      G1 <- function(z) G1.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
    }

    L.CDF <- function(z) { known.p * F1(z) + (1-known.p) * G1(z) }

  } else {

    L.CDF <- stats::ecdf(obs.data)
  }

  ## Use the Donsker correlation formula:
  res <- L.CDF(min(x,y)) * (1 - L.CDF(max(x,y)))

  return(res)
}
