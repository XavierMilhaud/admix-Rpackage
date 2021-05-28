#' Compute expansion coefficients in a given orthonormal polynomial basis.
#'
#' Compute the coefficients corresponding to the decomposition of some density in a given orthonormal polynomial basis.
#'
#' @param data Observed sample from which the coefficients are calculated. Can be NULL if 'comp.dist' and 'comp.param' are specified.
#' @param comp.dist (default to NULL) A list with two elements corresponding to component distributions (specified with
#'                  R native names for these distributions) involved in the admixture model.
#'                  Unknown elements must be specified as 'NULL' objects (for instance unknown 'f': list(f=NULL, g='norm')).
#' @param comp.param (default to NULL) A list with two elements corresponding to the parameters of the component distributions, each
#'                   element being a list itself. The names used in this list must correspond to the native R argument names
#'                   for these distributions. Unknown elements must be specified as 'NULL' objects.
#'                   For instance if 'f' is unknown: list(f = NULL, g = list(mean=0,sd=1)).
#' @param supp Support of the density considered.
#' @param degree Degree up to which the polynomial basis is built.
#' @param m (default to 3) Only used when support is 'Integer'. Corresponds to the mean of the reference measure, i.e. Poisson(m).
#' @param other (default to NULL) A list to precise bounds when the support is bounded, where the second and fourth elements give bounds.
#'
#' @return The list composed of 'degree' elements, each element being a numeric vector (with sample size) where each value represents
#'         the k-th order coefficient found when decomposing the density in the orthonormal polynomial basis.
#'
#' @examples
#' ## Simulate data:
#' sample1 <- rnorm(n = 7000, mean = 3, sd = 1)
#' ## Compute the expansion coefficients in the orthonormal polynomial basis:
#' coeff <- orthoBasis_coef(data = sample1, comp.dist = NULL, comp.param = NULL, supp = 'Real',
#'                          degree = 3, m = 3, other = NULL)
#' sapply(coeff, mean)
#' ## No observed data and decomposition of the known component of the admixture model:
#' coeff <- orthoBasis_coef(data = NULL, comp.dist = list(NULL, 'norm'),
#'             comp.param=list(NULL,list(mean=3,sd=1)), supp = 'Real', degree=3, m=3, other = NULL)
#' sapply(coeff, mean)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

orthoBasis_coef <- function(data, comp.dist = NULL, comp.param = NULL, supp = c('Real','Integer','Positive','Bounded.continuous'),
                            degree, m = 3, other = NULL)
{
  if (is.null(data)) {
    stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
    if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component of the admixture model must be specified.")
    ## Extracts the information on component distributions and stores in expressions:
    comp.dist.sim <- paste0("r", comp.dist[[2]])
    comp_ortho <- sapply(X = comp.dist.sim, FUN = get, pos = "package:stats", mode = "function")
    assign(x = names(comp_ortho)[1], value = comp_ortho[[1]])
    expr.sim <- paste(names(comp_ortho)[1],"(n=1000000,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
    data <- eval(parse(text = expr.sim))
  }

  ## Builds the orthonormal polynomial basis:
  poly_basis <- poly_orthonormal_basis(support = supp, deg = degree, x = data, m = m)

  ## Store the coefficients to be computed :
  coef.list <- vector(mode = "list", length = degree)

  ## Depending on the support :
  if (supp == "Real") {
    ## Reference measure N(0,1)
    for (i in 1:degree) coef.list[[i]] <- orthopolynom::polynomial.values(poly_basis, data)[[i+1]] / sqrt(factorial(i))

  } else if (supp == "Integer") {
    ## Reference measure P(3)
    for (i in 1:degree) coef.list[[i]] <- poly_basis[ ,i]

  } else if (supp == "Positive") {
    ## Reference measure Exp(1)
    for (i in 1:degree) coef.list[[i]] <- orthopolynom::polynomial.values(poly_basis, data)[[i+1]]

  } else if (supp == "Bounded.continuous") {
    ## Reference measure Unif(a,b)
    if (is.null(other)) { bounds <- c(min(data), max(data))
    } else { bounds <- other[[2]] }
    for (i in 1:degree) coef.list[[i]] <- (orthopolynom::polynomial.values(poly_basis, (2*data-bounds[1]-bounds[2])/(bounds[2]-bounds[1]))[[i+1]]) / sqrt(2*i+1)

  } else stop("Change the support since the choosen one is not considered!")

  return(coef.poly = coef.list)
}
