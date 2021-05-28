#' Build an orthonormal basis to decompose some given probability density function
#'
#' Build an orthonormal basis, needed to decompose the probability density function (pdf) of the unknown component
#' from the admixture, depending on the support under consideration.
#'
#' @param support Support of the random variables implied in the two-component mixture distribution.
#' @param deg Degree up to which the basis is built.
#' @param x (NULL by default) Only used when support is 'Integer'. The point at which the polynomial value will be evaluated.
#' @param m (NULL by default) Only used when support is 'Integer'. Corresponds to the mean of the reference measure, i.e. Poisson(m).
#'
#' @return the orthonormal polynomial basis used to decompose the density of the unknown component of the mixture distribution.
#'
#' @examples
#' poly_orthonormal_basis(support = 'Real', deg = 10, x = NULL, m = NULL)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

poly_orthonormal_basis <- function(support = c("Real","Integer","Positive","Bounded.continuous","Bounded.discrete"), deg, x, m)
{
  if (support == "Real") {
    ## Hermite polynomials :
    poly_basis <- orthopolynom::hermite.he.polynomials(n = deg, normalized = FALSE)

  } else if (support == "Integer") {
    ## Charlier polynomials : 'm' is the mean of the reference measure, i.e. P(m).
    poly_basis <- function(x, m, deg) {
      CH <- CHA <- matrix(NA, nrow = length(x), ncol = deg)
      CH[ ,1] <- (m-x) / m
      CH[ ,2] <- ((m+1-x) * (m-x) / m - 1) / m
      for (k in 3:deg) { CH[ ,k] <- ((m+k-1-x) * CH[ ,k-1] - (k-1) * CH[ ,k-2]) / m }
      for (k in 1:deg) { CHA[ ,k] <- CH[ ,k] / (sqrt(factorial(k)*m^{-k})) }
      return(CHA)
    }

  } else if (support == "Positive") {
    ## Laguerre polynomials :
    poly_basis <- orthopolynom::laguerre.polynomials(n = deg, normalized = FALSE)

  } else if (support == "Bounded.continuous") {
    ## Legendre polynomials :
    poly_basis <- orthopolynom::legendre.polynomials(n = deg, normalized = FALSE)

  } else stop("Please give a correct argument for the support of the distributions!")

  return(poly_basis)
}
