#' Kernel estimation
#'
#' Functions to perform the estimation of cumulative distribution function (cdf) by kernel estimators
#' (with a non-gaussian kernel).
#'
#' @param h window of the kernel estimation.
#' @param u the point at which the estimation is made.
#'
#' @return the estimated value of the cdf.
#'
#' @examples
#' kernel_cdf(0.4,0.5)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

kernel_cdf <- function(u, h)
{
  return( (0.5 + u/h + (u/h)^2 / 2) * ((u > -h) & (u <= 0)) +
            (0.5 + u/h - (u/h)^2 / 2) * ((0 < u) & (u < h)) + 1*(u >= h) )
}
