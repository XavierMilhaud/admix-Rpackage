#' Kernel estimation
#'
#' Functions to perform the estimation of probability density function (pdf) by kernel estimators (with a non-gaussian kernel).
#'
#' @param h window of the kernel estimation.
#' @param u the point at which the estimation is made.
#'
#' @return the estimated value of the pdf.
#'
#' @examples
#' kernel_density(0.4,0.5)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

kernel_density <- function(u, h)
{
  return( (1 - abs(u/h)) * ((u > -h) & (u < h)) / h )
}
