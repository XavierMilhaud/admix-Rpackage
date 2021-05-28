#' Simulation of a Gaussian process
#'
#' Simulate the trajectory of a Gaussian process, given a mean vector and a variance-covariance structure.
#'
#' @param mean_vec Vector (if multimensional) of means for the increments following gaussian distribution.
#' @param varCov_mat Corresponding variance-covariance structure.
#' @param from Initial time point at which the process is simulated.
#' @param to Last time point at which the process is simulated.
#' @param start Useful if the user wants to make the trajectory start from some given value.
#' @param nb.points Number of points at which the process is simulated.
#'
#' @return The trajectory of the Gaussian processes after simulating the multivariate Gaussian distributions with
#'         specified variance-covariance structure.
#'
#' @examples
#' list.comp <- list(f1 = "norm", g1 = "norm")
#' list.param <- list(f1 = list(mean = 12, sd = 0.4),
#'                    g1 = list(mean = 16, sd = 0.7))
#' sample1 <- rsimmix(n = 2000, unknownComp_weight = 0.5, comp.dist = list.comp,
#'                    comp.param = list.param)$mixt.data
#' ## First get the variance-covariance matrix of the empirical process (Donsker correlation):
#' cov_mat <- .Call('_admix_estimVarCov_empProcess_Rcpp', PACKAGE = 'admix',
#'                   seq(from = min(sample1), to = max(sample1), length.out = 100), sample1)
#' ## Plug it into the simulation of the gaussian process:
#' B1 <- sim_gaussianProcess(mean_vec=rep(0,nrow(cov_mat)), varCov_mat=cov_mat, from=min(sample1),
#'                           to = max(sample1), start = 0, nb.points = nrow(cov_mat))
#' plot(x = B1$dates, y = B1$traj1, type="l", xlim = c(min(sample1),max(sample1)), ylim = c(-1,1))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

sim_gaussianProcess <- function(mean_vec, varCov_mat, from = 0, to = 1, start = 0, nb.points = 10)
{
  t <- seq(from = from, to = to, length.out = nb.points)
  rownames(varCov_mat) <- paste("x=", round(t,2), sep = "")
  colnames(varCov_mat) <- paste("y=", round(t,2), sep = "")

  ## Simulation of the multivariate gaussian distributions, without trend :
  trajectory.L1 <- MASS::mvrnorm(n = 1, m = mean_vec, Sigma = varCov_mat)
  traj1 <- trajectory.L1 - trajectory.L1[1] + start

  return( list(dates = t, traj1 = traj1) )
}
