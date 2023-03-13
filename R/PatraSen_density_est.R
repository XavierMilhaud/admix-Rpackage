#' Compute the estimate of the density of the unknown component in an admixture model
#'
#' Compute by Patra and Sen technique the estimate of f.s (density corresponding to F.s) when f.s is known to
#' be either decreasing or increasing.
#'
#' @param input an R object of class 'cv.mixmodel' or 'mixmodel'.
#' @param dec.density a boolean indicating whether the density is increasing or decreasing.
#'
#' @details See Patra, R.K. and Sen, B. (2016); Estimation of a Two-component Mixture Model with Applications to Multiple Testing;
#'          JRSS Series B, 78, pp. 869--893.
#'
#' @return an estimator of the unknown component density.
#'
#' @examples
#' comp.dist <- list(f = 'norm', g = 'norm')
#' comp.param <- list(f = list(mean = 3, sd = 0.5), g = list(mean = 0, sd = 1))
#' data1 <- rsimmix(n = 2000, unknownComp_weight = 0.6, comp.dist, comp.param)[['mixt.data']]
#' data1_transfo <- knownComp_to_uniform(data = data1, comp.dist = list(comp.dist$f, comp.dist$g),
#'                                       comp.param = list(comp.param$f, comp.param$g))
#' res <- PatraSen_cv_mixmodel(data = data1_transfo, folds = 3, reps = 1, cn.s = NULL,
#'                             cn.length = 3, gridsize = 200)
#' PatraSen_density_est(res, dec.density = TRUE)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

PatraSen_density_est <- function(input, dec.density = TRUE)
{
	if (inherits(input, "cv.mixmodel") || inherits(input, "mixmodel")) {
		Fs.hat <- input$Fs.hat
	} else {
		stop("This function only works on objects of class 'cv.mixmodel' or 'mixmodel'. See functions 'est.mix.model' or 'cv.mixmodel'.")
	}
	if (dec.density == TRUE){
		ll <- fdrtool::gcmlcm(Fs.hat$x, Fs.hat$y, type = "lcm")
	} else if (dec.density == FALSE){
		ll <- fdrtool::gcmlcm(Fs.hat$x, Fs.hat$y, type = "gcm")
	}

	fs.hat <- NULL
	fs.hat$x <- rep(ll$x.knots, each = 2)                 # data points for density
	fs.hat$y <- c(0, rep(ll$slope.knots, each = 2), 0)    # value of density

	return(fs.hat)
}



