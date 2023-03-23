#' Estimate by Patra and Sen the unknown component weight as well as the unknown distribution in admixture models
#'
#' Estimation of unknown elements (by Patra and Sen method) under the admixture model with probability density function l:
#'          l = p*f + (1-p)*g,
#' where g is the known component of the two-component mixture, p is the unknown proportion of the unknown component distribution f.
#' More information in 'Details' below concerning the estimation method.
#'
#' @param data Sample where the known component density of the admixture model has been transformed into a Uniform(0,1) distribution.
#' @param method Either 'fixed' or 'cv', depending on whether compute the estimate based on the value of 'c.n' or
#'               use cross-validation for choosing 'c.n' (tuning parameter).
#' @param c.n A positive number, with default value equal to 0.1 log(log(n)), where 'n' is the length of the observed sample.
#' @param folds Number of folds used for cross-validation, default is 10.
#' @param reps Number of replications for cross-validation, default is 1.
#' @param cn.s A sequence of 'c.n' to be used for cross-validation (vector of values). Default is equally
#'            spaced grid of 100 values between .001 x log(log(n)) and 0.2 x log(log(n)).
#' @param cn.length (default to 100) Number of equally spaced tuning parameter (between .001 x log(log(n)) and 0.2 x log(log(n))).
#'                  Values to search from.
#' @param gridsize (default to 600) Number of equally spaced points (between 0 and 1) to evaluate the distance function.
#'                 Larger values are more computationally intensive but also lead to more accurate estimates.
#'
#' @details See Patra, R.K. and Sen, B. (2016); Estimation of a Two-component Mixture Model with Applications to Multiple Testing;
#'          JRSS Series B, 78, pp. 869--893.
#'
#' @return A list containing 'alp.hat' (estimate of the unknown component weight), 'Fs.hat' (list with elements 'x' and 'y' values for the function estimate
#'         of the unknown cumulative distribution function), 'dist.out' which is an object of the class 'dist.fun'
#'         using the complete data.gen, 'c.n' the value of the tuning parameter used to compute the final estimate,
#'         and finally 'cv.out' which is an object of class 'cv.mixmodel'. The object is NULL if method is "fixed".
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f = 'norm', g = 'norm')
#' list.param <- list(f = list(mean = 3, sd = 0.5),
#'                    g = list(mean = 0, sd = 1))
#' data1 <- rsimmix(n = 1500, unknownComp_weight = 0.8, list.comp, list.param)[['mixt.data']]
#' ## Transform the known component of the admixture model into a Uniform(O,1) distribution:
#' list.comp <- list(f = NULL, g = 'norm')
#' list.param <- list(f = NULL, g = list(mean = 0, sd = 1))
#' data1_transfo <- knownComp_to_uniform(data = data1, comp.dist=list.comp, comp.param=list.param)
#' PatraSen_est_mix_model(data = data1_transfo, method = 'fixed',
#'                        c.n = 0.1*log(log(length(data1_transfo))), gridsize = 1000)$alp.hat
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

PatraSen_est_mix_model <- function(data, method = c("lwr.bnd", "fixed", "cv"), c.n = NULL, folds = 10, reps = 1,
                                   cn.s = NULL, cn.length = 100, gridsize = 600)
{
	if (!is.vector(data)) stop(" 'data' has to be a vector.")
	n <- length(data)
	if (is.null(method)) stop("'method' can not be NULL")

	if (method == 'lwr.bnd') {
		dist.out <- PatraSen_dist_calc(data, gridsize = gridsize)
		q <- 0.6792
		alp.Lwr <- sum(dist.out$distance > q/ sqrt(n)) / gridsize
		alp.hat <- NULL
		Fs.hat.fun <- NULL
		c.n <- NULL
	} else {
		alp.Lwr <- NULL
	}

	if (method == "fixed"){
		if (is.null(c.n)) {
			warning("'c.n' is not given. Fixing it to be '0.1*log(log(n))")
			c.n <- 0.1 *log(log(n))
		}
		dist.out <- PatraSen_dist_calc(data, gridsize = gridsize)
		alp.hat <- sum(dist.out$distance > c.n/ sqrt(n)) / gridsize
		if (alp.hat > 0) {
			F.hat <- (dist.out$F.n - (1-alp.hat)*dist.out$F.b) / alp.hat    # computes the naive estimator of F_s
			Fs.hat <- Iso::pava(F.hat, dist.out$Freq, decreasing = FALSE)   # computes the Isotonic Estimator of F_s
			Fs.hat[which(Fs.hat <= 0)] <- 0
			Fs.hat[which(Fs.hat >= 1)] <- 1
		} else {
			Fs.hat <- dist.out$F.b
		}
		Fs.hat.fun <- NULL
		Fs.hat.fun$y <- Fs.hat
		Fs.hat.fun$x <-  dist.out$F.n.x

	} else if (method == "cv") {
		out.cv <- PatraSen_cv_mixmodel(data, folds = folds, reps = reps, cn.s = cn.s, cn.length = cn.length, gridsize = gridsize)
		alp.hat <- out.cv$alp.hat
		Fs.hat.fun <- out.cv$Fs.hat
		dist.out <- out.cv$dist.out
		c.n <- out.cv$cn.cv
	}

	ret <- list(alp.hat = alp.hat,
	            Fs.hat = Fs.hat.fun,
	            dist.out = dist.out,
	            c.n = c.n,
	            alp.Lwr =alp.Lwr,
	            n = n)

	if (method == "cv"){
		ret$cv.out <- out.cv
	} else {
		ret$cv.out <- NULL
	}
	ret$method = method
	ret$call <- match.call()

	class(ret) <- "mixmodel"
	return(ret)
}


print.mixmodel <- function(x, ...){
	cat("Call:")
	print(x$call)
	if(x$method != "lwr.bnd"){
		print(paste("Estimate of alp is" , x$alp.hat))
		print(paste(" The chosen value c_n is", x$c.n))
		if( !is.null(x$cv.out)){
		  old_par <- graphics::par()$mfrow
		  graphics::par(mfrow=c(1,2))
			plot(x$cv.out)
			on.exit(graphics::par(old_par))
		}
		plot(x$dist.out)
	} else if(x$method == 'lwr.bnd'){
		plot(x$dist.out)
		print (paste("The  '95%' lower confidence for alp_0 is ", x$alp.Lwr))
	}
}

plot.mixmodel <- function(x, ...){
	if(x$method != "lwr.bnd"){
		plot(x$dist.out)
	  graphics::abline(h= {x$c.n /sqrt(x$n)}, col="red", lwd= 1.3)
	  graphics::abline(v= x$alp.hat, lty=2, lwd =1, col="black")
	} else {
		plot(x$dist.out)
		graphics::abline(h = 0.6792/sqrt(x$n), col="red", lwd= 1.3)
		graphics::abline(v= x$alp.Lwr, lty=2, lwd =1, col="black")
	}
}

