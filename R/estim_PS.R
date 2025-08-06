#' Estimates in an admixture using Patra and Sen approach
#'
#' Estimation of both the weight and the distribution of the unknown component in an admixture model, by Patra and Sen approach.
#' Remind that the admixture probability density function (pdf) l is given by
#'          l = p*f + (1-p)*g,
#' where g is the known component of the two-component mixture, p is the unknown proportion of the unknown component distribution f.
#' More information in 'Details' below concerning the estimation method.
#'
#' @param samples Sample to be studied.
#' @param admixMod An object of class \link[admix]{admix_model}, containing information about the known component distribution and its parameter(s).
#' @param method One of 'lwr.bnd', fixed' or 'cv': depending on whether compute some lower bound of the mixing proportion, the estimate
#'               based on the value of 'c.n' or use cross-validation for choosing 'c.n' (tuning parameter).
#' @param c.n (default to NULL) A positive number for the penalization, see reference below. If NULL, equals to 0.1*log(log(n)).
#' @param folds (optional, default to 10) Number of folds used for cross-validation.
#' @param reps (optional, default to 1) Number of replications for cross-validation.
#' @param cn.s (optional) A sequence of 'c.n' to be used for cross-validation (vector of values). Default is equally
#'            spaced grid of 100 values between .001 x log(log(n)) and 0.2 x log(log(n)).
#' @param cn.length (optional, default to 100) Number of equally spaced tuning parameter (between .001 x log(log(n)) and 0.2 x log(log(n))).
#'                  Values to search from.
#' @param gridsize (default to 600) Number of equally spaced points (between 0 and 1) to evaluate the distance function.
#'                 Larger values are more computationally intensive but also lead to more accurate estimates.
#'
#' @references
#' \insertRef{PatraSen2016}{admix}
#'
#' @return An object of class \link[admix]{estim_PS}, containing 10 attributes: 1) the number of samples studied (1 in this case); 2) the sample
#'         size; 3) the information about component distributions of the admixture model; 4) the estimation method 5patra and Sen here);
#'         5) the estimated mixing weight (estimate of the unknown component proportion); 6) the estimated decontaminated CDF;
#'         7) an object of the class 'dist.fun' (that gives the distance); 8) the tuning parameter 'c.n'; 9) the lower bound of the
#'         estimated mixing proportion (if such an option has been chosen); 10) the number of observations.
#'
#' @seealso [print.estim_PS()] for printing a short version of the results from this estimation method,
#'          and [summary.estim_PS()] for more comprehensive results.
#'
#' @examples
#' \dontrun{
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 800, weight = 0.33,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 2, "scale" = 0.5),
#'                                         list("rate" = 0.25)))
#' data1 <- get_mixture_data(mixt1)
#' ## Define the admixture model:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Estimation step:
#' ex <- estim_PS(samples = data1, admixMod = admixMod1, method = 'fixed')
#' print.estim_PS(ex)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @keywords internal

estim_PS <- function(samples, admixMod, method = c("fixed", "lwr.bnd", "cv"),
                     c.n =  0.1*log(log(length(samples))), folds = 10, reps = 1,
                     cn.s = NULL, cn.length = 100, gridsize = 1200)
{
  if (!inherits(x = admixMod, what = "admix_model"))
    stop("Argument 'admixMod' is not correctly specified. See ?admix_model.")

	if (!is.vector(samples)) stop("'samples' has to be a numerical vector.")
	if (is.null(method)) stop("'method' cannot be NULL")
  n <- length(samples)
  ## Transform the known component distribution into a Uniform distribution:
  samples <- knownComp_to_uniform(samples, admixMod)

  method <- match.arg(method)
	if (method == 'lwr.bnd') {
		dist.out <- PS_dist_calc(samples, gridsize = gridsize)
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
			warning("\n'c.n' is not given. Fixing it to be '0.1*log(log(n))\n")
			c.n <- 0.1 *log(log(n))
		}
		dist.out <- PS_dist_calc(samples, gridsize = gridsize)
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
		out.cv <- estimCV_PS(samples, admixMod, folds = folds, reps = reps,
		                     cn.s = cn.s, cn.length = cn.length, gridsize = gridsize)
		alp.hat <- out.cv$alp.hat
		Fs.hat.fun <- out.cv$Fs.hat
		dist.out <- out.cv$dist.out
		c.n <- out.cv$cn.cv
	}

	ret <- list(
	  n_populations = 1,
	  population_sizes = length(samples),
	  admixture_models = admixMod,
	  estimation_method = "Patra and Sen",
	  estimated_mixing_weights = alp.hat,
	  Fs.hat = Fs.hat.fun,
	  dist.out = dist.out,
	  c.n = c.n,
	  alp.Lwr = alp.Lwr,
	  n = n,
	  data = samples,
	  data.name = deparse1(substitute(samples))
	)

	if (method == "cv"){ ret$cv.out <- out.cv
	} else { ret$cv.out <- NULL }

	ret$method <- method
	ret$call <- match.call()
	class(ret) <- c("estim_PS", "admix_estim")
	return(ret)
}

#' Print method for objects of class 'estim_PS'
#'
#' Print all the attributes of objects of class 'estim_PS'. Results of estimated quantities in an admixture
#' using Patra and Sen approach
#'
#' @param x An object of class 'estim_PS'.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @keywords internal

print.estim_PS <- function(x, ...){
  #cat("\n")
  #cat("Call:")
  #print(x$call)
  cat("\n")
  if(x$method != "lwr.bnd"){
    cat(paste(" Estimated mixing weight (of the unknown component):" , round(x$estimated_mixing_weights,3)))
    cat("\n", paste("The chosen value c_n is", round(x$c.n, 3)), "\n")
#    if( !is.null(x$cv.out)){
#      old_par <- graphics::par()$mfrow
#      on.exit(graphics::par(old_par))
#      graphics::par(mfrow=c(1,2))
#      plot(x$cv.out)
#    }
#    plot(x$dist.out)
  } else if(x$method == 'lwr.bnd'){
    plot(x$dist.out)
    print (paste("The  '95%' lower confidence for alp_0 is ", x$alp.Lwr))
  }
  cat("\n")
}


#' Summary method for objects 'estim_PS'
#'
#' Summarizes the results stored in an object of class 'estim_PS'.
#'
#' @param object An object of class 'estim_PS'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @keywords internal

summary.estim_PS <- function(object, ...)
{
  cat("Call:")
  print(object$call)
  cat("\n")
  cat("------- Sample characteristics -------\n")
  cat("Sample size: ", object$population_sizes, "\n")
  cat("-> Distribution of the known component:", object$admixture_models$comp.dist$known, "\n", sep="")
  cat("-> Parameter(s) of the known component:", paste(names(object$admixture_models$comp.param$known), object$admixture_models$comp.param$known, collapse="\t", sep="="), sep="")
  cat("\n")
  cat("\n------- Estimation results -------\n")
  cat(paste("Estimate of the mixing weight (proportion of the unknown component distribution) is" , round(object$estimated_mixing_weights,3)))
  cat("\n", paste(" The chosen value c_n is", round(object$c.n,3)))
  cat("\n")
}

#plot.estim_PS <- function(x, ...){
#  if(x$method != "lwr.bnd"){
#    plot(x$dist.out$grid.pts)
#    graphics::abline(h= {x$c.n /sqrt(x$n)}, col="red", lwd= 1.3)
#    graphics::abline(v= x$alp.hat, lty=2, lwd =1, col="black")
#  } else {
#    plot(x$dist.out$grid.pts)
#    graphics::abline(h = 0.6792/sqrt(x$n), col="red", lwd= 1.3)
#    graphics::abline(v= x$alp.Lwr, lty=2, lwd =1, col="black")
#  }
#}


#' Cross-validation estimates in an admixture using Patra and Sen approach
#'
#' Estimation of both the unknown component weight and the unknown component distribution in an admixture model,
#' by Patra and Sen (2016) method. Remind that an admixture probability density function l is given by
#'          l = p*f + (1-p)*g,
#' where l is observed, g is the known component, and p is the unknown proportion of the unknown component distribution f.
#' The estimated unknown component weight p is selected using a cross-validation technique that helps to choose the right
#' penalization, see 'Details' below for further information.
#'
#' @param data Sample where the known component density of the admixture model has been transformed into a Uniform(0,1) distribution.
#' @param admixMod An object of class 'admix_model', containing the information about the known component distribution and its parameter(s).
#' @param folds (default to 10) Number of folds used for cross-validation.
#' @param reps (default to 1) Number of replications for cross-validation.
#' @param cn.s (default to NULL) A sequence of 'c.n' to be used for cross-validation (vector of values).
#' @param cn.length (default to NULL) Number of equally spaced tuning parameter (between .001 x log(log(n)) and 0.2 x log(log(n))).
#'                  Values to search from.
#' @param gridsize (default to 200) Number of equally spaced points (between 0 and 1) to evaluate the distance function.
#'                 Larger values are more computationally intensive but also lead to more accurate estimates.
#'
#' @references
#' \insertRef{PatraSen2016}{admix}
#'
#' @return A list containing 'alp.hat' (estimate of the unknown component weight), 'Fs.hat' (list with elements 'x' and 'y' values for the function estimate
#'         of the unknown cumultaive distribution function), 'dist.out' which is an object of the class 'dist.fun'
#'         using the complete data.gen, 'c.n' the value of the tuning parameter used to compute the final estimate,
#'         and finally 'cv.out' which is an object of class 'cv.mixmodel'. The object is NULL if method is "fixed".
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1000, weight = 0.2,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#'
#' ## Define the admixture model:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Estimate the proportion of the unknown component of the admixture model:
#' estimCV_PS(data = data1, admixMod = admixMod1, folds = 10,
#'            reps = 1, cn.s = NULL, cn.length = 3, gridsize = 100)$alp.hat
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

estimCV_PS <- function(data, admixMod, folds = 10, reps = 1, cn.s = NULL, cn.length = NULL, gridsize = 200)
{
  if (!is.vector(data)) stop("'data' has to be a vector.")
  n <- length(data)
  # cvFolds(length(data), K = folds, R = reps, type = "random")

  ## Transform the known component distribution into a Uniform distribution:
  data <- knownComp_to_uniform(data, admixMod)

  if (is.null(cn.s)) {
    if (is.null(cn.length)) stop("Both 'cn.s' and 'cn.length' can not be null")
    cn.s <- seq(.001 * log(log(n)), .2* log(log(n)), length.out = cn.length)
  } else if (length(cn.s) != cn.length) {
    stop("Length of 'cn.s' is different  from cn.length")
  }

  score <-  rep(0, length(cn.s))
  for (rr in 1:reps) {
    cv.ind <- rep(1:folds, times = ceiling(n/folds) )
    cv.ind <- cv.ind[sample.int(length(cv.ind))]
    cv.ind <- cv.ind[1:n]
    for (kk in 1:folds) {
      t.1 <- Sys.time()
      t.data <- PS_dist_calc(data[cv.ind != kk], gridsize = gridsize)
      for (cc in 1:length(cn.s)) {
        score[cc] <- score[cc] +  cv.score(t.data, test.data = data[cv.ind == kk], c.n = cn.s[cc] )
      }
      t.2 <- Sys.time()
      #			if (kk == 1) {
      #				print("Expected time for completion")
      #				print( folds*reps*(t.2-t.1))
      #			}
    }
  }
  cn.cv <- cn.s[which.min(score)]
  tot.out <- PS_dist_calc(data, gridsize = gridsize)
  alp.hat.cv <- sum(tot.out$distance > cn.cv / sqrt(tot.out$n)) / tot.out$gridsize
  if (alp.hat.cv > 0) {
    F.hat <- (tot.out$F.n - (1-alp.hat.cv) * tot.out$F.b) / alp.hat.cv    # computes the naive estimator of F_s
    Fs.hat <- Iso::pava(F.hat,tot.out$Freq, decreasing=FALSE)             # computes the Isotonic Estimator of F_s
    Fs.hat[which(Fs.hat <= 0)] <- 0
    Fs.hat[which(Fs.hat >= 1)] <- 1
  } else {
    Fs.hat <- tot.out$F.b
  }
  Fs.hat.fun <- NULL
  Fs.hat.fun$y <- Fs.hat
  Fs.hat.fun$x <- tot.out$F.n.x

  ret <-  list(cn.cv = cn.cv,
               score = score,
               cn.s = cn.s,
               Fs.hat = Fs.hat.fun,
               alp.hat = alp.hat.cv,
               folds = folds,
               rep = rep,
               data = data,
               dist.out = tot.out)
  ret$call <- match.call()
  class(ret)  <- "estimCV_PS"

  return(ret)
}


#' Print object of class 'estimCV_PS'
#'
#' @param x An object of class 'estimCV_PS'.
#' @param ... Arguments to be passed to print default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

print.estimCV_PS <- function(x,...)
{
  cat("Call:")
  print(x$call)
  print(paste("Cross validated estimate of the mixing weight (proportion of the unknown component distribution) is" , x$alp.hat))
  print(paste(" The cross-validated choice of c_n is", x$cn.cv))
  cat("\n")
}


#' Plot object of class 'estimCV_PS'
#'
#' @param x An object of class 'estimCV_PS'.
#' @param ... Arguments to be passed to plot default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

plot.PS_estimCV <- function(x,...)
{
  plot(x$cn.s, x$score, ylab= "cross validation score", xlab= "c_n" )
}


#' Function that computes the cross-validation score
#'
#' @param tr.data An object of class 'PS_dist_fun'.
#' @param test.data The data used in the test sample.
#' @param c.n The coefficient in the penalization.
#'
#' @references
#' \insertRef{PatraSen2016}{admix}
#'
#' @return The loss corresponding to the cross-validation error.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

cv.score <- function(tr.data, test.data, c.n)
{
  if(!inherits(tr.data, "PS_dist_fun")) stop("'tr.data' is of the wrong type.")

  alp.hat <- sum(tr.data$distance>c.n/sqrt(tr.data$n))/tr.data$gridsize
  if (alp.hat > 0) {
    F.hat <- (tr.data$F.n - (1-alp.hat)*tr.data$F.b) / alp.hat      # computes the naive estimator of F_s
    F.is <- Iso::pava(F.hat, tr.data$Freq, decreasing=FALSE)        # computes the Isotonic Estimator of F_s
    F.is[which(F.is <= 0)] <- 0
    F.is[which(F.is >= 1)] <- 1
    # Fhat.train = alp.hat * F.is + (1-alp.hat) * tr.data$F.b
  } else {
    F.is <- tr.data$F.b
  }

  test.data <- sort(test.data) ## Sorts the data set
  test.data.1 <- unique(test.data) ## Finds the unique data points
  test.Fn <- stats::ecdf(test.data) ## Computes the empirical DF of the data
  test.Fn.1 <- test.Fn(test.data.1)
  test.Fb <- test.data.1
  test.Freq <- diff(c(0,test.Fn.1))

  F.is.interpolate <- stats::approx(x = tr.data$F.b, y = F.is, xout = test.Fb, yleft = 0, yright = 1)
  loss <- sum( {alp.hat* F.is.interpolate$y + (1- alp.hat) *test.Fb - test.Fn.1}^2* test.Freq)
  return(loss)
}


#' Distance to minimize for estimation in an admixture, using Patra and Sen approach
#'
#' Computes the distance to be minimized within the estimation process using Patra and Sen estimation technique by integrating along some given grid
#' the appropriate distance. For further developments, see 'Details' below.
#'
#' @param data Sample where the known component density of the admixture model has been transformed into a Uniform(0,1) distribution.
#' @param gridsize Gridsize to make the computations.
#'
#' @references
#' \insertRef{PatraSen2016}{admix}
#'
#' @return a list containing the evaluated distance and some additional information.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1000, weight = 0.6,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' ## Define the admixture model:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Transform the known component distribution into a Uniform distribution:
#' data <- knownComp_to_uniform(data1, admixMod1)
#' PS_dist_calc(data = data, gridsize = 200)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

PS_dist_calc <- function(data, gridsize = 200)
{
  ## Transform the known component distribution into a Uniform distribution:
  #data <- knownComp_to_uniform(data, admixMod)

  q <- 0.6792
  n <- length(data)         # length of the data set
  data <- sort(data)        # sorts the data set
  data.1 <- unique(data)    # finds the unique data points
  Fn <- stats::ecdf(data)   # computes the empirical DF of the data
  Fn.1 <- Fn(data.1)        # empirical DF of the data at the data points
  ## Calculate the known F_b at the data points
  ## Note: for Uniform(0,1) F_b(x) = x
  ## Usually would need to CHANGE this
  Fb <- data.1
  ## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
  Freq <- diff(c(0, Fn.1))
  distance <- rep(0, floor(gridsize * .12))
  distance[0] <- sqrt( t((Fn.1-Fb)^2) %*% Freq )

  grid.pts <- {1:floor(gridsize)} / gridsize

  for (ii in 1:length(grid.pts)) {
    aa <- grid.pts[ii]
    F.hat <- (Fn.1 - (1-aa)*Fb) / aa                # computes the naive estimator of F_s
    F.is <- Iso::pava(F.hat, Freq, decreasing=FALSE)     # computes the Isotonic Estimator of F_s
    F.is[which(F.is <= 0)] <- 0
    F.is[which(F.is >= 1)] <- 1
    distance[ii] <- aa * sqrt( t((F.hat-F.is)^2) %*% Freq )
  }
  # Lower.Cfd.Bound <- sum(distance>q/sqrt(n))/gridsize
  ret <- list(distance = distance,
              gridsize = gridsize,
              grid.pts = grid.pts,
              F.n = Fn.1,
              F.n.x = data.1,
              Freq = Freq,
              F.b = Fb,
              n = n)
  class(ret) <- "PS_dist_fun"
  ret$call <- match.call()

  return(ret)
}

#' Results of the computation of the contrast using Patra and Sen approach
#'
#' Plots all the attributes of objects of class 'PS_dist_fun'.
#'
#' @param x An object of class 'PS_dist_fun'.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

plot.PS_dist_fun <- function(x,...){
  plot(x$grid.pts, x$distance, ylab = "distance", xlab = "gamma" , type = "l")
}


#' Results of the computation of the contrast using Patra and Sen approach
#'
#' Print all the attributes of objects of class 'PS_dist_fun'.
#'
#' @param x An object of class 'PS_dist_fun'.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

print.PS_dist_fun <- function(x,...){
  cat("Call:\n")
  print(x$call)
  t.mat <- cbind(x$grid.pts , x$distance)
  colnames(t.mat) <- c("gamma", "distance")
  print(t(t.mat))
  cat("\n")
}


#' Estimates the density of the unknown component, using Patra and Sen approach
#'
#' Computes by Patra and Sen technique the estimate of the density function when the latter is known to be either decreasing or increasing, still in the framework of an admixture model.
#' Remind that an admixture probability density function (pdf) l is given by
#'          l = p*f + (1-p)*g,
#' where l is observed, g is the known component, and p is the unknown proportion of the unknown component distribution f.
#'
#' @param input an R object of class 'PS_estim' or 'PS_estimCV'.
#' @param dec.density a boolean indicating whether the density is increasing or decreasing.
#'
#' @references
#' \insertRef{PatraSen2016}{admix}
#'
#' @return an estimator of the unknown component density.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1000, weight = 0.6,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' ## Define the admixture model:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Estimate the proportion of the unknown component of the admixture model:
#' res <- estim_PS(data = data1, admixMod = admixMod1, method = 'fixed',
#'          c.n = 0.1*log(log(length(data1))), gridsize = 1000)
#' PS_unknownDensity_estim(res, dec.density = TRUE)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

PS_unknownDensity_estim <- function(input, dec.density = TRUE)
{
  if (inherits(input, "estim_PS") || inherits(input, "estimCV_PS")) {
    ## CDF:
    Fs.hat <- input$Fs.hat
  } else {
    stop("This function only works on objects of class 'estim_PS' or 'estimCV_PS'. See functions 'estim_PS' or 'estimCV_PS'.")
  }
  if (dec.density == TRUE){
    ll <- fdrtool::gcmlcm(Fs.hat$x, Fs.hat$y, type = "lcm")
  } else if (dec.density == FALSE){
    ll <- fdrtool::gcmlcm(Fs.hat$x, Fs.hat$y, type = "gcm")
  }

  ## Density function:
  fs.hat <- NULL
  fs.hat$x <- rep(ll$x.knots, each = 2)                 # data points for density
  fs.hat$y <- c(0, rep(ll$slope.knots, each = 2), 0)    # value of density

  return(fs.hat)
}
