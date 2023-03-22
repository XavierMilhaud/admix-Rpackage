#' Cross-validation estimate (by Patra and Sen) of the unknown component weight as well as the unknown distribution in an admixture model
#'
#' Estimation of unknown elements (by Patra and Sen method) under the admixture model with probability density function l:
#'          l = p*f + (1-p)*g,
#' where g is the known component of the two-component admixture, p is the unknown proportion of the unknown component distribution f.
#' The estimated unknown component weight p is selected using a cross-validation technique that helps to choose the right penalization, see
#' 'Details' below for further information.
#'
#' @param data Sample where the known component density of the admixture model has been transformed into a Uniform(0,1) distribution.
#' @param folds (default to 10) Number of folds used for cross-validation.
#' @param reps (default to 1) Number of replications for cross-validation.
#' @param cn.s (default to NULL) A sequence of 'c.n' to be used for cross-validation (vector of values).
#' @param cn.length (default to NULL) Number of equally spaced tuning parameter (between .001 x log(log(n)) and 0.2 x log(log(n))).
#'                  Values to search from.
#' @param gridsize (default to 200) Number of equally spaced points (between 0 and 1) to evaluate the distance function.
#'                 Larger values are more computationally intensive but also lead to more accurate estimates.
#'
#' @details See Patra, R.K. and Sen, B. (2016); Estimation of a Two-component Mixture Model with Applications to Multiple Testing;
#'          JRSS Series B, 78, pp. 869--893.
#'
#' @return A list containing 'alp.hat' (estimate of the unknown component weight), 'Fs.hat' (list with elements 'x' and 'y' values for the function estimate
#'         of the unknown cumultaive distribution function), 'dist.out' which is an object of the class 'dist.fun'
#'         using the complete data.gen, 'c.n' the value of the tuning parameter used to compute the final estimate,
#'         and finally 'cv.out' which is an object of class 'cv.mixmodel'. The object is NULL if method is "fixed".
#'
#' @examples
#' ## Simulate data:
#' comp.dist <- list(f = 'norm', g = 'norm')
#' comp.param <- list(f = list(mean = 3, sd = 0.5),
#'                    g = list(mean = 0, sd = 1))
#' data1 <- rsimmix(n = 2000, unknownComp_weight = 0.3, comp.dist, comp.param)[['mixt.data']]
#' ## Transform the known component of the admixture model into a Uniform(0,1) distribution:
#' comp.dist <- list(f = NULL, g = 'norm')
#' comp.param <- list(f = NULL, g = list(mean = 0, sd = 1))
#' data1_transfo <- knownComp_to_uniform(data = data1, comp.dist = list(comp.dist$f, comp.dist$g),
#'                                       comp.param = list(comp.param$f, comp.param$g))
#' plot(density(data1_transfo))
#' ## Estimate the proportion of the unknown component of the admixture model:
#' PatraSen_cv_mixmodel(data = data1_transfo, folds = 3, reps = 1, cn.s = NULL,
#'                                cn.length = 3, gridsize = 100)$alp.hat
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

PatraSen_cv_mixmodel <- function(data, folds = 10, reps = 1, cn.s = NULL, cn.length = NULL, gridsize = 200)
{
  if (!is.vector(data)) stop("'data' has to be a vector.")
  warning("Make sure that data is transformed such that the known component is Uniformly(0,1) distributed.")
  n <- length(data)
  # cvFolds(length(data), K = folds, R = reps, type = "random")

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
      t.data <- PatraSen_dist_calc(data[cv.ind != kk], gridsize = gridsize)
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
  tot.out <- PatraSen_dist_calc(data, gridsize = gridsize)
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
  class(ret)  <- "cv.mixmodel"

  return(ret)
}

print.cv.mixmodel <- function(x,...)
{
  cat("Call:")
  print(x$call)
  print(paste("Cross validated estimate of alp is" , x$alp.hat))
  print(paste(" The cross-validated choice of c_n is", x$cn.cv))
}

plot.cv.mixmodel <- function(x,...)
{
  plot(x$cn.s, x$score, ylab= "cross validation score", xlab= "c_n" )
}

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

