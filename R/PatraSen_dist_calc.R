#' Compute the distance to be minimized using Patra and Sen estimation technique in admixture models
#'
#' Compute the distance to be minimized using Patra and Sen estimation technique by integrating along some given grid
#' the appropriate distance. For further developments, see 'Details' below.
#'
#' @param data Sample where the known component density of the admixture model has been transformed into a Uniform(0,1) distribution.
#' @param gridsize Gridsize to make the computations.
#'
#' @details See Patra, R.K. and Sen, B. (2016); Estimation of a Two-component Mixture Model with Applications to Multiple Testing;
#'          JRSS Series B, 78, pp. 869--893.
#'
#' @return a list containing the evaluated distance and some additional information.
#'
#' @examples
#' comp.dist <- list(f = 'norm', g = 'norm')
#' comp.param <- list(f = list(mean = 3, sd = 0.5), g = list(mean = 0, sd = 1))
#' data1 <- rsimmix(n = 3000, unknownComp_weight = 0.6, comp.dist, comp.param)[['mixt.data']]
#' data1_transfo <- knownComp_to_uniform(data = data1, comp.dist = list(comp.dist$f, comp.dist$g),
#'                                       comp.param = list(comp.param$f, comp.param$g))
#' PatraSen_dist_calc(data = data1_transfo, gridsize = 200)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

PatraSen_dist_calc <- function(data, gridsize = 200)
{
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
		# print(i)
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
	class(ret) <- "dist.fun"
	ret$call <- match.call()

	return(ret)
}

print.dist.fun <- function(x,...){
	cat("Call:\n")
	print(x$call)
	t.mat <- cbind(x$grid.pts , x$distance)
	colnames(t.mat) <- c("gamma", "distance")
	print(t(t.mat))
}

plot.dist.fun <- function(x,...){
	plot(x$grid.pts, x$distance, ylab = "distance", xlab = "gamma" , type = "l")
}
