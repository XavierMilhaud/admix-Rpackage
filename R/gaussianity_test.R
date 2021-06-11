#' One-sample gaussianity test in admixture models using Bordes and Vandekerkhove estimation method
#'
#' Perform the hypothesis test to know whether the unknown mixture component is gaussian or not, knowing that the known one
#' has support on the real line (R). However, the case of non-gaussian known component can be overcome thanks to the basic
#' transformation by cdf. Recall that an admixture model has probability density function (pdf) l = p*f + (1-p)*g, where g is
#' the known pdf and l is observed (others are unknown). Requires optimization (to estimate the unknown parameters) as defined
#' by Bordes & Vandekerkhove (2010), which means that the unknown mixture component must have a symmetric density.
#'
#' @param sample1 Observed sample with mixture distribution given by l = p*f + (1-p)*g, where f and p are unknown and g is known.
#' @param comp.dist List with two elements corresponding to the component distributions involved in the admixture model. Unknown
#'                  elements must be specified as 'NULL' objects. For instance if 'f' is unknown: list(f = NULL, g = 'norm').
#' @param comp.param List with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R names for distributions.
#'                   Unknown elements must be specified as 'NULL' objects (e.g. if 'f' is unknown: list(f=NULL, g=list(mean=0,sd=1)).
#' @param K Number of coefficients considered for the polynomial basis expansion.
#' @param lambda Rate at which the normalization factor is set in the penalization rule for model selection (in ]0,1/2[). See 'Details' below.
#' @param support Support of the densities under consideration, useful to choose the polynomial orthonormal basis. One of 'Real',
#'                'Integer', 'Positive', or 'Bounded.continuous'.
#'
#' @details See the paper 'False Discovery Rate model Gaussianity test' (Pommeret & Vanderkerkhove, 2017).
#'
#' @return A list of 6 elements, containing: 1) the rejection decision; 2) the p-value of the test; 3) the test statistic; 4) the
#'         variance-covariance matrix of the test statistic; 5) the selected rank for testing; and 6) a list of the estimates
#'         (unknown component weight 'p', shift location parameter 'mu' and standard deviation 's' of the symmetric unknown distribution).
#'
#' @examples
#' ####### Under the null hypothesis H0.
#' ## Parameters of the gaussian distribution to be tested:
#' list.comp <- list(f = "norm", g = "norm")
#' list.param <- list(f = c(mean = 2, sd = 0.5),
#'                    g = c(mean = 0, sd = 1))
#' ## Simulate and plot the data at hand:
#' obs.data <- rsimmix(n = 150, unknownComp_weight = 0.9, comp.dist = list.comp,
#'                     comp.param = list.param)[['mixt.data']]
#' plot(density(obs.data))
#' ## Performs the test:
#' list.comp <- list(f = NULL, g = "norm")
#' list.param <- list(f = NULL, g = c(mean = 0, sd = 1))
#' gaussianity_test(sample1 = obs.data, comp.dist = list.comp, comp.param = list.param,
#'                        K = 3, lambda = 0.1, support = 'Real')
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

gaussianity_test <- function(sample1, comp.dist, comp.param, K = 3, lambda = 0.2,
                             support = c('Real','Integer','Positive','Bounded.continuous'))
{
  if ( (length(comp.dist) != 2) | (length(comp.param) != 2) ) stop("Arguments 'comp.dist' and/or 'comp.param' were not correctly specified")
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
  }
  ## Extracts the information on component distributions and stores in expressions:
  comp.dens <- sapply(X = paste0("d",comp.dist), FUN = get, pos = "package:stats", mode = "function")
  assign(x = names(comp.dens)[2], value = comp.dens[[2]])
  comp.sim <- sapply(X = paste0("r",comp.dist), FUN = get, pos = "package:stats", mode = "function")
  assign(x = names(comp.sim)[2], value = comp.sim[[2]])
  ## Creates the expression allowing further to compute the theoretical 2nd-order moment of the known component:
  expr.dens <- paste(names(comp.dens)[2],"(x,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  expr.sim <- paste(names(comp.sim)[2],"(n=1000000,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
  ## Accurate approximation of the second-order moment of the known component:
  m1_knownComp <- mean(eval(parse(text = expr.sim)))
  m2_knownComp <- mean(eval(parse(text = expr.sim))^2)

	##-------- Split data sample -----------##
	## Random sampling (p.6 Pommeret & Vandekerkhove) : create subsamples to get uncorrelated estimators of the different parameters
	n <- length(sample1)
	blocksize <- n %/% 3							                   # block size
	values <- stats::runif(n)									           # simulate random values
	rang <- rank(values)									               # associate a rank to each insured / observation
	block <- (rang-1) %/% blocksize + 1						       # associate to each individual a block number
	block <- as.factor(block)
	data.coef <- sample1[block == 1]
	n.coef <- length(data.coef)	                         # estimation of empirical moments
	data.p <- sample1[block == 2]
	n.p <- length(data.p)	                               # estimation of component weights param de poids des composantes 'p'
	data.loc <- sample1[block == 3]
	n.loc <- length(data.loc)	                           # estimation of localisation parameter 'mu' of the symmetric unknown density

	##-------- Estimation of parameters and corresponding variances -----------##
	## Focus on parameters (weight, localization and variance), consider independent subsamples of the original data:
	hat_p <- BVdk_estimParam(data = data.p, method = "L-BFGS-B", comp.dist = comp.dist, comp.param = comp.param)[1]
	hat_loc <- BVdk_estimParam(data = data.loc, method = "L-BFGS-B", comp.dist = comp.dist, comp.param = comp.param)[2]
  ## Plug-in method to estimate the variance:
	kernelDensity_est_obs <- stats::density(sample1)
	integrand_totweight <- function(x) {
	  w <- (1/hat_p) * stats::approxfun(kernelDensity_est_obs)(x) - ((1-hat_p)/hat_p) * eval(parse(text = expr.dens))
	  return( w * (w > 0) )
	}
	total_weight <- stats::integrate(integrand_totweight, lower = min(kernelDensity_est_obs$x), upper = max(kernelDensity_est_obs$x))
	integrand_unknownComp <- function(x) {
	  w <- (1/hat_p) * stats::approxfun(kernelDensity_est_obs)(x) - ((1-hat_p)/hat_p) * eval(parse(text = expr.dens))
	  return( (x-hat_loc)^2 * (w * (w > 0)) / total_weight$value )
	}
	hat_s2 <- stats::integrate(integrand_unknownComp, lower = min(kernelDensity_est_obs$x), upper = max(kernelDensity_est_obs$x))$value
	#hat_s2 <- (1/hat_p) * ( mean(sample1^2) - ((1-hat_p) * m2_knownComp) ) - (1/hat_p^2) * (mean(sample1)-(1-hat_p)*m1_knownComp)^2
	## Then on the variances of the estimators: semiparametric estimation (time-consuming), based on results by Bordes & Vandekerkhove (2010)
	varCov <- BVdk_varCov_estimators(data = sample1, loc = hat_loc, p = hat_p, comp.dist = comp.dist, comp.param = comp.param)
	var.hat_p <- varCov[["var_pEstim"]]
	var.hat_loc <- varCov[["var_muEstim"]]
	## Estimation of the variance of the variance estimator:
	#var.hat_s2 <- BVdk_ML_varCov_estimators(data = sample1, hat_w = hat_p, hat_loc = hat_loc, hat_var = hat_s2,
	#                                        comp.dist = comp.dist, comp.param = comp.param)

	##-------- Compute test statistics with data 'data.coef' -----------##
	stat.R <- matrix(rep(NA, n.coef*K), nrow = K, ncol = n.coef)
	## Plug-in (estimation of param. 'mu' et 's') to deduce coef. in Hermite polynomial orthonormal basis:
	coef.unknown <- unlist(lapply(X = orthoBasis_coef(data = stats::rnorm(100000, hat_loc, sqrt(hat_s2)), supp = support, degree = K, m = 3, other = NULL), FUN = mean))
	coef.known <- unlist(lapply(X = orthoBasis_coef(data = eval(parse(text = expr.sim)), supp = support, degree = K, m = 3, other = NULL), FUN = mean))

	## Cf definition of R_kn below formula (12) p.5 :
	poly_basis <- poly_orthonormal_basis(support = support, deg = K, x = data.coef, m = 3)
	var.coef <- numeric(length = K)
	for (i in 1:K) {
	  mixt.coefs <- orthopolynom::polynomial.values(poly_basis, data.coef)[[i+1]] / sqrt(factorial(i))
	  var.coef[i] <- stats::var(mixt.coefs)
	  stat.R[i, ] <- mixt.coefs - hat_p * (coef.unknown[i]-coef.known[i]) - coef.known[i]
	  mixt.coefs <- NULL
	}
	## Mean (at each order) of differences b/w 'theoretical' gaussian coef (hat_loc,schap) in ortho basis and computed coef from data:
	statistics.R <- colMeans(t(stat.R))

	##-------- Scaling (with variance) of the test statistic -----------##
	var.R <- matrix(data = NA, nrow = K, ncol = K)
	## Introduce the correction factor for the adjustment of the variance, given the initial split of the sample:
	w.coef <- sqrt( (n.p * n.loc) / (n.coef + n.p + n.loc)^2 )
	w.p <- sqrt( (n.coef * n.loc) / (n.coef + n.p + n.loc)^2 )
	w.loc <- sqrt( (n.p * n.coef) / (n.coef + n.p + n.loc)^2 )
	## FIXME: this formula has to be checked...
	var.R[1,1] <- w.coef * var.coef[1] + w.p * var.hat_p * (hat_loc - coef.known[1])^2 + w.loc * var.hat_loc * hat_p^2
	var.R[2,2] <- w.coef * var.coef[2] + w.p * var.hat_p * ((1/2)*(hat_s2+hat_loc^2-1) - coef.known[1])^2 + w.loc * var.hat_loc * hat_p^2 * hat_loc^2
	var.R[3,3] <- w.coef * var.coef[2] + w.p * var.hat_p * (1/36) * (3*hat_loc*hat_s2 + hat_loc^2 - 6*hat_loc)^2 +
	              w.loc * var.hat_loc * (1/36) * (hat_p*(3*hat_s2+hat_loc-6))^2

	##-------- Compute the test statistics at different orders -----------##
	test.statistic <- numeric(length = K)
	val <- 0
	for (k in 1:K) {
	  ## Equations (13) and (14) p.5 (except that we already scaled here, cf explanation top of p.6):
	  test.statistic[k] <- val + n^lambda * statistics.R[k]^2 * var.R[k,k]^(-1) - log(n)
		val <- test.statistic[k]
	}
	## Select the right order (optimal number of coefficients needed in the expansion) :
	selected.index <- which.max(test.statistic)
	## Then remove previously introduced penalty in the selection rule to get back to the test statistic (Equation (14) p.5):
	final.stat <- test.statistic[selected.index] + selected.index * log(n)
	p.value <- 1 - stats::pchisq(final.stat, 1)

	## If the test statistic is greater that the quantile of interest, reject the null hypothesis (otherwise do not reject):
	rej <- 0
  if (final.stat > stats::qchisq(0.95,1)) { rej <- rej+1 }

	list(decision = rej, p_value = p.value, test.stat = final.stat, var.stat = var.R, rank = selected.index,
	     estimates = list(p = hat_p, mu = hat_loc, s = sqrt(hat_s2)))
}

