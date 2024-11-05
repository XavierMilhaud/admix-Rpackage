#' Gaussianity test in an admixture model
#'
#' Performs an hypothesis test to check for the gaussianity of the unknown mixture component, given that the known component
#' has support on the real line. Recall that an admixture model has probability density function (pdf) l = p*f + (1-p)*g, where g is
#' the known pdf and l is observed (others are unknown). This test requires optimization (to estimate the unknown parameters) as
#' defined by Bordes & Vandekerkhove (2010), which means that the unknown mixture component must have a symmetric density.
#'
#' @param samples Sample under study.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#' @param conf_level (default to 0.95) The confidence level. Equals 1-alpha, where alpha is the level of the test (type-I error).
#' @param ask_poly_param (default to FALSE) If TRUE, ask the user to choose both the order 'K' of expansion coefficients in the
#'                        orthonormal polynomial basis, and the penalization rate 's' involved on the penalization rule for the test.
#' @param K (K > 0, default to 3) If not asked (see the previous argument), number of coefficients considered for the polynomial basis expansion.
#' @param s (in ]0,1/2[, default to 0.25) If not asked (see the previous argument), normalization rate involved in the penalization rule
#'          for model selection. See the reference below.
#' @param support Support of the probability distributions, useful to choose the appropriate polynomial orthonormal basis. One of 'Real',
#'                'Integer', 'Positive', or 'Bounded.continuous'.
#' @param ... Optional arguments to \link[admix]{estim_BVdk}.
#'
#' @details Extensions to the case of non-Gaussian known components can be overcome thanks to basic transformations using cdf.
#'
#' @references
#' \insertRef{PommeretVandekerkhove2019}{admix}
#'
#' @return An object of class 'gaussianity_test', containing 10 elements: 1) the number of populations under study (1 in this case);
#'         2) the sample size; 3) the information about the known component distribution; 4) the reject decision of the test; 5) the
#'         confidence level of the test, 6) the p-value of the test; 7) the value of the test statistic; 8) the variance of the test
#'         statistic at each order in the polynomial orthobasis expansion; 9) the selected rank (order) for the test statistic;
#'         10) a list of estimates (mixing weight, mean and standard deviation of the Gaussian unknown distribution).
#'
#' @examples
#' ####### Under the null hypothesis H0.
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 250, weight = 0.4,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("rate" = 1)))
#' data1 <- getmixtData(mixt1)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Performs the test:
#' gaussianity_test(samples = data1, admixMod = admixMod1,
#'                  conf_level = 0.95, K = 3, s = 0.1, support = "Real")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

gaussianity_test <- function(samples, admixMod, conf_level = 0.95, ask_poly_param = FALSE, K = 3, s = 0.25,
                             support = c('Real','Integer','Positive','Bounded.continuous'), ...)
{
  support <- match.arg(support)

  if (ask_poly_param) {
    K.user <- base::readline("Please enter 'K' (integer), the order for the polynomial expansion in the orthonormal basis: ")
    s.user <- base::readline("Please enter 's' in ]0,0.5[, involved in the penalization rule for model selection where lower values of 's' lead to more powerful tests: ")
  } else {
    K.user <- K
    s.user <- s
  }

  if (any(admixMod$comp.dist == "multinom")) stop("Gaussianity test for those mixture distribution is not supported.\n")

  ## Extract the information on component distributions:
  comp.dist.dens <- paste0("d", admixMod$comp.dist$known)
  comp.dens <- sapply(X = comp.dist.dens, FUN = get, mode = "function")
  assign(x = names(comp.dens)[1], value = comp.dens[[1]])
  comp.dist.sim <- paste0("r", admixMod$comp.dist$known)
  comp.sim <- sapply(X = comp.dist.sim, FUN = get, mode = "function")
  assign(x = names(comp.sim)[1], value = comp.sim[[1]])

  ## Creates the expression allowing further to compute the theoretical 2nd-order moment of the known component:
  expr.dens <- paste(names(comp.dens)[1],"(x,", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")
  expr.sim <- paste(names(comp.sim)[1],"(n=1000000,", paste(names(admixMod$comp.param$known),
                    "=", admixMod$comp.param$known, sep = "", collapse = ","), ")", sep="")

  ## Accurate approximation of the second-order moment of the known component:
  m1_knownComp <- mean(eval(parse(text = expr.sim)))
  m2_knownComp <- mean(eval(parse(text = expr.sim))^2)

	##-------- Split data sample -----------##
	## Random sampling (p.6 Pommeret & Vandekerkhove) : create subsamples to get uncorrelated estimators of the different parameters
	n <- length(samples)
	blocksize <- n %/% 2							                   # block size
	values <- stats::runif(n)									           # simulate random values
	rang <- rank(values)									               # associate a rank to each insured / observation
	block <- (rang-1) %/% blocksize + 1						       # associate to each individual a block number
	block <- as.factor(block)
	data.coef <- samples[block == 1]
	n.coef <- length(data.coef)	                         # for the estimation of empirical moments
	data.BVdk <- samples[block == 2]
	n.BVdk <- length(data.BVdk)	                        # for the estimation of mixing weight / location shift

	##-------- Estimation of parameters and corresponding variances -----------##
	## Focus on parameters (weight, localization and variance), consider independent subsamples of the original data:
	BVdk <- estim_BVdk(samples = data.BVdk, admixMod = admixMod, ...)
	hat_p <- BVdk$estimated_mixing_weights
	hat_loc <- BVdk$estimated_locations
	## Plug-in method to estimate the variance:
	kernelDensity_est_obs <- stats::density(samples)
	integrand_totweight <- function(x) {
	  w <- (1/hat_p) * stats::approxfun(kernelDensity_est_obs)(x) - ((1-hat_p)/hat_p) * eval(parse(text = expr.dens))
	  return( w * (w > 0) )
	}
	integrand_unknownComp <- function(x) {
	  w <- (1/hat_p) * stats::approxfun(kernelDensity_est_obs)(x) - ((1-hat_p)/hat_p) * eval(parse(text = expr.dens))
	  return( (x-hat_loc)^2 * (w * (w > 0)) / total_weight )
	}
	total_weight <- try(stats::integrate(integrand_totweight, lower = min(kernelDensity_est_obs$x), upper = max(kernelDensity_est_obs$x))$value, silent = TRUE)
	hat_s2 <- try(stats::integrate(integrand_unknownComp, lower = min(kernelDensity_est_obs$x), upper = max(kernelDensity_est_obs$x))$value, silent = TRUE)

	if (inherits(x = total_weight, what = "try-error", which = FALSE)) {
	  total_weight <- cubature::cubintegrate(integrand_totweight,
	                                         lower = min(kernelDensity_est_obs$x),
	                                         upper = max(kernelDensity_est_obs$x),
	                                         method = "pcubature")$integral
	}
	if (inherits(x = hat_s2, what = "try-error", which = FALSE)) {
	  hat_s2 <- cubature::cubintegrate(integrand_unknownComp,
	                                   lower = min(kernelDensity_est_obs$x),
	                                   upper = max(kernelDensity_est_obs$x),
	                                   method = "pcubature")$integral
	}
  #hat_s2 <- (1/hat_p) * ( mean(samples^2) - ((1-hat_p) * m2_knownComp) ) - (1/hat_p^2) * (mean(samples)-(1-hat_p)*m1_knownComp)^2

	## Then on the variances of the estimators: semiparametric estimation (time-consuming), based on results by Bordes & Vandekerkhove (2010)
	varCov <- BVdk_varCov_estimators(estim = BVdk, data = samples, admixMod = admixMod)
	var.hat_p <- varCov$var.estim_prop
	var.hat_loc <- varCov$var.estim_location
	## Estimation of the variance of the variance estimator:
	#var.hat_s2 <- BVdk_ML_varCov_estimators(data = samples, hat_w = hat_p, hat_loc = hat_loc, hat_var = hat_s2,
	#                                        comp.dist = comp.dist, comp.param = comp.param)

	##-------- Compute test statistics with data 'data.coef' -----------##
	stat.R <- matrix(rep(NA, n.coef*K.user), nrow = K.user, ncol = n.coef)
	## Plug-in (estimation of param. 'mu' et 's') to deduce coef. in Hermite polynomial orthonormal basis:
	coef.known <- unlist(lapply(X = orthoBasis_coef(data = eval(parse(text = expr.sim)), supp = support, degree = K.user, m = 3, other = NULL), FUN = mean))
	## Gaussianity test implies to assume a gaussian distribution of the unknown component:
	coef.unknown <- unlist(lapply(X = orthoBasis_coef(data = stats::rnorm(100000, hat_loc, sqrt(hat_s2)), supp = support, degree = K.user, m = 3, other = NULL), FUN = mean))

	## Cf definition of R_kn below formula (12) p.5 :
	poly_basis <- poly_orthonormal_basis(support = support, deg = K.user, x = data.coef, m = 3)
	var.coef <- numeric(length = K.user)
	for (i in 1:K.user) {
	  mixt.coefs <- orthopolynom::polynomial.values(poly_basis, data.coef)[[i+1]] / sqrt(factorial(i))
	  var.coef[i] <- stats::var(mixt.coefs)
	  stat.R[i, ] <- mixt.coefs - hat_p * (coef.unknown[i]-coef.known[i]) - coef.known[i]
	  mixt.coefs <- NULL
	}
	## Mean (at each order) of differences b/w 'theoretical' gaussian coef (hat_loc,schap) in ortho basis and computed coef from data:
	statistics.R <- colMeans(t(stat.R))

	##-------- Scaling (with variance) of the test statistic -----------##
	var.R <- matrix(data = NA, nrow = K.user, ncol = K.user)
	## Introduce the correction factor for the adjustment of the variance, given the initial split of the sample:
	#w.coef <- sqrt( (n.p * n.loc) / (n.coef + n.p + n.loc)^2 )
	#w.p <- sqrt( (n.coef * n.loc) / (n.coef + n.p + n.loc)^2 )
	#w.loc <- sqrt( (n.p * n.coef) / (n.coef + n.p + n.loc)^2 )
	w.coef <- sqrt( (n.BVdk^2) / (n.coef + n.BVdk)^2 )
	w.p <- sqrt( ((n.coef/2)^2) / (n.coef + n.BVdk)^2 )
	w.loc <- sqrt( ((n.coef/2)^2) / (n.coef + n.BVdk)^2 )
	## FIXME: check this formula!
	var.R[1,1] <- w.coef * var.coef[1] + w.p * var.hat_p * (hat_loc - coef.known[1])^2 + w.loc * var.hat_loc * hat_p^2
	var.R[2,2] <- w.coef * var.coef[2] + w.p * var.hat_p * ((1/2)*(hat_s2+hat_loc^2-1) - coef.known[1])^2 + w.loc * var.hat_loc * hat_p^2 * hat_loc^2
	var.R[3,3] <- w.coef * var.coef[2] + w.p * var.hat_p * (1/36) * (3*hat_loc*hat_s2 + hat_loc^2 - 6*hat_loc)^2 +
	              w.loc * var.hat_loc * (1/36) * (hat_p*(3*hat_s2+hat_loc-6))^2

	##-------- Compute the test statistics at different orders -----------##
	test.statistic <- numeric(length = K.user)
	val <- 0
	for (k in 1:K.user) {
	  ## Equations (13) and (14) p.5 (except that we already scaled here, cf explanation top of p.6):
	  test.statistic[k] <- val + n^s.user * statistics.R[k]^2 * var.R[k,k]^(-1) - log(n)
		val <- test.statistic[k]
	}
	## Select the right order (optimal number of coefficients needed in the expansion) :
	selected.index <- which.max(test.statistic)
	## Then remove previously introduced penalty in the selection rule to get back to the test statistic (Equation (14) p.5):
	final.stat <- test.statistic[selected.index] + selected.index * log(n)
	p.value <- 1 - stats::pchisq(final.stat, 1)

	## If the test statistic is greater that the quantile of interest, reject the null hypothesis (otherwise do not reject):
	rej <- FALSE
  if (final.stat > stats::qchisq(conf_level,1)) { rej <- TRUE }

	obj <- list(
	  n_populations = 1,
	  population_sizes = length(samples),
	  admixture_models = admixMod,
	  reject_decision = rej,
	  confidence_level = conf_level,
	  p_value = p.value,
	  test_statistic_value = final.stat,
	  var_statistic = var.R,
	  selected_rank = selected.index,
	  estimates = list(p = hat_p, mu = hat_loc, s = sqrt(hat_s2))
	  )
	class(obj) <- c("gaussianity_test", "admix_test")
	obj$call <- match.call()
	return(obj)
}



#' Print method for objects 'gaussianity_test'
#'
#' @param x An object of class 'gaussianity_test'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.gaussianity_test <- function(x, ...)
{
  cat("Call:")
  print(x$call)
  cat("\n")
  cat("Is the null hypothesis (gaussian unknown component distribution) rejected? ",
      ifelse(x$reject_decision, "Yes", "No"), sep="")
  cat("\nTest p-value: ", round(x$p_value,3), "\n", sep="")
}


#' Summary method for objects 'gaussianity_test'
#'
#' @param object An object of class 'gaussianity_test'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.gaussianity_test <- function(object, ...)
{
  cat("Call:")
  print(object$call)
  cat("\n--------- About samples ---------\n")
  cat(paste("Size of sample ", 1:object$n_populations, ": ", object$population_sizes, sep = ""), sep = "\n")
  cat("\n-------- About contamination (admixture) models -------")
  cat("-> Distribution and parameters of the known component \n for the admixture model: ", sep="")
  cat(object$admixture_models$comp.dist$known, "\n")
  print(unlist(object$admixture_models$comp.param$known, use.names = TRUE))
  cat("\n------- Test decision -------\n")
  cat("Is the null hypothesis (gaussian unknown component distribution) rejected? ",
      ifelse(object$reject_decision, "Yes", "No"), sep="")
  cat("\nConfidence level of the test: ", object$confidence_level, sep="")
  cat("\nTest p-value: ", round(object$p_value,3), sep="")
  cat("\n\n------- Test statistic -------\n")
  cat("Selected rank of the test statistic (following the penalization rule): ", object$selected_rank, sep="")
  cat("\nValue of the test statistic ", round(object$test_statistic_value,2), sep="")
  cat("\nVariance of the test statistic (for each order of the expansion):", paste(round(diag(object$var_statistic),2), sep=" "), sep=" ")
  cat("\n\n------- Estimates -------\n")
  cat(paste(c(" Estimation of the mixing weight (proportion of the unknown component distribution): ",
          "Estimation of the mean of the unknown gaussian distribution: ",
          "Estimation of the standard deviation of the unknown gaussian distribution: "),
          sapply(object$estimates, round, 2), "\n"))
}
