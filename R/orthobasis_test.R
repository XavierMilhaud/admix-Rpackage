#' Equality test of two unknown component distributions using polynomial expansions
#'
#' Tests the null hypothesis (H0: f1=f2) using the decomposition of unknown component densities of two admixture distributions in
#' an adequate orthonormal polynomial basis. Recall that we have two admixture models with respective probability density
#' functions (pdf) l1 = p1*f1 + (1-p1)*g1 and l2 = p2*f2 + (1-p2)*g2, where g1 and g2 are the only known elements and l1 and l2
#' are observed. The admixture weights p1 and p2 thus have to be estimated. For further information on this method, see 'Details' below.
#'
#' @param samples List of the two samples, each one following the mixture distribution given by l = p*f + (1-p)*g,
#'                with f and p unknown and g known.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#' @param conf_level The confidence level, default to 95 percent. Equals 1-alpha, where alpha is the level of the test (type-I error).
#' @param est_method Estimation method to get the component weights, either 'PS' (Patra and Sen estimation) or 'BVdk'
#'                   (Bordes and Vendekerkhove estimation). Choosing 'PS' requires to specify the number of bootstrap samples.
#' @param ask_poly_param (default to FALSE) If TRUE, ask the user to choose both the order 'K' of expansion coefficients in the
#'                        orthonormal polynomial basis, and the penalization rate 's' involved on the penalization rule for the test.
#' @param K (K > 0, default to 3) If not asked (see the previous argument), number of coefficients considered for the polynomial basis expansion.
#' @param s (in ]0,1/2[, default to 0.49) If not asked (see the previous argument), rate at which the normalization factor is set in
#'           the penalization rule for model selection (in ]0,1/2[). Low values of 's' favors the detection of alternative hypothesis.
#'           See reference below.
#' @param nb_echBoot (default to 100) Number of bootstrap samples, useful when choosing 'PS' estimation method.
#' @param support Support of the probability distributions, useful to choose the appropriate polynomial orthonormal basis. One of 'Real',
#'                'Integer', 'Positive', or 'Bounded.continuous'.
#' @param bounds_supp (default to NULL) Useful if support = 'Bounded.continuous', a list of minimum and maximum bounds, specified as
#'                     following: list( list(min.f1,min.g1,min.f2,min.g2) , list(max.f1,max.g1,max.f2,max.g2) )
#' @param ... Optional arguments to \link[admix]{estim_BVdk} or \link[admix]{estim_PS}, depending on the chosen argument 'est_method' (see above).
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2022}{admix}
#'
#' @return An object of class 'orthobasis_test', containing ten attributes: 1) the number of populations under study (2 in this case);
#'         2) the sizes of samples; 3) the information about the known component distribution; 4) the reject decision of the test; 5) the
#'         confidence level of the test, 6) the p-value of the test; 7) the value of the test statistic; 8) the variance of the test
#'         statistic at each order in the polynomial orthobasis expansion; 9) the selected rank (order) for the test statistic;
#'         10) a vector of estimates, related to the estimated mixing proportions in the two samples.
#'
#' @examples
#' \donttest{
#' #### Under the null hypothesis H0.
#' mixt1 <- twoComp_mixt(n = 300, weight = 0.77,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean" = 1, "sd" = 1),
#'                                         list("rate" = 0.33)))
#' data1 <- getmixtData(mixt1)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                         knownComp_param = mixt1$comp.param[[2]])
#' mixt2 <- twoComp_mixt(n = 500, weight = 0.62,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 1, "sd" = 1),
#'                                         list("mean" = -2, "sd" = 0.5)))
#' data2 <- getmixtData(mixt2)
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                         knownComp_param = mixt2$comp.param[[2]])
#' ## Test procedure:
#' orthobasis_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'                 conf_level = 0.95, est_method = 'BVdk', support = 'Real')
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

orthobasis_test <- function(samples, admixMod, conf_level = 0.95, est_method = c("BVdk","PS"),
                            ask_poly_param = FALSE, K = 3, s = 0.49, nb_echBoot = 100,
                            support = c("Real","Integer","Positive","Bounded.continuous","Bounded.discrete"),
                            bounds_supp = NULL, ...)
{
  old_options_warn <- base::options()$warn
  base::options(warn = -1)

  support <- match.arg(support)
  meth <- match.arg(est_method)

  if ((meth == "PS") & (nb_echBoot <= 1)) stop("Patra & Sen estimator is not square-root n consistent: bootstrap
                                               is necessary to assess the variance of the statistic for the test
                                               to be performed. Please specify a number of bootstrap samples > 1.")

  if (ask_poly_param) {
    K.user <- as.numeric(base::readline("Please enter 'K' (integer), the order for the polynomial expansion in the orthonormal basis: "))
    s.user <- as.numeric(base::readline("Please enter 's' in ]0,0.5[, involved in the penalization rule for model selection where lower values of 's' lead to more powerful tests: "))
  } else {
    K.user <- K
    s.user <- s
  }

  ## Extract useful information about known component distribution:
  knownComp.sim <- paste0("r", sapply(admixMod, "[[", "comp.dist")["known", ])
  comp.sim <- sapply(X = knownComp.sim, FUN = get, mode = "function")
  expr.sim <- vector(mode = "character", length = length(comp.sim))
  sim.knownComp <- vector(mode = "list", length = length(comp.sim))
  for (i in 1:length(comp.sim)) {
    assign(x = names(comp.sim)[i], value = comp.sim[[i]])
    expr.sim[i] <- paste(names(comp.sim)[i],"(n=100000,", paste(names(admixMod[[i]]$comp.param$known),
                                    "=", admixMod[[i]]$comp.param$known, sep = "", collapse = ","), ")", sep="")
    sim.knownComp[[i]] <- eval(parse(text = expr.sim[i]))
  }

  ## Compute expansion coefficients in the orthonormal basis of the known components of the two admixture models:
  known.coef <- list(g1 = NULL, g2 = NULL)
  coef.g1 <- orthoBasis_coef(data = sim.knownComp[[1]], supp = support, degree = K.user, m = 3, other = NULL)
  known.coef$g1 <- sapply(coef.g1, mean)
  coef.g2 <- orthoBasis_coef(data = sim.knownComp[[2]], supp = support, degree = K.user, m = 3, other = NULL)
  known.coef$g2 <- sapply(coef.g2, mean)

  ##---- Splitting the original data for future uncorrelated estimations ----##
  ## Perform this step twice: once for the estimation of component weights, and the other for polynomial coefficients.
  n1 <- length(samples[[1]])
  n2 <- length(samples[[2]])
  indices1 <- sample.int(n = n1, size = floor(n1 / 2), replace = FALSE, prob = NULL)
  data.coef1 <- samples[[1]][indices1]
  data.p1 <- samples[[1]][-indices1]
  n.p1 <- length(data.p1)
  indices2 <- sample.int(n = n2, size = floor(n2 / 2), replace = FALSE, prob = NULL)
  data.coef2 <- samples[[2]][indices2]
  data.p2 <- samples[[2]][-indices2]
  n.p2 <- length(data.p2)

  ##---- Estimate the expansion coefficients of the admixture sample in the orthonormal polynomial basis ----##
  coef.h1 <- coef.h2 <- moy.coef1 <- moy.coef2 <- var.coef1 <- var.coef2 <- NULL
  coef.h1 <- orthoBasis_coef(data = data.coef1, supp = support, degree = K.user, m = 3, other = bounds_supp)
  moy.coef1 <- sapply(coef.h1, mean)
  var.coef1 <- sapply(coef.h1, stats::var)
  coef.h2 <- orthoBasis_coef(data = data.coef2, supp = support, degree = K.user, m = 3, other = bounds_supp)
  moy.coef2 <- sapply(coef.h2, mean)
  var.coef2 <- sapply(coef.h2, stats::var)

  ##---- Estimate nonparametrically the mixture weights (and the estimators variance) ----##
  if (meth == "PS") {
    ## Despite Patra & Sen estimator is not square-root-n consistent, it can be useful when working with
    ## asymetric unknown distributions:
    PS_est1 <- estim_PS(samples = samples[[1]], admixMod = admixMod[[1]], ...)
    PS_est2 <- estim_PS(samples = samples[[2]], admixMod = admixMod[[2]], ...)
    hat.p1 <- getmixingWeight(PS_est1)
    hat.p2 <- getmixingWeight(PS_est2)
  } else {
    BVdk_est1 <- estim_BVdk(samples = data.p1, admixMod = admixMod[[1]], ...)
    BVdk_est2 <- estim_BVdk(samples = data.p2, admixMod = admixMod[[2]], ...)
    hat.p1 <- getmixingWeight(BVdk_est1)
    hat.p2 <- getmixingWeight(BVdk_est2)
    ## Estimation of the variances of the estimators :
    varCov.p1 <- BVdk_varCov_estimators(estim = BVdk_est1, data = data.p1, admixMod = admixMod[[1]])
    var_hat.p1 <- varCov.p1$var.estim_prop
    varCov.p2 <- BVdk_varCov_estimators(estim = BVdk_est2, data = data.p2, admixMod = admixMod[[2]])
    var_hat.p2 <- varCov.p2$var.estim_prop
  }

  ##-------- Computation of the test statistic U -----------##
  ## We get the vector U representing the test statistic for each development order, composed of terms R_kn (cf p5) :
  statU <- (hat.p2 * (moy.coef1 - (1-hat.p1) * known.coef$g1)) - (hat.p1 * (moy.coef2 - (1-hat.p2) * known.coef$g2))

  ##-------- Estimation of variances of the test statistic -----------##
  var.T <- matrix(data = NA, nrow = K.user, ncol = K.user)

  ## With Patra & Sen estimator, we do not have the explicit expression of the variance of the estimators.
  ## With Bordes & Vandekerkhove estimator, we already computed them above.
  if (meth == "PS") {
    message("Testing using Patra & Sen estimator is implemented,
    but is likely to lead to wrong conclusions since the estimator
    variance remains unknown (and so the variance of the test
    statistic). The variance of the test statistic is thus recovered
    using bootstrap, which do not guarantee to provide consistent results.")
    vect.p1 <- vect.p2 <- NULL
    ## Create bootstrap samples: differentiates the cases to estimate weights and the one to estimate the polynomial coefficients.
    indices.bootstrap1 <- t(replicate(n = nb_echBoot, sample.int(n.p1, size = n.p1, replace = TRUE), simplify = "array"))
    indices.bootstrap2 <- t(replicate(n = nb_echBoot, sample.int(n.p2, size = n.p2, replace = TRUE), simplify = "array"))
    bootstrap.samples.p1 <- matrix(data.p1[indices.bootstrap1], nrow = nb_echBoot, ncol = n.p1)
    bootstrap.samples.p2 <- matrix(data.p2[indices.bootstrap2], nrow = nb_echBoot, ncol = n.p2)
    bootstrap.samples.coef1 <- matrix(data.coef1[indices.bootstrap1], nrow = nb_echBoot, ncol = n.p1)
    bootstrap.samples.coef2 <- matrix(data.coef2[indices.bootstrap2], nrow = nb_echBoot, ncol = n.p2)

    ##---------- Estimate the unknown component weight on each bootstrap sample -----------##
	  ##--- by Patra-Sen estimator :
    liste.p1 <- apply(bootstrap.samples.p1, 1, estim_PS, admixMod[[1]], ...)
    vect.p1 <- sapply(liste.p1, "[[", "estimated_mixing_weights")
  	##--- Same task for the second sample Y, still by Patra-Sen estimator :
    liste.p2 <- apply(bootstrap.samples.p2, 1, estim_PS, admixMod[[2]], ...)
    vect.p2 <- sapply(liste.p2, "[[", "estimated_mixing_weights")

    ## Calcul de la statistique de test pour chacun des echantillons bootstrap et pour chaque ordre du developpement:
    coef.h1 <- coef.h2 <- moy.coef1 <- moy.coef2 <- NULL
    statU.boot <- moy.coef1 <- moy.coef2 <- matrix(NA, nrow = nb_echBoot, ncol = K.user)
    for (j in 1:nb_echBoot) {
      coef.h1 <- orthoBasis_coef(data = bootstrap.samples.coef1[j, ], supp = support, degree = K.user, m = 3, other = bounds_supp)
      coef.h2 <- orthoBasis_coef(data = bootstrap.samples.coef2[j, ], supp = support, degree = K.user, m = 3, other = bounds_supp)
      moy.coef1[j, ] <- unlist( lapply(coef.h1, mean, na.rm = TRUE) )
      moy.coef2[j, ] <- unlist( lapply(coef.h2, mean, na.rm = TRUE) )
      ## We get the vector U representing the test statistic for each development order, composed of terms R_kn :
      statU.boot[j, ] <- (vect.p2[j] * (moy.coef1[j, ] - (1-vect.p1[j])*known.coef$g1))  - (vect.p1[j] * (moy.coef2[j, ] - (1-vect.p2[j])*known.coef$g2))
      coef.h1 <- coef.h2 <- NULL
    }
    ## Compute the vector of variances of R_k at each development order (to standardize the test statistic) :
    diag(var.T) <- apply(statU.boot, 2, stats::var)	# diag(var.T) <- pmax(apply(statU.boot, 2, var), 1/log(log(n.X)))
    ## Save memory :
    rm(indices.bootstrap1) ; rm(indices.bootstrap2) ; rm(bootstrap.samples.coef1) ; rm(bootstrap.samples.coef2)
    rm(bootstrap.samples.p1) ; rm(bootstrap.samples.p2) ; rm(liste.p1) ; rm(liste.p2) ; rm(statU.boot)

  } else {
    ##---- Explicit computation of variances of diagonal terms (R_k) ----##
    a <- n1 / (n1 + n2)
    diag(var.T) <- 2 * (1-a) * var_hat.p1 * (hat.p2 * known.coef$g1 - (moy.coef2 - (1-hat.p2) * known.coef$g2))^2 +
                    2 * a * var_hat.p2 * (hat.p1 * known.coef$g2 - (moy.coef1 - (1-hat.p1) * known.coef$g1))^2 +
                    2 * (1-a) * var.coef1 * hat.p2^2 + 2 * a * var.coef2 * hat.p1^2
  }

  ##---- Define the test statistic T ----##
  n.tilde <- (n1*n2) / (n1+n2)
  T.stat <- n.tilde * cumsum( statU^2 / diag(var.T) )

  ##---- Model selection: apply the rule to select the right development order ----##
  ## Remind that under H0, K = 1 asymptotically [s(n) between 1/sqrt(n) and 1 (tuning parameter)]
  normalization.rate <- n.tilde^(s.user-1)
  penalty.importance <- 1:K.user
  penalty <- log(n.tilde)
  penalized.T.stat <- normalization.rate * T.stat - penalty.importance * penalty
  ## Selection de l'ordre jusqu'auquel aller:
  indice.opt <- which.max(penalized.T.stat)
  ## Final assessment of the test statistic: select the right statistic value
  stat.test.final <- T.stat[indice.opt]

  ##---- Decision to reject the null hypothesis (or not) ----##
  rej <- FALSE
  if (stat.test.final > stats::qchisq(conf_level,1)) rej <- TRUE
  ## p-value of the test:
  pvalu <- 1 - stats::pchisq(stat.test.final, 1)

  ## Save memory :
  rm(data.coef1) ; rm(data.coef2) ; rm(data.p1) ; rm(data.p2)
  rm(moy.coef1) ; rm(moy.coef2) ; rm(var.coef1) ; rm(var.coef2)

  obj <- list(
    n_populations = 2,
    population_sizes = c(n1,n2),
    admixture_models = admixMod,
    reject_decision = rej,
    confidence_level = conf_level,
    p_value = pvalu,
    test_statistic_value = stat.test.final,
    varCov_matrix = var.T,
    selected_rank = indice.opt,
    estimated_mixing_weights = c(hat.p1,hat.p2)
    )
  obj$call <- match.call()
  class(obj) <- c("orthobasis_test", "admix_test")

  on.exit(base::options(warn = old_options_warn))
  return(obj)
}


#' Print method for objects of class 'orthobasis_test'
#'
#' @param x An object of class 'orthobasis_test'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.orthobasis_test <- function(x, ...)
{
  cat("Call:")
  print(x$call)
  cat("\n")
  cat("Is the null hypothesis (gaussian unknown component distribution) rejected? ",
      ifelse(x$reject_decision, "Yes", "No"), sep="")
  cat("\nTest p-value: ", round(x$p_value,3), "\n", sep="")
}


#' Summary method for objects of class 'orthobasis_test'
#'
#' @param object An object of class 'orthobasis_test'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.orthobasis_test <- function(object, ...)
{
  cat("Call:")
  print(object$call)
  cat("\n--------- About samples ---------\n")
  cat(paste("Size of sample ", 1:object$n_populations, ": ", object$population_sizes, sep = ""), sep = "\n")
  cat("\n-------- About contamination (admixture) models -------")
  cat("\n")
  for (k in 1:object$n_populations) {
    cat("-> Distribution and parameters of the known component \n for admixture model #", k, ": ", sep="")
    cat(paste(sapply(object$admixture_models[[k]], "[[", "known")[1:2], collapse = " - "))
    cat("\n")
  }
  cat("------- Test decision -------\n")
  cat("Is the null hypothesis (equality of unknown component distributions) rejected? ",
      ifelse(object$reject_decision, "Yes", "No"), sep="")
  cat("\nConfidence level of the test: ", object$confidence_level, sep="")
  cat("\nTest p-value: ", round(object$p_value,3), sep="")
  cat("\n\n------- Test statistic -------\n")
  cat("Selected rank of the test statistic (following the penalization rule): ", object$selected_rank, sep="")
  cat("\nValue of the test statistic: ", round(object$test_statistic_value,4), "\n", sep="")
  cat("Variance-covariance matrix of the test statistic (at each order of expansion):\n", sep = "")
  print(object$varCov_matrix)
  cat("\n------- Estimates -------\n")
  cat("Estimated mixing proportion (listed in the same order as samples): ",
      paste(round(object$estimated_mixing_weights,3), collapse = " "), "\n", sep = "")
}
