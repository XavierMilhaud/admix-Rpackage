#' Equality test of two unknown component distributions using polynomial expansions
#'
#' Tests the null hypothesis (H0: f1=f2) using the decomposition of unknown component densities of two admixture distributions in
#' an adequate orthonormal polynomial basis. Recall that we have two admixture models with respective probability density
#' functions (pdf) l1 = p1*f1 + (1-p1)*g1 and l2 = p2*f2 + (1-p2)*g2, where g1 and g2 are the only known elements and l1 and l2
#' are observed. The admixture weights p1 and p2 thus have to be estimated. For further information on this method, see 'Details' below.
#'
#' @param samples A list of the two observed samples, where each sample follows the mixture distribution given by l = p*f + (1-p)*g,
#'                with f and p unknown and g known.
#' @param known.p (default to NULL) Numeric vector with two elements, respectively the component weight for the unknown component
#'                 in the first and in the second samples.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param known.coef Coefficients in the polynomial basis expansion, corresponding to the known component densities g1 and g2.
#' @param K Number of coefficients considered for the polynomial basis expansion.
#' @param nb.ssEch Number of subsamples created from the original data to decorrelate the estimation of the different parameters.
#' @param s Rate at which the normalization factor is set in the penalization rule for model selection (in ]0,1/2[), see 'Details'.
#' @param var.explicit Boolean that allows to choose between explicit assessment of the variance of the test statistic or not (FALSE=bootstrap),
#'                     FIXME : it seems that bootstrap procedure does not work in the context of admixtures.
#' @param nb.echBoot number of bootstrap samples if 'var.explicit' is set to FALSE.
#' @param support support of the densities under consideration, useful to choose the polynomial orthonormal basis.
#' @param bounds.supp (default to NULL) useful if support = 'bounded', a list of minimum and maximum bounds, specified as
#'                     following: list( list(min.f1,min.g1,min.f2,min.g2) , list(max.f1,max.g1,max.f2,max.g2) )
#' @param est.method Estimation method to get the component weights, either 'PS' (Patra and Sen estimation) or 'BVdk' (Bordes and Vendekerkhove estimation).
#' @param uniformized.knownComp_data (default to NULL) Only useful if 'est.method' has been set to 'PS', and for real-life applications
#'                                   where the distribution of the known component of the admixture model is also unknown. In this
#'                                   case, this known component is previously made uniformly(0,1)-distributed by applying the empirical
#'                                   cumulative distribution of the known component function on the data. This means that
#'                                   all 'comp.dist' and 'comp.param' must be set to NULL.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2022}{admix}
#'
#' @return A list with six elements containing: 1) the rejection decision; 2) the p-value of the test; 3) the test statistic; 4) the
#'         variance-covariance matrix of the test statistic; 5) selected rank for testing, and 6) estimates of the two component weights.
#'
#' @examples
#' \donttest{
#' ###### Using Bordes and Vandekerkhove estimation (valid if symmetric unknown component densities).
#' #### Under the null hypothesis H0.
#' ## Simulate data:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm")
#' list.param <- list(f1 = c(mean = 1, sd = 1), g1 = c(mean = 4, sd = 1),
#'                    f2 = c(mean = 1, sd = 1), g2 = c(mean = 5, sd = 0.5))
#' sim.X <- rsimmix(n = 250, unknownComp_weight=0.9, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' sim.Y <- rsimmix(n = 300, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' plot.2comp_mixt(samples = list(sim.X, sim.Y), support = "continuous")
#' ## Perform the hypothesis test in real-life conditions:
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm")
#' list.param <- list(f1 = NULL, g1 = c(mean = 4, sd = 1),
#'                    f2 = NULL, g2 = c(mean = 5, sd = 0.5))
#' test <- orthoBasis_test_H0(samples = list(sim.X, sim.Y),
#'              known.p=NULL, comp.dist = list.comp, comp.param = list.param, known.coef=NULL, K=3,
#'              nb.ssEch = 2, s = 0.25, var.explicit=TRUE, nb.echBoot=NULL, support = 'Real',
#'              bounds.supp = NULL, est.method = 'BVdk', uniformized.knownComp_data = NULL)
#' test$rejection_rule
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

orthoBasis_test_H0 <- function(samples, known.p = NULL, comp.dist = NULL, comp.param = NULL, known.coef = NULL,
                                K = 3, nb.ssEch = 2, s = 0.49, var.explicit = TRUE, nb.echBoot = NULL,
                                support = c("Real","Integer","Positive","Bounded.continuous","Bounded.discrete"), bounds.supp = NULL,
                                est.method = c("BVdk","PS"), uniformized.knownComp_data = NULL)
{
  meth <- match.arg(est.method)
  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )

  if (all(sapply(comp.dist, is.null)) & all(sapply(comp.param, is.null))) {
    stopifnot(est.method == 'PS')
  } else {
    if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
      stop("Known components must be specified in the two admixture models.")
    }
    if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
      comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
      comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
      if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
    }
  }

  ## Compute expansion coefficients in the orthonormal basis of the known components of the two admixture models:
  if (is.null(known.coef)) {
    known.coef <- list(g1 = NULL, g2 = NULL)
    coef.g1 <- orthoBasis_coef(data = NULL, comp.dist = list(comp.dist[[1]],comp.dist[[2]]), comp.param = list(comp.param[[1]],comp.param[[2]]),
                               supp = support, degree = K, m = 3, other = NULL)
    known.coef$g1 <- sapply(coef.g1, mean)
    coef.g2 <- orthoBasis_coef(data = NULL, comp.dist = list(comp.dist[[3]],comp.dist[[4]]), comp.param = list(comp.param[[3]],comp.param[[4]]),
                               supp = support, degree = K, m = 3, other = NULL)
    known.coef$g2 <- sapply(coef.g2, mean)
  }


  ##---- Splitting the original data for future uncorrelated estimations ----##
  ## Perform this step twice: once for the estimation of component weights, and the other for polynomial coefficients.
  n1 <- length(samples[[1]])
  n2 <- length(samples[[2]])
  indices1 <- sample.int(n = n1, size = floor(n1 / nb.ssEch), replace = FALSE, prob = NULL)
  data.coef1 <- samples[[1]][indices1]
  data.p1 <- samples[[1]][-indices1]
  n.p1 <- length(data.p1)
  indices2 <- sample.int(n = n2, size = floor(n2 / nb.ssEch), replace = FALSE, prob = NULL)
  data.coef2 <- samples[[2]][indices2]
  data.p2 <- samples[[2]][-indices2]
  n.p2 <- length(data.p2)

  ##---- Estimate the expansion coefficients of the admixture sample in the orthonormal polynomial basis ----##
  coef.h1 <- coef.h2 <- moy.coef1 <- moy.coef2 <- var.coef1 <- var.coef2 <- NULL
  coef.h1 <- orthoBasis_coef(data = data.coef1, supp = support, degree = K, m = 3, other = bounds.supp)
  moy.coef1 <- sapply(coef.h1, mean)
  var.coef1 <- sapply(coef.h1, stats::var)
  coef.h2 <- orthoBasis_coef(data = data.coef2, supp = support, degree = K, m = 3, other = bounds.supp)
  moy.coef2 <- sapply(coef.h2, mean)
  var.coef2 <- sapply(coef.h2, stats::var)

  ##---- Estimate nonparametricly the mixture weights (and the estimators variance) ----##
  if (is.null(known.p)) {

    if (meth == "PS") {
      ## Despite the Patra & Sen estimator is not square-root consistent, it can be useful with asymetric unknown distributions:
      ##--- Data transformation necessary to use Patra & Sen estimator ---##
      if (is.null(uniformized.knownComp_data)) {
        data.p1_transfo <- knownComp_to_uniform(data = data.p1, comp.dist = list(comp.dist[[1]], comp.dist[[2]]), comp.param = list(comp.param[[1]], comp.param[[2]]))
        data.p2_transfo <- knownComp_to_uniform(data = data.p2, comp.dist = list(comp.dist[[3]], comp.dist[[4]]), comp.param = list(comp.param[[3]], comp.param[[4]]))
      } else {
        data.p1_transfo <- uniformized.knownComp_data[[1]]
        data.p2_transfo <- uniformized.knownComp_data[[2]]
      }
      ##--- Estimation of component weights by Patra & Sen ---##
      hat.p1 <- PatraSen_est_mix_model(data.p1_transfo, method = "cv", folds = 10, reps = 5, gridsize = 2000)$alp.hat
      hat.p2 <- PatraSen_est_mix_model(data.p2_transfo, method = "cv", folds = 10, reps = 5, gridsize = 2000)$alp.hat

    } else {

      ## Mandatory to have a symmetric unknown density in the admixture model:
      p1.estim <- BVdk_estimParam(data = data.p1, method = "L-BFGS-B", comp.dist = list(comp.dist[[1]], comp.dist[[2]]),
                                  comp.param = list(comp.param[[1]], comp.param[[2]]))
      p2.estim <- BVdk_estimParam(data = data.p2, method = "L-BFGS-B", comp.dist = list(comp.dist[[3]], comp.dist[[4]]),
                                  comp.param = list(comp.param[[3]], comp.param[[4]]))
      hat.p1 <- p1.estim[1]
      hat.p2 <- p2.estim[1]
      ## Estimation of the variances of the estimators :
      varCov.p1 <- BVdk_varCov_estimators(data = data.p1, loc = p1.estim[2], p = hat.p1,
                                          comp.dist = list(comp.dist$f1, comp.dist$g1), comp.param = list(comp.param$f1, comp.param$g1))
      var_hat.p1 <- varCov.p1[["var_pEstim"]]
      varCov.p2 <- BVdk_varCov_estimators(data = data.p2, loc = p2.estim[2], p = hat.p2,
                                          comp.dist = list(comp.dist$f2, comp.dist$g2), comp.param = list(comp.param$f2, comp.param$g2))
      var_hat.p2 <- varCov.p2[["var_pEstim"]]
    }

  } else {
    ## Already known admixture weights:
    hat.p1 <- known.p[1]
    hat.p2 <- known.p[2]
  }

  ##-------- Computation of the test statistic U -----------##
  ## We get the vector U representing the test statistic for each development order, composed of terms R_kn (cf p5) :
  statU <- (hat.p2 * (moy.coef1 - (1-hat.p1) * known.coef$g1)) - (hat.p1 * (moy.coef2 - (1-hat.p2) * known.coef$g2))

  ##-------- Estimation of variances of the test statistic -----------##
  var.T <- matrix(data = NA, nrow = K, ncol = K)

  ## With Patra & Sen estimator, we do not have the explicit expression of the variance of the estimators.
  ## With Bordes & Vandekerkhove estimator, we already computed them above.
  if (!var.explicit) {
    warning("Hypothesis test using Patra & Sen estimator is implemented but is very likely to lead to wrong conclusions
            since the estimator variance remains unknown, and so the variance of the test statistic. The latter is recovered
            by an inappropriate bootstrap procedure that seems to provide strange results.")
    vect.p1 <- vect.p2 <- NULL
    ## Create bootstrap samples: differentiates the cases to estimate weights and the one to estimate the polynomial coefficients.
    indices.bootstrap1 <- t(replicate(n = nb.echBoot, sample.int(n.p1, size = n.p1, replace = TRUE), simplify = "array"))
    indices.bootstrap2 <- t(replicate(n = nb.echBoot, sample.int(n.p2, size = n.p2, replace = TRUE), simplify = "array"))
    bootstrap.samples.p1 <- matrix(data.p1_transfo[indices.bootstrap1], nrow = nb.echBoot, ncol = n.p1)
    bootstrap.samples.p2 <- matrix(data.p2_transfo[indices.bootstrap2], nrow = nb.echBoot, ncol = n.p2)
    bootstrap.samples.coef1 <- matrix(data.coef1[indices.bootstrap1], nrow = nb.echBoot, ncol = n.p1)
    bootstrap.samples.coef2 <- matrix(data.coef2[indices.bootstrap2], nrow = nb.echBoot, ncol = n.p2)

    ##---------- Estimate the unknown component weight on each bootstrap sample -----------##
	  ##--- by Patra-Sen estimator :
    liste.p1 <- apply(bootstrap.samples.p1, 1, PatraSen_est_mix_model, method = "fixed", gridsize = 1000)
    vect.p1 <- sapply(liste.p1, "[[", "alp.hat")
  	##--- Same task for the second sample Y, still by Patra-Sen estimator :
    liste.p2 <- apply(bootstrap.samples.p2, 1, PatraSen_est_mix_model, method = "fixed", gridsize = 1000)
    vect.p2 <- sapply(liste.p2, "[[", "alp.hat")

    ## Calcul de la statistique de test pour chacun des echantillons bootstrap et pour chaque ordre du developpement:
    coef.h1 <- coef.h2 <- moy.coef1 <- moy.coef2 <- NULL
    statU.boot <- moy.coef1 <- moy.coef2 <- matrix(NA, nrow = nb.echBoot, ncol = K)
    for (j in 1:nb.echBoot) {
      coef.h1 <- orthoBasis_coef(data = bootstrap.samples.coef1[j, ], supp = support, degree = K, m = 3, other = bounds.supp)
      coef.h2 <- orthoBasis_coef(data = bootstrap.samples.coef2[j, ], supp = support, degree = K, m = 3, other = bounds.supp)
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
  normalization.rate <- n.tilde^(s-1)
  penalty.importance <- 1:K
  penalty <- log(n.tilde)
  penalized.T.stat <- normalization.rate * T.stat - penalty.importance * penalty
  ## Selection de l'ordre jusqu'auquel aller:
  indice.opt <- which.max(penalized.T.stat)
  ## Final assessment of the test statistic: select the right statistic value
  stat.test.final <- T.stat[indice.opt]

  ##---- Decision to reject the null hypothesis (or not) ----##
  rej <- FALSE
  if (stat.test.final > stats::qchisq(0.95,1)) rej <- TRUE
  ## P-value of the test:
  pvalu <- 1 - stats::pchisq(stat.test.final, 1)

  ## Save memory :
  rm(data.coef1) ; rm(data.coef2) ; rm(data.p1) ; rm(data.p2)
  rm(moy.coef1) ; rm(moy.coef2) ; rm(var.coef1) ; rm(var.coef2)

  list(rejection_rule = rej, p_value = pvalu, test.stat = stat.test.final, varCov.matrix = var.T,
       rank = indice.opt, p1 = hat.p1, p2 = hat.p2)
}


#' Expansion coefficients for some given orthonormal polynomial basis.
#'
#' Compute the coefficients of the decomposition of some density in a given orthonormal polynomial basis.
#'
#' @param data Observed sample from which the coefficients are calculated. Can be NULL if 'comp.dist' and 'comp.param' are specified.
#' @param comp.dist (default to NULL) A list with two elements corresponding to component distributions (specified with
#'                  R native names for these distributions) involved in the admixture model.
#'                  Unknown elements must be specified as 'NULL' objects (for instance unknown 'f': list(f=NULL, g='norm')).
#' @param comp.param (default to NULL) A list with two elements corresponding to the parameters of the component distributions, each
#'                   element being a list itself. The names used in this list must correspond to the native R argument names
#'                   for these distributions. Unknown elements must be specified as 'NULL' objects.
#'                   For instance if 'f' is unknown: list(f = NULL, g = list(mean=0,sd=1)).
#' @param supp Support of the density considered.
#' @param degree Degree up to which the polynomial basis is built.
#' @param m (default to 3) Only used when support is 'Integer'. Corresponds to the mean of the reference measure, i.e. Poisson(m).
#' @param other (default to NULL) A list to precise bounds when the support is bounded, where the second and fourth elements give bounds.
#'
#' @return The list composed of 'degree' elements, each element being a numeric vector (with sample size) where each value represents
#'         the k-th order coefficient found when decomposing the density in the orthonormal polynomial basis.
#'
#' @examples
#' ## Simulate data:
#' sample1 <- rnorm(n = 7000, mean = 3, sd = 1)
#' ## Compute the expansion coefficients in the orthonormal polynomial basis:
#' coeff <- orthoBasis_coef(data = sample1, comp.dist = NULL, comp.param = NULL, supp = 'Real',
#'                          degree = 3, m = 3, other = NULL)
#' sapply(coeff, mean)
#' ## No observed data and decomposition of the known component of the admixture model:
#' coeff <- orthoBasis_coef(data = NULL, comp.dist = list(NULL, 'norm'),
#'             comp.param=list(NULL,list(mean=3,sd=1)), supp = 'Real', degree=3, m=3, other = NULL)
#' sapply(coeff, mean)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

orthoBasis_coef <- function(data, comp.dist = NULL, comp.param = NULL, supp = c('Real','Integer','Positive','Bounded.continuous'),
                            degree, m = 3, other = NULL)
{
  if (is.null(data)) {
    stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
    if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component of the admixture model must be specified.")
    ## Extracts the information on component distributions and stores in expressions:
    comp.dist.sim <- paste0("r", comp.dist[[2]])
    #    comp_ortho <- sapply(X = comp.dist.sim, FUN = get, pos = "package:stats", mode = "function")
    comp_ortho <- sapply(X = comp.dist.sim, FUN = get, mode = "function")
    assign(x = names(comp_ortho)[1], value = comp_ortho[[1]])
    expr.sim <- paste(names(comp_ortho)[1],"(n=1000000,", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep = "", collapse = ","), ")", sep="")
    data <- eval(parse(text = expr.sim))
  }

  ## Builds the orthonormal polynomial basis:
  poly_basis <- poly_orthonormal_basis(support = supp, deg = degree, x = data, m = m)

  ## Store the coefficients to be computed :
  coef.list <- vector(mode = "list", length = degree)

  ## Depending on the support :
  if (supp == "Real") {
    ## Reference measure N(0,1)
    for (i in 1:degree) coef.list[[i]] <- orthopolynom::polynomial.values(poly_basis, data)[[i+1]] / sqrt(factorial(i))

  } else if (supp == "Integer") {
    ## Reference measure P(3)
    for (i in 1:degree) coef.list[[i]] <- poly_basis[ ,i]

  } else if (supp == "Positive") {
    ## Reference measure Exp(1)
    for (i in 1:degree) coef.list[[i]] <- orthopolynom::polynomial.values(poly_basis, data)[[i+1]]

  } else if (supp == "Bounded.continuous") {
    ## Reference measure Unif(a,b)
    if (is.null(other)) { bounds <- c(min(data), max(data))
    } else { bounds <- other[[2]] }
    for (i in 1:degree) coef.list[[i]] <- (orthopolynom::polynomial.values(poly_basis, (2*data-bounds[1]-bounds[2])/(bounds[2]-bounds[1]))[[i+1]]) / sqrt(2*i+1)

  } else stop("Change the support since the choosen one is not considered!")

  return(coef.poly = coef.list)
}
