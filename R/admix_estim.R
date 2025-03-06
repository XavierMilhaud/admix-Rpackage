#' Estimate the unknown parameters of the admixture model(s)
#'
#' Estimate the component weights, the location shift parameter (in case of a symmetric unknown component density),
#' and the unknown component distribution using different estimation techniques. We remind that the i-th admixture
#' model has probability density function (pdf) l_i such that:
#'    l_i = p_i * f_i + (1-p_i) * g_i, where g_i is the known component density.
#' The unknown quantities p_i and f_i then have to be estimated.
#'
#' @param samples A list of the K (K>0) samples to be studied, all following admixture distributions.
#' @param admixMod A list of objects of class \link[admix]{admix_model}, containing useful information about distributions and parameters.
#' @param est_method The estimation method to be applied. Can be one of 'BVdk' (Bordes and Vandekerkhove estimator), 'PS' (Patra and Sen
#'         estimator), or 'IBM' (Inversion Best-Matching approach) in the continuous case (continuous random variable). Only 'IBM' for
#'         discrete random variables. The same estimation method is performed on each sample if several samples are provided.
#' @param ... Optional arguments to \link[admix]{estim_PS}, \link[admix]{estim_BVdk} or \link[admix]{estim_IBM} depending on the
#'            choice made by the user for the estimation method.
#'
#' @details For further details on the different estimation techniques, see references below on i) Patra and Sen estimator ;
#'          ii) Bordes and Vandekerkhove estimator ; iii) Inversion Best-Matching approach. Important note: estimation by 'IBM'
#'          requires at least two samples at hand, and provides unbiased estimators only if the distributions of unknown components
#'          are equal (meaning that it requires to perform previously this test between the pairs of samples, see ?admix_test).
#'
#' @references
#' \insertRef{PatraSen2016}{admix}
#' \insertRef{BordesDelmasVandekerkhove2006}{admix}
#' \insertRef{BordesVandekerkhove2010}{admix}
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return An object of class \link[admix]{admix_estim}, containing at least 5 attributes: 1) the number of samples under study; 2) the information
#'         about the mixture components (distributions and parameters); 3) the sizes of the samples; 4) the chosen estimation technique
#'         (one of 'BVdk', 'PS' or 'IBM'); 5) the estimated mixing proportions (weights of the unknown component distributions in the
#'         mixture model). In case of 'BVdk' estimation, one additional attribute corresponding to the estimated location shift parameter
#'         is included.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 300, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 250, weight = 0.85,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("rate" = 1)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' # Estimation by different methods:
#' admix_estim(samples=list(data1), admixMod=list(admixMod1), est_method = "BVdk")
#' admix_estim(samples=list(data1,data2), admixMod=list(admixMod1,admixMod2), est_method = "PS")
#' admix_estim(samples=list(data1,data2), admixMod=list(admixMod1,admixMod2), est_method = "IBM")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_estim <- function(samples, admixMod, est_method = c("PS","BVdk","IBM"), ...)
{
  meth <- match.arg(est_method)
  supp <- detect_support_type(unlist(samples))
  if ((supp != "Continuous") & (meth != "IBM"))
    stop("'BVdk' or 'PS' estimation methods are not suitable to discrete random variables.")

  n_samples <- length(samples)
  ## Check right specification of arguments:
  if ((n_samples == 1) & (meth == "IBM")) stop("Estimation by 'IBM' requires (at least) two samples.")

  n_obs <- sapply(X = samples, FUN = length)
  estimate <- vector(mode = "list", length = n_samples)
  if (meth == "BVdk") {
    for (k in 1:n_samples) {
      estimate[[k]] <- estim_BVdk(samples = samples[[k]], admixMod = admixMod[[k]], ...)
    }
    estim_weight <- sapply(X = estimate, "[[", "estimated_mixing_weights")
    estim_loc <- sapply(X = estimate, "[[", "estimated_locations")
  } else if (meth == "PS") {
    for (k in 1:n_samples) {
      estimate[[k]] <- estim_PS(samples = samples[[k]], admixMod = admixMod[[k]], ...)
    }
    estim_weight <- sapply(X = estimate, "[[", "estimated_mixing_weights")
  } else if (meth == "IBM") {
    for (k in 2:n_samples) {
      estimate[[k]] <- estim_IBM(samples = list(samples[[1]], samples[[k]]),
                                 admixMod = list(admixMod[[1]], admixMod[[k]]), ...)
    }
    estimate <- estimate[-sapply(X = estimate, FUN = is.null)]
    estim_weight <- sapply(X = estimate, FUN = "[[", "estimated_mixing_weights")
    equal_knownComps <- sapply(X = estimate, FUN = "[[", "equal.knownComp")
    fixed_prop <- sapply(X = estimate, FUN = "[[", "p.X.fixed")
  } else stop("Please choose appropriately the arguments of the function.")

  estimators <- list(
    n_populations = n_samples,
    admixture_models = admixMod,
    population_sizes = n_obs,
    estimation_method = switch(meth, "BVdk" = "Bordes and Vandekerkhove (BVdk)",
                               "PS" = "Patra and Sen (PS)",
                               "IBM" = "Inversion Best-Matching (IBM)"),
    estimated_mixing_weights = estim_weight
    )
  if (meth == "BVdk") estimators$estimated_locations = estim_loc
  if (meth == "IBM") {
    estimators$equal_knownComp = equal_knownComps
    estimators$fixed_prop = fixed_prop
  }
  class(estimators) <- "admix_estim"
  estimators$call <- match.call()

  return(estimators)
}


#' Print the estimated parameters from K admixture models
#'
#' @param x An object of class 'admix_estim' (see ?admix_estim).
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_estim <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  if (x$estimation_method != "Inversion Best-Matching (IBM)") {
    cat("\n")
    cat(paste("Estimated weight of the unknown distribution in Sample ", 1:x$n_populations,
              ": ", round(x$estimated_mixing_weights,2), sep = ""), sep = "\n")
  } else {
    cat("\nPairwise estimation performed (IBM estimation method). \n")
    cat("N.B.: estimated weights are reliable only if the unknown
component distributions have previously been tested equal.
See ?admix_test. \n\n")
    for (k in 1:(x$n_populations-1)) {
      cat("----- Pair 1 and ", k+1, " -----\n", sep = "")
      if (!x$equal_knownComp[k]) {
        if (inherits(x$estimated_mixing_weights, "matrix")) {
          cat("Estimated weight of the unknown distribution in Sample 1: ",
              round(x$estimated_mixing_weights[ ,k][1],2), "\n", sep = "")
          cat("Estimated weight of the unknown distribution in Sample ", k+1, ": ",
              round(x$estimated_mixing_weights[ ,k][2],2), sep = "")
          cat("\n")
        } else {
          cat("Estimated weight of the unknown distribution in Sample 1: ",
              round(x$estimated_mixing_weights[[k]][1],2), "\n", sep = "")
          cat("Estimated weight of the unknown distribution in Sample ", k+1, ": ",
              round(x$estimated_mixing_weights[[k]][2],2), sep = "")
          cat("\n")
        }
      } else {
        cat("Estimated weights of the unknown distributions in Sample 1 and ", k+1, ": ",
            paste(c(x$fixed_prop[k],round(x$estimated_mixing_weights[[k]],2)), collapse=", ", sep = ""), sep="")
        cat("\nFirst proportion was fixed because of equal known components.
The estimated proportion for the 2nd sample makes sense through the ratio of
the two proportions, which should roughly be similar than the true actual ratio.\n")
        cat("\n")
      }
    }
  }
  #cat("\n")
  if (x$estimation_method == "Bordes and Vandekerkhove (BVdk)") {
    cat(paste("Estimated location of the unknown distribution in Sample ", 1:x$n_populations,
              ": ", round(x$estimated_locations,2), sep = ""), sep = "\n")
  }
}


#' Results of estimated parameters from K admixture models
#'
#' Summarize the estimated weight(s) of the unknown component(s) in the admixture model(s) under study.
#' Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with G a known CDF, and p and F unknown quantities.
#'
#' @param object An object of class 'admix_estim' (see ?admix_estim).
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.admix_estim <- function(object, ...)
{
  cat("Call:\n")
  print(object$call)
  cat("\n------- About samples ------\n")
  cat(paste("Size of sample ", 1:object$n_populations, ": ", object$population_sizes, sep = ""), sep = "\n")
  cat("\n------ About contamination (admixture) models -----")
  cat("\n")
  if (object$n_populations == 1) {
    cat("-> Distribution and parameters of the known component \n for the admixture model: ", sep="")
    cat(object$admixture_models[[1]]$comp.dist$known, "\n")
    print(unlist(object$admixture_models[[1]]$comp.param$known, use.names = TRUE))
  } else {
    for (k in 1:object$n_populations) {
      cat("-> Distribution and parameters of the known component \n for admixture model #", k, ": ", sep="")
      cat(paste(sapply(object$admixture_models[[k]], "[[", "known")[1:2], collapse = " - "))
      cat("\n")
    }
  }
  cat("\n------ About estimation -----")
  cat("\nEstimation approach: ", object$estimation_method,  "\n", sep = "")
  if (object$estimation_method != "Inversion Best-Matching (IBM)") {
    #cat("\n")
    cat(paste("Estimated weight of the unknown distribution in Sample ", 1:object$n_populations,
              ": ", round(object$estimated_mixing_weights,2), sep = ""), sep = "\n")
  } else {
    cat("\nPairwise estimation performed (IBM estimation method). \n")
    cat("N.B.: estimated weights are reliable only if the unknown
component distributions have previously been tested equal.
See ?admix_test. \n\n")
    for (k in 1:(object$n_populations-1)) {
      cat("--- Pair 1 and ", k+1, " ---\n", sep = "")
      if (!object$equal_knownComp[k]) {
        if (inherits(object$estimated_mixing_weights, "matrix")) {
          cat("Estimated weight of the unknown distribution in Sample 1: ",
              round(object$estimated_mixing_weights[ ,k][1],2), "\n", sep = "")
          cat("Estimated weight of the unknown distribution in Sample ", k+1, ": ",
              round(object$estimated_mixing_weights[ ,k][2],2), sep = "")
          cat("\n")
        } else {
          cat("Estimated weight of the unknown distribution in Sample 1: ",
              round(object$estimated_mixing_weights[[k]][1],2), "\n", sep = "")
          cat("Estimated weight of the unknown distribution in Sample ", k+1, ": ",
              round(object$estimated_mixing_weights[[k]][2],2), sep = "")
          cat("\n")
        }
      } else {
        cat("Estimated weights of the unknown distributions in Sample 1 and ", k+1, ": ",
            paste(c(object$fixed_prop[k],round(object$estimated_mixing_weights[[k]],2)), collapse=", ", sep = ""), sep="")
        cat("\nFirst proportion was fixed because of equal known components.
The estimated proportion for the 2nd sample makes sense through the ratio of
the two proportions, which should roughly be similar than the true actual ratio.\n")
        cat("\n")
      }
    }
  }
  cat("\n")
  if (object$estimation_method == "Bordes and Vandekerkhove") {
    cat("-> Location parameters (symmetric density assumption):\n")
    cat(paste("Estimated location of the unknown distribution in Sample ", 1:object$n_populations,
              ": ", round(object$estimated_locations,2), sep = ""), sep = "\n")
  }
}



#' Extractor for object of class 'admix_estim'
#'
#' Get the estimated unknown mixing proportion related to the weight of
#' the unknown component distribution of the admixture model.
#'
#' @param x An object of class 'admix_estim'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

getmixingWeight <- function(x)
{
  if (!inherits(x, "admix_estim") & !inherits(x, "admix_test") )
    stop("This function must be used with objects of class 'admix_estim' or 'admix_test'")
  x$estimated_mixing_weights
}
