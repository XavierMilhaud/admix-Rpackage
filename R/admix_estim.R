#' Estimate the unknown parameters of the admixture model(s)
#'
#' Estimate the component weights, the location shift parameter (in case of a symmetric unknown component density),
#' and the unknown component distribution using different estimation techniques. We remind that the i-th admixture
#' model has probability density function (pdf) l_i such that:
#'    l_i = p_i * f_i + (1-p_i) * g_i, where g_i is the known component density.
#' The unknown quantities p_i and f_i then have to be estimated.
#'
#' @param samples (list) A list of the K (K>0) samples to be studied, all following admixture distributions.
#' @param admixMod (list) A list of objects of class 'admix_model', containing useful information about distributions and parameters.
#' @param est.method The estimation method to be applied. Can be one of 'BVdk' (Bordes and Vandekerkhove estimator), 'PS' (Patra and Sen
#'         estimator), or 'IBM' (Inversion Best-Matching approach). The same estimation method is performed on each sample.
#'         Important note: estimation by 'IBM' is unbiased only under H0, meaning that choosing this method requires to perform
#'         previously the test hypothesis between the pairs of samples. For further details, see section 'Details' below.
#' @param sym.f A boolean indicating whether the unknown component densities are assumed to be symmetric or not.
#'
#' @details For further details on the different estimation techniques, see references below i) Patra and Sen estimator ;
#'          ii) BVdk estimator ; iii) IBM approach.
#'
#' @references
#' \insertRef{PatraSen2016}{admix}
#' \insertRef{BordesDelmasVandekerkhove2006}{admix}
#' \insertRef{BordesVandekerkhove2010}{admix}
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return An object of class 'admix_estim', containing at least 5 attributes: 1) the number of samples under study; 2) the information
#'         about the mixture components (distributions and parameters); 3) the sizes of the samples; 4) the chosen estimation technique
#'         (one of 'BVdk', 'PS' or 'IBM'); 5) the estimated mixing proportions (weights of the unknown component distributions in the
#'         mixture model). In case of 'BVdk' estimation, one additional attribute corresponding to the estimated location shift parameter
#'         is included.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 250, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 200, weight = 0.85,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean" = 3, "sd" = 1),
#'                                         list("rate" = 1)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admix_estim(samples = list(data1), admixMod = list(admixMod1),
#'             est.method = 'BVdk', sym.f = TRUE)
#' admix_estim(samples = list(data1,data2),
#'             admixMod = list(admixMod1,admixMod2),
#'             est.method = 'PS')
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_estim <- function(samples, admixMod, est.method = c("BVdk","PS","IBM"), sym.f = FALSE)
{
  meth <- match.arg(est.method)
  supp <- detect_support_type(unlist(samples))
  if ((supp != "Continuous") & (meth != "IBM"))
    stop("'BVdk' or 'PS' estimation methods are not suitable to discrete random variables.")

  n_samples <- length(samples)
  ## Check right specification of arguments:
  if ((n_samples == 1) & (meth == "IBM")) stop("Estimation by IBM requires (at least) two samples.")
  if ((sym.f == FALSE) & (meth == "BVdk")) stop("Estimation by BVdk assumes the unknown component density
                                                to be symmetric, thus set argument 'sym.f' to TRUE.")

  n_obs <- sapply(X = samples, FUN = length)
  estimate <- vector(mode = "list", length = n_samples)
  if (meth == "BVdk") {
    for (k in 1:n_samples) {
      estimate[[k]] <- estim_BVdk(data = samples[[k]], admixMod = admixMod[[k]], method = 'L-BFGS-B')
    }
    estim_weight <- sapply(X = estimate, "[[", "estimated_mixing_weights")
    estim_loc <- sapply(X = estimate, "[[", "estimated_locations")
  } else if (meth == "PS") {
    for (k in 1:n_samples) {
      estimate[[k]] <- estim_PS(data = samples[[k]], admixMod = admixMod[[k]],
                                method = 'fixed', c.n = 0.1*log(log(n_obs[k])), gridsize = 2000)
    }
    estim_weight <- sapply(X = estimate, "[[", "estimated_mixing_weights")
  } else if (meth == "IBM") {
    for (k in 2:n_samples) {
      estimate[[k]] <- estim_IBM(samples = list(samples[[1]],samples[[k]]), admixMod = list(admixMod[[1]],admixMod[[k]]), n.integ = 1000)
    }
    estimate <- estimate[-sapply(X = estimate, FUN = is.null)]
    estim_weight <- sapply(X = estimate, FUN = "[[", "estimated_mixing_weights")
    colnames(estim_weight) <- paste("Samples 1 and ", 2:n_samples, "/", sep="")
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
  if (x$estimation_method != "IBM") {
    cat("\n")
    cat(paste("Estimated mixing weight of the unknown component distribution in Sample ", 1:x$n_populations,
              ": ", round(x$estimated_mixing_weights,2), sep = ""), sep = "\n")
  } else {
    cat("\nPairwise estimation performed (between samples listed as column names): \n")
    print(round(x$estimated_mixing_weights,2))
    cat("First row provides estimated mixing weight in the first of the two samples listed,
    and second row the one in the second sample.\n")
  }
  cat("\n")
  if (x$estimation_method == "BVdk") {
    cat(paste("Estimated location parameters of the unknown component distribution in Sample ", 1:x$n_populations,
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
  cat("\n--------- About samples ---------\n")
  cat(paste("Size of sample ", 1:object$n_populations, ": ", object$population_sizes, sep = ""), sep = "\n")
  cat("\n-------- About contamination (admixture) models -------")
  cat("\n")
  if (object$n_populations == 1) {
    cat("-> Distribution and parameters of the known component \n for the admixture model: ", sep="")
    cat(object$admixture_models$comp.dist$known, "\n")
    print(unlist(object$admixture_models$comp.param$known, use.names = TRUE))
  } else {
    for (k in 1:object$n_populations) {
      cat("-> Distribution and parameters of the known component \n for admixture model #", k, ": ", sep="")
      cat(paste(sapply(object$admixture_models[[k]], "[[", "known")[1:2], collapse = " - "))
      cat("\n")
    }
  }
  cat("\n-------- About estimation -------")
  cat("\nEstimation approach: ", object$estimation_method,  "\n\n", sep = "")
  cat("-> Mixing weights:\n")
  cat(paste("Estimated proportion of the unknown component in Sample ", 1:object$n_populations,
            ": ", round(object$estimated_mixing_weights,2), sep = ""), sep = "\n")
  cat("\n")
  if (object$estimation_method == "Bordes and Vandekerkhove") {
    cat("-> Location parameters (symmetric density assumption):\n")
    cat(paste("Estimated location parameters of the unknown distribution in Sample ", 1:object$n_populations,
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
  if (!inherits(x, "admix_estim")) stop("This function must be used with objects of class 'admix_estim'")
  x$estimated_mixing_weights
}
