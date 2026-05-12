#' Estimate the unknown weight in an admixture model
#'
#' Estimate the unknown component weight (and possibly a location shift parameter
#' in case of a symmetric unknown component density), using different estimation
#' techniques. We recall that the \eqn{i}-th admixture model has probability density function
#' \eqn{\ell_i} such that:
#' \deqn{
#'   \ell_i = p_i f_i + (1 - p_i) g_i,
#' }
#' where \eqn{g_i} is the known component density. The unknown quantities
#' \eqn{p_i} and \eqn{f_i} then have to be estimated.
#'
#' @param samples A list of the K (K>0) samples to be studied, all following admixture distributions.
#' @param admixMod A list of objects of class \link[admix]{admix_model}, with information about known distributions and known parameters.
#' @param est_method The estimation method to be applied. Can be one of 'BVdk' (Bordes and Vandekerkhove estimator), 'PS' (Patra and Sen
#'         estimator), or 'IBM' (Inversion Best-Matching approach) in the continuous case (continuous random variable). Only 'IBM' for
#'         discrete random variables. The same estimation method is performed on each sample if several samples are provided.
#' @param ... Optional arguments to \link[admix]{estim_PS}, \link[admix]{estim_BVdk} or \link[admix]{estim_IBM} depending on the
#'            choice made by the user for the estimation method.
#'
#' @details For further details on the different estimation techniques, see references below on i) Patra and Sen estimator ;
#'          ii) Bordes and Vandekerkhove estimator ; iii) Inversion Best-Matching approach. Important note: estimation by 'IBM'
#'          requires at least two samples at hand, and provides unbiased estimators only if the distributions of unknown components
#'          are equal (meaning that it requires to perform previously this test between the pairs of samples, see \link[admix]{admix_test}).
#'
#' @return An object of class \code{estim_BVdk}, \code{estim_PS} or \code{estim_IBM} (that inherits from class \link[admix]{admix_estim}),
#'         with two attributes, 'class' and 'names'. The latter contains three elements, among which 'estim_objects' that lists for each
#'         sample under study all the information of the estimation procedure.
#'
#' @seealso [get_mixing_weights()] to access the estimated mixing weight(s), [get_known_component()] to access the known component(s), [print.admix_estim()] for a brief description of the results,
#'          and [summary.admix_estim()] for an overview of the estimation process. More precisely, 1) the number of samples under study;
#'          2) the information about the known mixture components (distributions and parameters); 3) the sizes of the samples;
#'          4) the chosen estimation technique (one of 'BVdk', 'PS' or 'IBM'); 5) the estimated mixing proportions (weights of the
#'          unknown component distributions in the mixture model). In case of 'BVdk' estimation, one additional attribute corresponding
#'          to the estimated location shift parameter is included.
#'
#' @references
#' \insertRef{PatraSen2016}{admix}
#' \insertRef{BordesDelmasVandekerkhove2006}{admix}
#' \insertRef{BordesVandekerkhove2010}{admix}
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
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
#' mixt3 <- twoComp_mixt(n = 500, weight = 0.5,
#'                       comp.dist = list("pois", "pois"),
#'                       comp.param = list(list("lambda" = 2),
#'                                         list("lambda" = 7)))
#' mixt4 <- twoComp_mixt(n = 1500, weight = 0.2, comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size"=1, "prob" = c(0.8,0.1,0.1)),
#'                                    list("size"=1, "prob" = c(0.1,0.2,0.7))))
#' data1 <- get_mixture_data(mixt1)
#' data2 <- get_mixture_data(mixt2)
#' data3 <- get_mixture_data(mixt3)
#' data4 <- get_mixture_data(mixt4)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
#'                          knownComp_param = mixt3$comp.param[[2]])
#' admixMod4 <- admix_model(knownComp_dist = mixt4$comp.dist[[2]],
#'                          knownComp_param = mixt4$comp.param[[2]])
#' # Estimation by different methods:
#' admix_estim(samples = list(data1), admixMod = list(admixMod1), est_method = "BVdk")
#' admix_estim(samples = list(data1, data2, data3, data4),
#'             admixMod = list(admixMod1, admixMod2, admixMod3, admixMod4), est_method = "PS")
#' admix_estim(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2), est_method = "IBM")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_estim <- function(samples, admixMod, est_method = c("PS","BVdk","IBM"), ...)
{
  if (!is.list(samples) | !is.list(admixMod))
    stop("Please provide sample(s) AND admixture model(s) in a list, also with only one sample!")
  if (!all(sapply(X = admixMod, FUN = inherits, what = "admix_model")))
    stop("Argument 'admixMod' is not correctly specified. See ?admix_model.")

  meth <- match.arg(est_method)
  supp <- detect_support_type(unlist(samples))
  if ((supp != "Continuous") & (meth == "BVdk"))
    stop("'BVdk' estimation method is not suitable to discrete random variables.")

  n_samples <- length(samples)
  ## Check right specification of arguments:
  if ((n_samples == 1) & (meth == "IBM")) stop("Estimation by 'IBM' requires (at least) two samples.")

  n_obs <- sapply(X = samples, FUN = length)
  estimate <- vector(mode = "list", length = n_samples)
  if (meth == "BVdk") {
    message("Mixing weight estimation using 'BVdk' assumes the unknown component
distribution to have a symmetric probability density function.")
    for (k in 1:n_samples) {
      estimate[[k]] <- estim_BVdk(samples = samples[[k]], admixMod = admixMod[[k]], ...)
    }
  } else if (meth == "PS") {
    for (k in 1:n_samples) {
      estimate[[k]] <- estim_PS(samples = samples[[k]], admixMod = admixMod[[k]], ...)
      #estimate[[k]]$data.name <- sample_names[k]
    }
  } else if (meth == "IBM") {
    message(" IBM estimators of two unknown proportions are reliable only if the two corresponding
 unknown component distributions have previously been tested equal (see ?admix_test).")
    any_knownComp_equal <- vector(mode = "logical", length = (n_samples-1))
    for (k in 2:n_samples) { any_knownComp_equal[k-1] <- is_equal_knownComp(admixMod[[1]], admixMod[[k]]) }
    if (any(any_knownComp_equal == TRUE)) {
      message("/n When both the known and unknown component distributions of the mixture models are
 identical, IBM provides an estimated ratio of the mixing weights (and not the weights).\n")
    }
    for (k in 2:n_samples) {
      estimate[[k]] <- estim_IBM(samples = list(samples[[1]], samples[[k]]),
                                 admixMod = list(admixMod[[1]], admixMod[[k]]), ...)
    }
    estimate[[1]] <- NULL
  } else stop("Please choose appropriately the arguments of the function.")

  estimators <- list(estim_objects = estimate)
  specific_class <- switch(meth, "BVdk" = "estim_BVdk",
                           "PS" = "estim_PS",
                           "IBM" = "estim_IBM")
  class(estimators) <- c(specific_class, "admix_estim")
  estimators$call <- match.call()
  ## Retrieve names of objects
  sample_names <- NULL
  sample_expr <- match.call()$samples
  if (is.call(sample_expr) && sample_expr[[1]] == as.name("list")) { sample_names <- as.character(sample_expr)[-1] }
  ## fallback if not retrievable
  if (is.null(sample_names) || length(sample_names) != n_samples) { sample_names <- paste0("Sample_", seq_len(n_samples)) }
  estimators$sample_names <- sample_names

  return(estimators)
}


#' Print method for object of class \code{admix_estim}
#'
#' @param x An object of class \code{admix_estim} (see ?admix_estim).
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_estim <- function(x, ...) {

  cat("\nCall:\n")
  print(x$call)
  method <- class(x)[1]
  method <- sub("estim_", "", method)
  cat("\nMethod:", method, " - ")

  if (inherits(x, "estim_IBM")) {
    valid_objects <- Filter(Negate(is.null), x$estim_objects)
    n_samples     <- length(valid_objects) + 1
    cat(" Pairwise estimation\n\n")
    sample_names <- x$sample_names
    if (is.null(sample_names)) { sample_names <- names(valid_objects) }
    if (is.null(sample_names) || any(sample_names == "")) { sample_names <- paste0("Sample_", seq_len(n_samples)) }
    pairs <- paste0(sample_names[1], " vs ", sample_names[-1])

    rows <- lapply(seq_along(valid_objects), function(k) {
      obj <- valid_objects[[k]]
      w   <- obj$estimated_mixing_weights
      variance.p1 <- obj$variance_est_p1
      variance.p2 <- obj$variance_est_p2
      if (isTRUE(obj$equal.knownComp)) {
        data.frame(pair = pairs[k], size_1st = obj$population_sizes[1], size_2nd = obj$population_sizes[2],
                   `mix_weight_1st (fixed)` = format(round(obj$p.X.fixed, 3), nsmall=3), var_1st = format(round(variance.p1,5), nsmall=5),
                   mix_weight_2nd = format(round(w,3), nsmall=3), var_2nd = format(round(variance.p2,5), nsmall=5), check.names = FALSE)
      } else {
        data.frame(pair = pairs[k], size_1st = obj$population_sizes[1], size_2nd = obj$population_sizes[2],
                   `mix_weight_1st` = format(round(w[1], 3), nsmall = 3), var_1st = format(round(variance.p1,5), nsmall=5),
                   mix_weight_2nd = format(round(w[2], 3), nsmall = 3), var_2nd = format(round(variance.p2,5), nsmall=5), check.names = FALSE)
      }
    })
    ## Equivalent of dplyr::bind_rows (fills missing columns with NA)
    has_fixed <- any(sapply(valid_objects, function(obj) isTRUE(obj$equal.knownComp)))
    all_cols <- if (has_fixed) {
      c("pair", "size_1st", "size_2nd", "mix_weight_1st (fixed)", "mix_weight_1st", "var_1st", "mix_weight_2nd", "var_2nd")
    } else {
      c("pair", "size_1st", "size_2nd", "mix_weight_1st", "var_1st", "mix_weight_2nd", "var_2nd")
    }
    df <- do.call(rbind, lapply(rows, function(r) {
      missing <- setdiff(all_cols, names(r))
      r[missing] <- NA
      r[all_cols]
    }))
    print(df, row.names = FALSE, right = TRUE)

  } else {
    n_samples <- length(x$estim_objects)
    cat(" Number of samples:", n_samples, "\n\n")
    sample_names <- x$sample_names
    if (is.null(sample_names)) { sample_names <- names(x$estim_objects) }
    if (is.null(sample_names) || any(sample_names == "")) { sample_names <- paste0("Sample_", seq_len(n_samples)) }
    weights <- sapply(x$estim_objects, function(obj) {
      format(round(obj$estimated_mixing_weights, 3), nsmall = 3)
    })
    sizes <- sapply(x$estim_objects, function(obj) { obj$population_sizes })
    df <- data.frame(sample = sample_names, size = sizes, `mix_weight` = weights, check.names = FALSE)
    if (inherits(x, "estim_BVdk")) {
      df$location <- sapply(x$estim_objects, function(obj) { format(round(obj$estimated_locations, 2), nsmall = 2) })
      if (!is.na(x$estim_objects[[1]]$mix_weight_variance) && !is.na(x$estim_objects[[1]]$location_variance)) {
        df$var_mix_weight <- sapply(x$estim_objects, function(obj) { format(round(obj$mix_weight_variance, 5), nsmall = 5) })
        df$var_location <- sapply(x$estim_objects, function(obj) { format(round(obj$location_variance, 5), nsmall = 5) })
      }
    }
    print(df, row.names = FALSE, right = TRUE)
  }

#  if (inherits(x, "estim_PS")) {
#    cat("\n Use `?estim_PS` for details on the penalization term.\n")
#  } else if (inherits(x, "estim_BVdk")) {
#    cat("\n Use `?estim_BVdk` for details on the optimization method.\n")
#  } else {
#    cat("\n Use `?estim_IBM` for further details.\n")
#  }
  invisible(x)
}


#' Summary method for object of class \code{admix_estim}
#'
#' Summarize the estimated weight(s) of the unknown component(s), and admixture model(s) under study.
#' Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with G a known CDF, and p and F unknown quantities.
#'
#' @param object An object of class \code{admix_estim} (see ?admix_estim).
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.admix_estim <- function(object, ...) {

  cat("\nCall:\n")
  print(object$call)

  method <- sub("estim_", "", class(object)[1])
  n_samples <- length(object$estim_objects)
  sample_names <- object$sample_names
  if (is.null(sample_names) || any(sample_names == "")) {
    sample_names <- paste0("Sample_", seq_len(n_samples))
  }

  if (inherits(object, "estim_IBM")) {
    ## Pairwise:
    valid_names   <- sample_names[-1]
    for (k in seq_along(object$estim_objects)) {
      cat("\n==============================================\n")
      cat("Pair:", sample_names[1], "vs", valid_names[k], "\n")
      summary(object$estim_objects[[k]], show.call = FALSE, ...)
    }

  } else {
    for (k in seq_along(object$estim_objects)) {
      cat("\n==============================================\n")
      cat("Sample:", sample_names[k], "\n")
      summary(object$estim_objects[[k]], show.call = FALSE, ...)
    }
  }
  invisible(object)
}
