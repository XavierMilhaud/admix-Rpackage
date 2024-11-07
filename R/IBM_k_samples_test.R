#' Equality test of K unknown component distributions
#'
#' Equality test of the unknown component distributions coming from K (K > 1) admixture models, based on the Inversion - Best
#' Matching (IBM) approach. Recall that we have K populations following admixture models, each one with probability
#' density functions (pdf) l_k = p_k*f_k + (1-p_k)*g_k, where g_k is the known pdf and l_k corresponds to the
#' observed sample. Perform the following hypothesis test:
#'    H0 : f_1 = ... = f_K  against  H1 : f_i differs from f_j (i different from j, and i,j in 1,...,K).
#'
#' @param samples A list of the K samples to be studied, all following admixture distributions.
#' @param admixMod A list of objects of class class \link[admix]{admix_model}, containing useful information about distributions and parameters.
#' @param conf_level The confidence level of the K-sample test.
#' @param sim_U (default to NULL) Random draws of the inner convergence part of the contrast as defined in the IBM approach (see 'Details' below).
#' @param tune_penalty (default to TRUE) A boolean that allows to choose between a classical penalty term or an optimized penalty (embedding
#'                     some tuning parameters, automatically optimized). Optimized penalty is very useful for low or unbalanced sample sizes
#'                     to detect alternatives to the null hypothesis (H0).
#' @param n_sim_tab (default to 100) Number of simulated Gaussian processes when tabulating the inner convergence distribution
#'                  in the 'icv' testing method using the IBM estimation approach.
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when paralleling computations.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024b}{admix}
#'
#' @return An object of class \link[admix]{admix_test}, containing 17 attributes: 1) the number of samples for the test; 2) the sizes of each sample;
#'         3) the information about component distributions for each sample; 4) the reject decision of the test; 5) the confidence level
#'         of the test (1-alpha, where alpha refers to the first-type error); 6) the test p-value; 7) the 95th-percentile of the contrast
#'         tabulated distribution; 8) the test statistic value; 9) the selected rank (number of terms involved in the test statistic);
#'         10) the terms (pairwise contrasts) involved in the test statistic; 11) A boolean indicating whether the penalization corresponds
#'         to the null hypothesis has been considered; 12) the value of penalized test statistics; 13) the selected optimal 'gamma' parameter
#'         (see reference below); 14) the selected optimal constant involved in the penalization process (see also the reference); 15) the
#'         tabulated distribution of the contrast; 16) the estimated mixing proportions (not implemented yet, since that makes sense only
#'         in case of equal unknown component distributions); 17) the matrix of pairwise contrasts (distance between two samples).
#'
#' @examples
#' \donttest{
#' ####### Under the null hypothesis H0 (with K=3 populations):
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 450, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 1, "sd" = 1)))
#' mixt3 <- twoComp_mixt(n = 400, weight = 0.8,
#'                       comp.dist = list("norm", "exp"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("rate" = 1)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' data3 <- getmixtData(mixt3)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
#'                          knownComp_param = mixt3$comp.param[[2]])
#' ## Perform the 3-samples test:
#' IBM_k_samples_test(samples = list(data1, data2, data3),
#'                    admixMod = list(admixMod1, admixMod2, admixMod3),
#'                    conf_level = 0.95, parallel = FALSE, n_cpu = 2,
#'                    sim_U = NULL, n_sim_tab = 8, tune_penalty = FALSE)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_k_samples_test <- function(samples, admixMod, conf_level = 0.95, sim_U = NULL,
                               tune_penalty = TRUE, n_sim_tab = 100, parallel = FALSE, n_cpu = 2)
{
  old_options_warn <- base::options()$warn
  base::options(warn = -1)

  ## Control whether parallel computations were asked for or not:
  if (parallel) {
    `%fun%` <- doRNG::`%dorng%`
    doParallel::registerDoParallel(cores = n_cpu)
  } else {
    `%fun%` <- foreach::`%do%`
  }

  if (length(samples) == 2) {
    return(IBM_2samples_test(samples = samples, admixMod = admixMod, conf_level = conf_level,
                             sim_U = sim_U, n_sim_tab = n_sim_tab, parallel = parallel, n_cpu = n_cpu))
  } else {

    ## Get the minimal size among all sample sizes, useful for contrast tabulation (adjustment of covariance terms):
    minimal_size <- min(sapply(X = samples, FUN = length))
    if (tune_penalty) {

      ##------ Calibration of tuning parameter 'gamma' --------##
      ## Split each population 6 times in two subsamples:
      S_ii <- matrix(NA, nrow = length(samples), ncol = 6)
      for (k in 1:length(samples)) {

        #sapply(X = admixMod[as.numeric(model.list[[k]])], "[[", "comp.dist")["known", ]
        #sapply(X = admixMod[as.numeric(model.list[[k]])], "[[", "comp.param")["known", ]
        ## Artificially create 6 times two subsamples from one given original sample (thus under the null H0):
        subsample1_index <- matrix(NA, nrow = floor(length(samples[[k]])/2), ncol = 6)
        subsample1_index <- replicate(n = 6, expr = sample(x = 1:length(samples[[k]]),
                                                           size = floor(length(samples[[k]])/2), replace = FALSE, prob = NULL))
        subsample2_index <- matrix(NA, nrow = (length(samples[[k]])-floor(length(samples[[k]])/2)), ncol = 6)
        x_val <- seq(from = min(samples[[k]]), to = max(samples[[k]]), length.out = 1000)
        ## Parallel (or not!) computations of the supremum of the difference b/w decontaminated cdf:
        l <- NULL
        S <-
          foreach::foreach(l = 1:ncol(subsample2_index), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
            subsample2_index[ ,l] <- c(1:length(samples[[k]]))[-subsample1_index[ ,l]]
            ## Comparison between the two populations :
            XY <- estim_IBM(samples = list(samples[[k]][subsample1_index[ ,l]], samples[[k]][subsample2_index[ ,l]]),
                            admixMod = list(admixMod[[k]], admixMod[[k]]), n.integ = 1000)
            F_i1 <- decontaminated_cdf(sample1 = samples[[k]][subsample1_index[ ,l]], estim.p = XY$estimated_mixing_weights[1],
                                       admixMod = admixMod[[k]])
            F_i2 <- decontaminated_cdf(sample1 = samples[[k]][subsample2_index[ ,l]], estim.p = XY$estimated_mixing_weights[1],
                                       admixMod = admixMod[[k]])

            ## Assessment of the difference of interest at each specified x-value:
            max( sqrt(length(samples[[k]])/2) * abs(F_i1(x_val) - F_i2(x_val)) )
          }
        S_ii[k, which(sapply(X = S, FUN = class) != "numeric")] <- NA
        S_ii[k, which(sapply(X = S, FUN = class) == "numeric")] <- unlist(S[which(sapply(X = S, FUN = class) == "numeric")])
      }
      ## Extraction of the optimal gamma :
      gamma_opt <- max( apply(S_ii, 1, max, na.rm = TRUE) / log(sapply(X = samples, FUN = length)/2) )

      ##------ Calibration of tuning parameter 'C' --------##
      couples.expr <- couples.param <- vector(mode = "list", length = length(samples))
      ## Split each population into three subsets (thus leading to be under the null):
      U_k <- vector(mode = "list", length = length(samples))
      ## Parallel (or not!) computations of the embedded test statistics for each original sample:
      U_k <-
        foreach::foreach (i = 1:length(samples), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
          ## Artificially create three subsamples from one given sample, thus under the null H0:
          subsample1_index <- sample(x = 1:length(samples[[i]]), size = floor(length(samples[[i]])/3), replace = FALSE, prob = NULL)
          subsample2_index <- sample(x = c(1:length(samples[[i]]))[-subsample1_index],
                                     size = length(subsample1_index), replace = FALSE, prob = NULL)
          subsample3_index <- c(1:length(samples[[i]]))[-sort(c(subsample1_index, subsample2_index))]
          ## Compute the test statistics for each combination of subsets :
          XY <- estim_IBM(samples = list(samples[[i]][subsample1_index], samples[[i]][subsample2_index]),
                          admixMod = list(admixMod[[i]], admixMod[[i]]), n.integ = 1000)
          if (is.numeric(XY$estimated_mixing_weights)) {
            T_12 <- (length(samples[[i]])/3) *
              IBM_empirical_contrast(par = XY$estimated_mixing_weights, samples = list(samples[[i]][subsample1_index], samples[[i]][subsample2_index]),
                                     admixMod = list(admixMod[[i]], admixMod[[i]]),
                                     G = XY$integ.supp, fixed.p.X = XY$p.X.fixed)
          } else T_12 <- NA
          XZ <- estim_IBM(samples = list(samples[[i]][subsample1_index], samples[[i]][subsample3_index]),
                          admixMod = list(admixMod[[i]], admixMod[[i]]), n.integ = 1000)
          if (is.numeric(XZ$estimated_mixing_weights)) {
            T_13 <- (length(samples[[i]])/3) *
              IBM_empirical_contrast(par = XZ$estimated_mixing_weights, samples = list(samples[[i]][subsample1_index], samples[[i]][subsample3_index]),
                                     admixMod = list(admixMod[[i]], admixMod[[i]]),
                                     G = XZ$integ.supp, fixed.p.X = XZ[["p.X.fixed"]])
          } else T_13 <- NA
          YZ <- estim_IBM(samples = list(samples[[i]][subsample2_index], samples[[i]][subsample3_index]),
                          admixMod = list(admixMod[[i]], admixMod[[i]]), n.integ = 1000)
          if (is.numeric(YZ$estimated_mixing_weights)) {
            T_23 <- (length(samples[[i]])/3) *
              IBM_empirical_contrast(par = YZ$estimated_mixing_weights, samples = list(samples[[i]][subsample2_index], samples[[i]][subsample3_index]),
                                     admixMod = list(admixMod[[i]], admixMod[[i]]),
                                     G = YZ$integ.supp, fixed.p.X = YZ[["p.X.fixed"]])
          } else T_23 <- NA
          cumsum(x = c(T_12, T_13, T_23))
        }
      ## Remove elements of the list where an error was stored:
      U_k[which(sapply(X = U_k, FUN = function(x) any(is.na(x))) == TRUE)] <- NULL
      U_k[which(sapply(X = U_k, FUN = class) != "numeric")] <- NULL
      if (length(U_k) == 0) stop("The calibration of tuning parameters in the k-sample test procedure failed. Please try again.")
      ## We are artificially under the null H0, thus fix epsilon close to 1:
      epsilon_opt <- 0.99
      ## Determine the optimal constant, applying the rule leading to select the first term of the embedded statistic under H0 :
      constant_list <- vector(mode = "numeric", length = length(U_k))
      for (i in 1:length(U_k)) {
        constant_list[i] <- max( c( (U_k[[i]][2]-U_k[[i]][1])/((length(samples[[i]])/length(U_k[[i]]))^epsilon_opt),
                                    (U_k[[i]][3]-U_k[[i]][1])/(2*(length(samples[[i]])/length(U_k[[i]]))^epsilon_opt) ) )
      }
      #cst_selected <- max(constant_list)
      cst_selected <- min(constant_list)

    } else {
      gamma_opt <- NA
      cst_selected <- NA
    }

    #################################################################
    ##----------- Back to the general k-sample procedure ----------##
    couples.list <- NULL
    for (i in 1:(length(samples)-1)) { for (j in (i+1):length(samples)) { couples.list <- rbind(couples.list,c(i,j)) } }

    S_ij <- vector(mode = "list", length = nrow(couples.list))
    ## Compute the supremum between decontaminated cdfs over the different couples of populations under study:
    S_ij <-
      foreach::foreach (k = 1:nrow(couples.list), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
        ## Comparison between the two populations :
        XY <- estim_IBM(samples = samples[as.numeric(couples.list[k, ])], admixMod = admixMod[as.numeric(couples.list[k, ])], n.integ = 1000)
        emp.contr <- minimal_size * IBM_empirical_contrast(XY$estimated_mixing_weights,
                                                           samples = samples[as.numeric(couples.list[k, ])],
                                                           admixMod = admixMod[as.numeric(couples.list[k, ])],
                                                           G = XY$integ.supp, fixed.p.X = XY$p.X.fixed)
        x_val <- seq(from = min(samples[[couples.list[k, ][1]]], samples[[couples.list[k, ][2]]]),
                     to = max(samples[[couples.list[k, ][1]]], samples[[couples.list[k, ][2]]]), length.out = 1000)
        if (length(XY$estimated_mixing_weights) == 2) {
          F_1 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][1]]], estim.p = XY$estimated_mixing_weights[1],
                                    admixMod = admixMod[[couples.list[k, ][1]]])
          F_2 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][2]]], estim.p = XY$estimated_mixing_weights[2],
                                    admixMod = admixMod[[couples.list[k, ][2]]])
        } else {
          F_1 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][1]]], estim.p = XY$p.X.fixed,
                                    admixMod = admixMod[[couples.list[k, ][1]]])
          F_2 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][2]]], estim.p = XY$estimated_mixing_weights,
                                    admixMod = admixMod[[couples.list[k, ][2]]])
        }
        ## Assessment of the supremum:
        list( sup = max( sqrt(minimal_size) * abs(F_1(x_val) - F_2(x_val)) ),
              contrast = emp.contr)
      }
    list_sup <- sapply(X = S_ij, "[[", "sup")
    list_empirical.contr <- sapply(X = S_ij, "[[", "contrast")
    ## Remove elements of the list where an error was stored:
    to_remove <- which(sapply(X = list_sup, FUN = is.null))
    list_empirical.contr[to_remove] <- list_sup[to_remove] <- NULL
    if ((length(list_empirical.contr) == 0) | (length(list_sup) == 0)) {
      stop("The calibration of tuning parameters in the k-sample test procedure has failed. Please try again.")
    }
    max.S_ij <- max(unlist(list_sup))

    ## Choice of the penalty rule:
    if (tune_penalty) {
      ## Define the two possible exponents to apply to the sample size :
      epsilon_null <- 0.99
      epsilon_alt <- 0.75
      penalty_rule <- (max.S_ij <= (gamma_opt * log(minimal_size)))
      ## Define the penalty function, depending on the penalty rule:
      penalty <- penalty_rule * cst_selected * minimal_size^epsilon_null + (1-penalty_rule) * cst_selected * minimal_size^epsilon_alt
    } else {
      ## Apply the simple n^epsilon, taking epsilon right in the]0.5,0.9[:
      penalty_rule <- NA
      penalty <- minimal_size^0.87
    }

    ## Give a simple and useful representation of results for each couple:
    contrast.matrix <- discrepancy.id <- discrepancy.rank <- matrix(NA, nrow = length(samples), ncol = length(samples))
    ranks <- 1:nrow(couples.list)
    if (length(to_remove) != 0) {
      list.iterators <- c(1:nrow(couples.list))[-to_remove]
      ranks[to_remove] <- NA
    } else {
      list.iterators <- 1:nrow(couples.list)
    }
    ranks[which(!is.na(ranks))] <- rank(unlist(list_empirical.contr))
    for (k in list.iterators) {
      contrast.matrix[couples.list[k, ][1], couples.list[k, ][2]] <- sapply(X = S_ij, "[[", "contrast")[[k]]
      discrepancy.id[couples.list[k, ][1], couples.list[k, ][2]] <- paste(couples.list[k, ][1], couples.list[k, ][2], sep = "-")
      discrepancy.rank[couples.list[k, ][1], couples.list[k, ][2]] <- ranks[k]
    }

    ## Tabulate the contrast distribution to retrieve the quantile of interest, needed further to perform the test:
    if (is.null(sim_U)) {
      which_row <- unlist(apply(discrepancy.rank, 2, function(x) which(x == 1)))
      which_col <- unlist(apply(discrepancy.rank, 1, function(x) which(x == 1)))
      U <- IBM_tabul_stochasticInteg(samples = list(samples[[which_row]], samples[[which_col]]),
                                     admixMod = list(admixMod[[which_row]], admixMod[[which_col]]),
                                     min_size = minimal_size, n.varCovMat = 80, n_sim_tab = n_sim_tab,
                                     parallel = parallel, n_cpu = n_cpu)
      sim_U <- U[["U_sim"]]
    } else {
      sim_U <- sim_U
    }
    if (all(is.na(sim_U))) {
      stop("Impossible to tabulate the distribution of the contrast!")
      return( list(rejection_rule = TRUE, p_val = NA, stat_name = NA, stat_value = NA, stat_rank = NA,
                   penalized_stat = NA, ordered_contrasts = NA, quantiles = NA, stat_terms = NA,
                   contr_mat = contrast.matrix) )
    } else {
      q_H <- stats::quantile(sim_U, conf_level)
    }

    ##--------- Perform the hypothesis testing -----------##
    ## Sort the discrepancies and corresponding information:
    sorted_contrasts <- sort(unlist(list_empirical.contr))
    penalized.stats <- numeric(length = (length(sorted_contrasts)+1))
    stat.terms <- stat.terms_final <- character(length = (nrow(couples.list)-length(to_remove)))
    for (k in 1:length(stat.terms)) {
      stat.terms[k] <- paste("d_", discrepancy.id[which(discrepancy.rank == k)], sep = "")
      stat.terms_final[k] <- paste(stat.terms[1:k], collapse = " + ", sep = "")
      penalized.stats[k+1] <- penalized.stats[k] + sorted_contrasts[k] - penalty
    }
    penalized.stats <- round(penalized.stats[-1], 1)
    ## Selection of the number of embedded terms in the statistic, and compute the statistic itself:
    selected.index <- which.max(penalized.stats)
    finalStat_name <- stat.terms_final[selected.index]
    finalStat_value <- cumsum(sorted_contrasts)[selected.index]
    ## Test decision:
    final_test <- finalStat_value > q_H[1]
    ## Corresponding p-value:
    CDF_U <- stats::ecdf(sim_U)
    p_value <- 1 - CDF_U(finalStat_value)

    obj <- list(
      n_populations = length(samples),
      population_sizes = sapply(X = samples, FUN = length),
      admixture_models = admixMod,
      reject_decision = final_test,
      confidence_level = conf_level,
      p_value = p_value,
      extreme_quantile_tabul = unique(q_H),
      test_statistic_value = finalStat_value,
      selected_rank = selected.index,
      statistic_name = finalStat_name,
      penalty_nullHyp = penalty_rule,
      penalized_stat = penalized.stats,
      tuned_gamma = gamma_opt,
      tuned_constant = cst_selected,
      tabulated_dist = sim_U,
      estimated_mixing_weights = NA,
      contrast_matrix = contrast.matrix
    )
    class(obj) <- c("IBM_test", "admix_test")
    obj$call <- match.call()

    on.exit(base::options(warn = old_options_warn))
    return(obj)
  }

}

#' Print method for objects 'IBM_test'
#'
#' @param x An object of class 'IBM_test'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.IBM_test <- function(x, ...)
{
  cat("Call:")
  print(x$call)
  cat("\n")
  cat("Is the null hypothesis (equal unknown component distributions) rejected? ",
      ifelse(x$reject_decision, "Yes", "No"), sep="")
  cat("\nTest p-value: ", round(x$p_value,3), "\n", sep="")
}


#' Summary method for objects 'IBM_test'
#'
#' @param object An object of class 'IBM_test'.
#' @param ... A list of additional parameters belonging to the default method.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.IBM_test <- function(object, ...)
{
  cat("Call:")
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
  cat("\n------- Test decision -------\n")
  cat("Is the null hypothesis (equality of unknown component distributions) rejected? ",
      ifelse(object$reject_decision, "Yes", "No"), sep="")
  cat("\nConfidence level of the test: ", object$confidence_level, sep="")
  cat("\nTest p-value: ", round(object$p_value,3), sep="")
  cat("\n\n------- Test statistic -------\n")
  cat("Selected rank of the test statistic (following the penalization rule): ", object$selected_rank, sep="")
  cat("\nValue of the test statistic: ", round(object$test_statistic_value,2), "\n", sep="")
  cat("Discrepancy terms involved in the statistic: ", paste(object$statistic_name, sep = ""), "\n", sep = "")
  cat("Optimal tuning parameters involved in the test statistic (if argument 'tune.penalty' is true):\n")
  cat("Gamma: ", object$tuned_gamma, "\n", sep = "")
  cat("Constant: ", object$tuned_constant, "\n", sep = "")
  cat("Chosen penalty rule: ", ifelse(object$penalty_nullHyp, "H0", "H1"), sep = "")
  cat("\n\n------- Tabulated test statistic distribution -------\n")
  cat("Quantile at level ", object$confidence_level*100, "%: ", round(object$extreme_quantile_tabul, 3), "\n", sep = "")
  cat("Tabulated distribution: ", paste(utils::head(round(sort(object$tabulated_dist),2),3), collapse = " "), "....",
      paste(utils::tail(round(sort(object$tabulated_dist),2),3), collapse = " "), "\n", sep = "")
  cat("\n")
}


#' Test equality of 2 unknown component distributions using IBM
#'
#' Two-sample equality test of the unknown component distributions coming from two admixture models, using Inversion - Best Matching
#' (IBM) approach. Recall that we have two admixture models with respective probability density functions (pdf)
#' l1 = p1 f1 + (1-p1) g1 and l2 = p2 f2 + (1-p2) g2, where g1 and g2 are known pdf and l1 and l2 are observed.
#' Perform the following hypothesis test: H0 : f1 = f2 versus H1 : f1 differs from f2.
#'
#' @param samples A list of the two samples under study.
#' @param admixMod A list of two objects of class 'admix_model', containing useful information about distributions and parameters.
#' @param conf_level (default to 0.95) The confidence level of the 2-samples test, i.e. the quantile level to which the test statistic is compared.
#' @param sim_U Random draws of the inner convergence part of the contrast as defined in the IBM approach (see 'Details' below).
#' @param n_sim_tab (default to 100) Number of simulated Gaussian processes when tabulating the inner convergence distribution
#'                  in the 'icv' testing method using the IBM estimation approach.
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when paralleling computations.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return A list of 6 elements, containing : 1) the rejection decision; 2) the confidence level of the test (1-alpha, where
#'         alpha stands for the first-type error); 3) the p-value of the test, 4) the test statistic value; 5) the simulated
#'         distribution of the inner convergence regime (useful to perform the test when comparing to the extreme quantile of
#'         this distribution; 6) the estimated weights of the unknown components for each of the two admixture models.
#'
#' @examples
#' ####### Under the alternative H1:
#' mixt1 <- twoComp_mixt(n = 450, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 1, "sd" = 0.5),
#'                                         list("mean" = 1, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' IBM_2samples_test(samples = list(data1, data2),
#'                   admixMod = list(admixMod1, admixMod2), conf_level = 0.95,
#'                   parallel = FALSE, n_cpu = 2, sim_U = NULL, n_sim_tab = 10)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

IBM_2samples_test <- function(samples, admixMod, conf_level = 0.95, parallel = FALSE,
                              n_cpu = 2, sim_U = NULL, n_sim_tab = 100)
{
  min_size <- min(length(samples[[1]]), length(samples[[2]]))
  ## Estimate the proportions of the mixtures:
  estim <- try(suppressWarnings(estim_IBM(samples = samples, admixMod = admixMod, n.integ = 1000)), silent = TRUE)
  count_error <- 0
  while ((inherits(x = estim, what = "try-error", which = FALSE)) & (count_error < 3)) {
    estim <- NULL
    estim <- try(suppressWarnings(estim_IBM(samples = samples, admixMod = admixMod, n.integ = 1000)), silent = TRUE)
    count_error <- count_error + 1
  }

  if (inherits(x = estim, what = "try-error", which = FALSE)) {
    estim.weights <- 100
    contrast_val <- NA
  } else {
    estim.weights <- estim$estimated_mixing_weights
    min_size <- min(length(samples[[1]]), length(samples[[2]]))
    contrast_val <- min_size * IBM_empirical_contrast(par = estim.weights, samples = samples, admixMod = admixMod,
                                                      G = estim$integ.supp, fixed.p.X = estim$p.X.fixed)
    ## Identify boolean 'green light' criterion to known whether we need to perform the test with stochastic integral tabulation:
    # green_light <- IBM_greenLight_criterion(estim_obj = estim, samples = samples, admixMod = admixMod, alpha = (1-conf_level))
    # if (green_light) {
    #   if (is.null(sim_U)) {
    #     tab_dist <- IBM_tabul_stochasticInteg(samples = samples, admixMod = admixMod, min_size = min_size, n.varCovMat = 80,
    #                                           n_sim_tab = n_sim_tab, parallel = parallel, n_cpu = n_cpu)
    #     sim_U <- tab_dist$U_sim
    #     contrast_val <- tab_dist$contrast_value
    #   }
    #   reject.decision <- contrast_val > quantile(sim_U, 0.95)
    #   CDF_U <- ecdf(sim_U)
    #   p_value <- 1 - CDF_U(contrast_val)
    # } else {
    #   reject.decision <- TRUE
    #   p_value <- 1e-16
    # }
  }

  ## Earn computation time using this soft version of the green light criterion:
  if (any(abs(estim.weights) > 1)) {
    reject <- TRUE
    p_value <- 1e-16
    sim_U <- extreme_quantile <- NA
  } else {
    if (is.null(sim_U)) {
      tab_dist <- IBM_tabul_stochasticInteg(samples = samples, admixMod = admixMod, min_size = min_size, n.varCovMat = 80,
                                            n_sim_tab = n_sim_tab, parallel = parallel, n_cpu = n_cpu)
      sim_U <- tab_dist$U_sim
      contrast_val <- tab_dist$contrast_value
    }
    extreme_quantile <- stats::quantile(sim_U, conf_level)
    reject <- contrast_val > extreme_quantile
    CDF_U <- stats::ecdf(sim_U)
    p_value <- 1 - CDF_U(contrast_val)
  }

  obj <- list(
    reject_decision = reject,
    confidence_level = conf_level,
    p_value = p_value,
    extreme_quantile_tabul = extreme_quantile,
    test_statistic_value = contrast_val,
    selected_rank = NA,
    statistic_name = NA,
    penalty_nullHyp = NA,
    penalized_stat = NA,
    tuned_gamma = NA,
    tuned_constant = NA,
    tabulated_dist = sim_U,
    estimated_mixing_weights = estim.weights,
    contrast_matrix = NA
  )
  return(obj)
}


#' Green-light criterion for equality tests in admixture models
#'
#' Indicates whether one needs to perform the full version of the statistical test of equality between unknown component
#' distributions in two samples following admixture models (it is sometimes obvious that they differ), based on the IBM approach.
#' See more in 'Details' below.
#'
#' @param estim_obj Object of class 'estim_IBM', obtained from the estimation of the component weights related to the
#'                  proportions of the unknown component in each of the two admixture models studied.
#' @param samples (list) A list of the two samples under study.
#' @param admixMod (list) A list of two objects of class 'admix_model', containing useful information about distributions and parameters.
#' @param alpha Confidence region is defined by the probability (1-alpha), used to compute the confidence bands of the estimators
#'              of the unknown component weights.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return A boolean indicating whether it is useful or not to tabulate the contrast distribution in order to answer
#'         the testing problem (H0: f1 = f2).
#'
#' @examples
#' \donttest{
#' ####### Under the null hypothesis H0 (with K=3 populations):
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1200, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 1000, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 1, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
#' estim <- estim_IBM(samples = list(data1, data2),
#'                    admixMod = list(admixMod1, admixMod2), n.integ = 1000)
#' IBM_greenLight_criterion(estim_obj = estim, samples = list(data1, data2),
#'                          admixMod = list(admixMod1, admixMod2), alpha = 0.05)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

IBM_greenLight_criterion <- function(estim_obj, samples, admixMod, alpha = 0.05)
{
  if (length(samples) != 2) stop("There should be exactly TWO samples under study.")

  min_sample_size <- min(length(samples[[1]]), length(samples[[2]]))
  length.support <- length(estim_obj$integ.supp)
  z <- estim_obj$integ.supp[round(floor(length.support/2))]
  varCov_estim <- IBM_estimVarCov_gaussVect(x = z, y = z, IBMestim.obj = estim_obj, samples = samples, admixMod = admixMod)

  if (length(estim_obj$estimated_mixing_weights) == 2) {
    inf_bound.p1 <- estim_obj$estimated_mixing_weights[1] - sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p1 <- estim_obj$estimated_mixing_weights[1] + sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p1 <- c(inf_bound.p1, sup_bound.p1)
    inf_bound.p2 <- estim_obj$estimated_mixing_weights[2] - sqrt(varCov_estim[2,2]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p2 <- estim_obj$estimated_mixing_weights[2] + sqrt(varCov_estim[2,2]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p2 <- c(inf_bound.p2, sup_bound.p2)
    green_light_crit <- max(conf_interval.p1[1],conf_interval.p2[1]) <= 1
  } else {
    inf_bound.p <- estim_obj$estimated_mixing_weights[1] - sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p <- estim_obj$estimated_mixing_weights[1] + sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p2 <- c(inf_bound.p, sup_bound.p)
    conf_interval.p1 <- NULL
    green_light_crit <- conf_interval.p2[1] <= 1
  }

  ret <- list(
    green_light = green_light_crit,
    conf_interval_p1 = conf_interval.p1,
    conf_interval_p2 = conf_interval.p2
  )
  return(ret)
}


#' Simulated distribution of the contrast using IBM
#'
#' Tabulate the distribution related to the inner convergence part of the contrast, by simulating trajectories of Gaussian
#' processes and applying some transformations. Useful to perform the test hypothesis, by retrieving the (1-alpha)-quantile
#' of interest. See 'Details' below and the cited paper therein for further information.
#'
#' @param samples A list of the two samples under study.
#' @param admixMod A list of two objects of class 'admix_model', with information about distributions and parameters.
#' @param min_size (optional, NULL by default) In the k-sample case, useful to provide the minimal size among all samples
#'                 (needed to take into account the correction factor for variance-covariance assessment). Otherwise, useless.
#' @param n.varCovMat (default to 80) Number of time points at which the Gaussian processes are simulated.
#' @param n_sim_tab (default to 100) Number of simulated Gaussian processes when tabulating the inner convergence distribution
#'                  in the 'icv' testing method using the IBM estimation approach.
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when paralleling computations.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return A list with four elements, containing: 1) random draws of the contrast as defined in the reference given here;
#'         2) estimated unknown component weights for the two admixture models; 3) the value of the the empirical contrast;
#'         4) support that was used to evaluate the variance-covariance matrix of the empirical processes.
#'
#' @examples
#' \donttest{
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1200, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 1000, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 1, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' IBM_tabul_stochasticInteg(samples = list(data1, data2), admixMod = list(admixMod1, admixMod2),
#'                           min_size=NULL, n.varCovMat=20, n_sim_tab=2, parallel=FALSE, n_cpu=2)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_tabul_stochasticInteg <- function(samples, admixMod, min_size = NULL, n.varCovMat = 80,
                                      n_sim_tab = 100, parallel = FALSE, n_cpu = 2)
{
  stopifnot("Wrong number of samples... Must be 2!" = length(samples) == 2)

  i <- NULL
  if (parallel) {
    `%fun%` <- doRNG::`%dorng%`
    doParallel::registerDoParallel(cores = n_cpu)
  } else {
    `%fun%` <- foreach::`%do%`
  }

  ## Extract the information on component distributions:
  knownCDF_comp.dist <- paste0("p", unlist(sapply(admixMod, '[[', 'comp.dist')["known", ]))
  if (any(knownCDF_comp.dist == "pmultinom")) knownCDF_comp.dist[which(knownCDF_comp.dist == "pmultinom")] <- "stepfun"
  comp_emp <- sapply(X = knownCDF_comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_emp)) assign(x = names(comp_emp)[i], value = comp_emp[[i]])

  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_emp)[i],"(x = 1:", length(admixMod[[i]]$comp.param$known$prob),
                                      paste(", y = ", paste("cumsum(c(0,", paste(admixMod[[i]]$comp.param$known$prob, collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_emp)[i],"(z,", paste(names(admixMod[[i]]$comp.param$known),
                                                                 "=", admixMod[[i]]$comp.param$known, sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(knownCDF_comp.dist))
  expr[which(knownCDF_comp.dist == "stepfun")] <- sapply(which(knownCDF_comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  if (any(knownCDF_comp.dist == "stepfun")) {
    G1.fun <- eval(parse(text = expr[1]))
    G2.fun <- eval(parse(text = expr[2]))
    G1 <- function(z) G1.fun(z)
    G2 <- function(z) G2.fun(z)
  } else {
    G1 <- function(z) { eval(parse(text = expr[1])) }
    G2 <- function(z) { eval(parse(text = expr[2])) }
  }

  ## Empirical cumulative distribution function from the two observed samples:
  L1 <- stats::ecdf(samples[[1]])
  L2 <- stats::ecdf(samples[[2]])

  ## Estimate the weights of the unknown component distributions in first and second samples:
  estim <- estim_IBM(samples = samples, admixMod = admixMod, n.integ = 1000)
  if (is.null(min_size)) {
    sample.size <- min(length(samples[[1]]), length(samples[[2]]))
  } else {
    sample.size <- min_size
  }
  contrast_val <- sample.size *
    IBM_empirical_contrast(par = estim$estimated_mixing_weights, samples = samples, admixMod = admixMod,
                           G = estim$integ.supp, fixed.p.X = estim$p.X.fixed)
  ## Integration support:
  support <- detect_support_type(samples[[1]], samples[[2]])
  if (support == "Continuous") {
    densite.G <- stats::density(estim$integ.supp, bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integ <- c(min(estim$integ.supp), max(estim$integ.supp))
    t_seq <- seq(from = supp.integ[1], to = supp.integ[2], length.out = n.varCovMat)
  } else {
    ## Case of multinomial distribution :
    supp.integ <- estim$integ.supp
    t_seq <- sort(unique(supp.integ))
  }

  ## Compute the normalization matrix M(.) at each point, to be used further when determining the full simulated trajectories:
  normalization_factors <-
    foreach::foreach (i = 1:length(t_seq), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
      IBM_normalization_term(t_seq[i], estim, samples, admixMod)
    }

  ## Estimate the variance-covariance functions from the empirical processes:
  cov_mat_L1 <- estimVarCov_empProcess_Rcpp(t = t_seq, obs_data = samples[[1]])
  cov_mat_L2 <- estimVarCov_empProcess_Rcpp(t = t_seq, obs_data = samples[[2]])

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(admixMod[[1]], admixMod[[2]])
  if (G1equalG2) {
    psi1 <- function(z) 2*( ((2-estim[["estimated_mixing_weights"]])/estim[["estimated_mixing_weights"]]^3) * G1(z) -
                              (2/estim[["estimated_mixing_weights"]]^3)*L2(z) + (1/(estim[["estimated_mixing_weights"]]^2*estim[["p.X.fixed"]]))*L1(z) -
                              ((1-estim[["p.X.fixed"]])/(estim[["estimated_mixing_weights"]]^2*estim[["p.X.fixed"]])) * G1(z) )
    psi2 <- function(z) 2*( (1/(estim[["estimated_mixing_weights"]]^2*estim[["p.X.fixed"]])) * (L2(z) - G1(z)) )
  } else {
    psi1.1 <- function(z) 2*( ((2-estim[["estimated_mixing_weights"]][1])/estim[["estimated_mixing_weights"]][1]^3) * G1(z) -
                                (2/estim[["estimated_mixing_weights"]][1]^3)*L1(z) + (1/(estim[["estimated_mixing_weights"]][1]^2*estim[["estimated_mixing_weights"]][2]))*L2(z) -
                                ((1-estim[["estimated_mixing_weights"]][2])/(estim[["estimated_mixing_weights"]][1]^2*estim[["estimated_mixing_weights"]][2])) * G2(z) )
    psi1.2 <- function(z) 2*( (1/(estim[["estimated_mixing_weights"]][1]^2*estim[["estimated_mixing_weights"]][2])) * (L1(z) - G1(z)) )
    psi2.1 <- function(z) 2*( ((2-estim[["estimated_mixing_weights"]][2])/estim[["estimated_mixing_weights"]][2]^3) * G2(z) -
                                (2/estim[["estimated_mixing_weights"]][2]^3)*L2(z) + (1/(estim[["estimated_mixing_weights"]][2]^2*estim[["estimated_mixing_weights"]][1]))*L1(z) -
                                ((1-estim[["estimated_mixing_weights"]][1])/(estim[["estimated_mixing_weights"]][2]^2*estim[["estimated_mixing_weights"]][1])) * G1(z) )
    psi2.2 <- function(z) 2*( (1/(estim[["estimated_mixing_weights"]][2]^2*estim[["estimated_mixing_weights"]][1])) * (L2(z) - G2(z)) )
  }

  U_sim <-
    foreach::foreach (i = 1:n_sim_tab, .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
      ## Simulate random gaussian processes using the variance-covariance functions determined from empirical processes:
      B1_path <- sim_gaussianProcess(mean_vec = rep(0,nrow(cov_mat_L1)), varCov_mat = cov_mat_L1, from = supp.integ[1],
                                     to = supp.integ[length(supp.integ)], start = 0, nb.points = nrow(cov_mat_L1))
      B1_traj <- stats::approxfun(x = t_seq, y = B1_path$traj1, method = "linear")
      B2_path <- sim_gaussianProcess(mean_vec = rep(0,nrow(cov_mat_L2)), varCov_mat = cov_mat_L2, from = supp.integ[1],
                                     to = supp.integ[length(supp.integ)], start = 0, nb.points = nrow(cov_mat_L2))
      B2_traj <- stats::approxfun(x = t_seq, y = B2_path$traj1, method = "linear")

      ## Introdude vector 'Z' :
      if (G1equalG2) {

        integrand.phi2 <- function(z) {
          if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
          }
          psi2(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi1 <- function(z) {
          if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
          }
          psi1(z) * B2_traj(z) * densite.G.dataPoint
        }
        Z <- matrix(NA, nrow = 4, ncol = length(t_seq))
        if (support == "Continuous") {
          Z[1, ] <- stats::integrate(Vectorize(integrand.phi2), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[2, ] <- B1_path$traj1
          Z[3, ] <- stats::integrate(Vectorize(integrand.phi1), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[4, ] <- B2_path$traj1
        } else {
          Z[1, ] <- sum( unlist(sapply(supp.integ, integrand.phi2)) )
          Z[2, ] <- B1_path$traj1
          Z[3, ] <- sum( unlist(sapply(supp.integ, integrand.phi1)) )
          Z[4, ] <- B2_path$traj1
        }
        estim.vector <- matrix(NA, nrow = 2, ncol = length(t_seq))

      } else {

        integrand.phi1.1 <- function(z) {
          if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
          }
          psi1.1(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi2.2 <- function(z) {
          if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
          }
          psi2.2(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi2.1 <- function(z) {
          if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
          }
          psi2.1(z) * B2_traj(z) * densite.G.dataPoint
        }
        integrand.phi1.2 <- function(z) {
          if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
          }
          psi1.2(z) * B2_traj(z) * densite.G.dataPoint
        }
        Z <- matrix(NA, nrow = 6, ncol = length(t_seq))
        if (support == "Continuous") {
          Z[1, ] <- stats::integrate(Vectorize(integrand.phi1.1), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[2, ] <- stats::integrate(Vectorize(integrand.phi2.2), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[3, ] <- B1_path$traj1
          Z[4, ] <- stats::integrate(Vectorize(integrand.phi2.1), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[5, ] <- stats::integrate(Vectorize(integrand.phi1.2), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[6, ] <- B2_path$traj1
        } else {
          Z[1, ] <- sum( unlist(sapply(supp.integ, integrand.phi1.1)) )
          Z[2, ] <- sum( unlist(sapply(supp.integ, integrand.phi2.2)) )
          Z[3, ] <- B1_path$traj1
          Z[4, ] <- sum( unlist(sapply(supp.integ, integrand.phi2.1)) )
          Z[5, ] <- sum( unlist(sapply(supp.integ, integrand.phi1.2)) )
          Z[6, ] <- B2_path$traj1
        }
        estim.vector <- matrix(NA, nrow = 3, ncol = length(t_seq))
      }

      ## Get the trajectory of function 'sqrt(n)*(Dn(.)-D(.))' :
      for (j in 1:length(t_seq)) estim.vector[ ,j] <- normalization_factors[[j]] %*% Z[ ,j]
      D_function <- stats::approxfun(x = t_seq, y = estim.vector[nrow(estim.vector), ], method = "linear")

      integrand <- function(z) {
        if (support == "Continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
        } else {
          densite.G.dataPoint <- 1 / length(table(c(samples[[1]],samples[[2]])))
        }
        D_function(z)^2 * densite.G.dataPoint
      }

      if (support == "Continuous") {
        U <- stats::integrate(integrand, lower = supp.integ[1], upper = supp.integ[2], subdivisions = 10000L, rel.tol = 1e-04)$value
      } else {
        U <- sum( unlist(sapply(supp.integ, integrand)) )
      }
      U
    }

  indexes.toRemove <- which( (substr(U_sim, start = 1, stop = 5) == "Error") | (substr(U_sim, start = 1, stop = 5) == "list(") )
  if (length(indexes.toRemove) != 0) U_sim <- U_sim[-indexes.toRemove]

  obj <- list(
    U_sim = as.numeric(U_sim),
    estimator = estim,
    contrast_value = contrast_val,
    integ.points = t_seq
  )
  return(obj)
}
