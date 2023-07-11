#' Equality test of unknown component distributions in K admixture models, with IBM approach
#'
#' Test hypothesis on the unknown component of K (K > 1) admixture models using Inversion - Best Matching method.
#' K-samples test of the unknown component distribution in admixture models using Inversion - Best Matching
#' (IBM) method. Recall that we have K populations following admixture models, each one with probability
#' density functions (pdf) l_k = p_k*f_k + (1-p_k)*g_k, where g_k is the known pdf and l_k corresponds to the
#' observed sample. Perform the following hypothesis test:
#'    H0 : f_1 = ... = f_K  against  H1 : f_i differs from f_j (i different from j, and i,j in 1,...,K).
#'
#' @param samples A list of the K samples to be studied, all following admixture distributions.
#' @param sim_U (default to NULL) Random draws of the inner convergence part of the contrast as defined in the IBM approach (see 'Details' below).
#' @param n_sim_tab Number of simulated gaussian processes when tabulating the inner convergence distribution in the IBM approach.
#' @param comp.dist A list with 2*K elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the K admixture models. Elements, grouped by 2, refer to the unknown and known components of each admixture model,
#'                  If there are unknown elements, they must be specified as 'NULL' objects. For instance, 'comp.dist' could be specified
#'                  as follows with K = 3: list(f1 = NULL, g1 = 'norm', f2 = NULL, g2 = 'norm', f3 = NULL, g3 = 'rnorm').
#' @param comp.param A list with 2*K elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Elements, grouped by 2, refer to the parameters of unknown and known components of each admixture model.
#'                   If there are unknown elements, they must be specified as 'NULL' objects. For instance, 'comp.param' could
#'                   be specified as follows (with K = 3):
#'                   list(f1 = NULL, g1 = list(mean=0,sd=1), f2 = NULL, g2 = list(mean=3,sd=1.1), f3 = NULL, g3 = list(mean=-2,sd=0.6)).
#' @param conf.level The confidence level of the K-sample test.
#' @param tune.penalty A boolean that allows to choose between a classical penalty term or an optimized penalty embedding some tuning parameters
#'                     (automatically optimized). Optimized penalty is particularly useful for low sample size.
#' @param parallel (default to FALSE) Boolean indicating whether parallel computations are performed.
#' @param n_cpu (default to 2) Number of cores used when parallelizing.
#'
#' @details See the paper at the following HAL weblink: https://hal.science/hal-04129130
#'
#' @return A list of ten elements, containing: 1) the rejection decision; 2) the test p-value; 3) the terms involved in
#'         the test statistic; 4) the test statistic value; 5) the selected rank (number of terms involved in the test statistic);
#'         6) the value of the penalized test statistic; 7) a boolean indicating whether the applied penalty rule is that under
#'         the null H0; 8) the sorted contrast values; 9) the 95th-quantile of the contrast distribution; 10) the final terms of
#'         the statistic; and 11) the contrast matrix.
#'
#' @examples
#' \donttest{
#' ####### Under the null hypothesis H0 (with K=3 populations):
#' ## Specify the parameters of the mixture models for simulation:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm",
#'                   f3 = "norm", g3 = "norm")
#' list.param <- list(f1 = list(mean = 0, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 0, sd = 1), g2 = list(mean = 4, sd = 1.1),
#'                    f3 = list(mean = 0, sd = 1), g3 = list(mean = -3, sd = 0.8))
#' ## Simulate the data:
#' sim1 <- rsimmix(n = 1000, unknownComp_weight = 0.8, comp.dist = list(list.comp$f1,list.comp$g1),
#'                 comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' sim2 <- rsimmix(n= 1300, unknownComp_weight = 0.6, comp.dist = list(list.comp$f2,list.comp$g2),
#'                 comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' sim3 <- rsimmix(n = 1100, unknownComp_weight = 0.7, comp.dist = list(list.comp$f3,list.comp$g3),
#'                 comp.param = list(list.param$f3, list.param$g3))$mixt.data
#' ## Back to the context of admixture models, where one mixture component is unknown:
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm",
#'                   f3 = NULL, g3 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 4, sd = 1.1),
#'                    f3 = NULL, g3 = list(mean = -3, sd = 0.8))
#' ## Perform the 3-samples test:
#' IBM_k_samples_test(samples = list(sim1,sim2,sim3), sim_U= NULL, n_sim_tab = 20,
#'                    comp.dist = list.comp, comp.param = list.param, conf.level = 0.95,
#'                    tune.penalty = FALSE, parallel = FALSE, n_cpu = 2)
#'
#' ####### Now under the alternative H1:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm",
#'                   f3 = "norm", g3 = "norm")
#' list.param <- list(f1 = list(mean = 0, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 0, sd = 1), g2 = list(mean = 4, sd = 1.1),
#'                    f3 = list(mean = 2, sd = 0.7), g3 = list(mean = 3, sd = 0.8))
#' sim1 <- rsimmix(n = 3000, unknownComp_weight = 0.8, comp.dist = list(list.comp$f1,list.comp$g1),
#'                 comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' sim2 <- rsimmix(n= 3300, unknownComp_weight = 0.6, comp.dist = list(list.comp$f2,list.comp$g2),
#'                 comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' sim3 <- rsimmix(n = 3100, unknownComp_weight = 0.7, comp.dist = list(list.comp$f3,list.comp$g3),
#'                 comp.param = list(list.param$f3, list.param$g3))$mixt.data
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm",
#'                   f3 = NULL, g3 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 4, sd = 1.1),
#'                    f3 = NULL, g3 = list(mean = 3, sd = 0.8))
#' IBM_k_samples_test(samples = list(sim1,sim2,sim3), sim_U= NULL, n_sim_tab = 20,
#'                    comp.dist = list.comp, comp.param = list.param, conf.level = 0.95,
#'                    tune.penalty = FALSE, parallel = FALSE, n_cpu = 2)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_k_samples_test <- function(samples = NULL, sim_U = NULL, n_sim_tab = 100, comp.dist = NULL, comp.param = NULL,
                               conf.level = 0.95, tune.penalty = TRUE, parallel = FALSE, n_cpu = 2)
{
  ## Control whether parallel computations were asked for or not:
  if (parallel) {
    `%fun%` <- foreach::`%dopar%`
    doParallel::registerDoParallel(cores = n_cpu)
  } else {
    `%fun%` <- foreach::`%do%`
  }

  ## Check the validity of specified arguments:
  stopifnot( length(comp.dist) == (2*length(samples)) )
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    if ( (!all(sapply(comp.param, is.null)[seq.int(from = 2, to = length(comp.dist), by = 2)] == FALSE)) |
         (!all(sapply(comp.param, is.null)[seq.int(from = 1, to = length(comp.dist), by = 2)] == TRUE)) ) {
      stop("Component distributions and/or parameters must have been badly specified in the admixture models.")
    }
  }

  ## Get the minimal size among all sample sizes, useful for contrast tabulation (adjustment of covariance terms):
  minimal_size <- min(sapply(X = samples, FUN = length))

  ## Look for all possible couples on which the discrepancy will be computed :
  model.list <- lapply(X = seq.int(from = 1, to = length(comp.dist), by = 2), FUN = seq.int, length.out = 2)
  couples.list <- NULL
  for (i in 1:(length(samples)-1)) {
    for (j in (i+1):length(samples)) { couples.list <- rbind(couples.list,c(i,j)) }
  }

  if (length(samples) == 2) {

    return(IBM_2samples_test(samples = samples, comp.dist = comp.dist, comp.param = comp.param, sim_U = sim_U, n_sim_tab = n_sim_tab,
                             min_size = minimal_size, conf.level = conf.level, parallel = parallel, n_cpu = n_cpu))

  } else {

    if (tune.penalty) {

      ##------ Calibration of tuning parameter gamma --------##
      ## Studying each population to tune parameters used in the selection rule:
      couples.expr <- couples.param <- vector(mode = "list", length = length(samples))
      ## Split each population 6 times in two subsamples:
      S_ii <- matrix(NA, nrow = length(samples), ncol = 6)
      for (k in 1:length(samples)) {
        ## Retrieves the right components for distributions and corresponding parameters under study:
        couples.expr[[k]] <- append(comp.dist[model.list[[k]]], comp.dist[model.list[[k]]])
        couples.param[[k]] <- append(comp.param[model.list[[k]]], comp.param[model.list[[k]]])
        names(couples.expr[[k]]) <- names(couples.param[[k]]) <- c("f1","g1","f2","g2")
        ## Artificially create 6 times two subsamples from one given original sample (thus under the null H0):
        subsample1_index <- matrix(NA, nrow = floor(length(samples[[k]])/2), ncol = 6)
        subsample1_index <- replicate(n = 6, expr = sample(x = 1:length(samples[[k]]), size = floor(length(samples[[k]])/2), replace = FALSE, prob = NULL))
        subsample2_index <- matrix(NA, nrow = (length(samples[[k]])-floor(length(samples[[k]])/2)), ncol = 6)
        x_val <- seq(from = min(samples[[k]]), to = max(samples[[k]]), length.out = 1000)
        ## Parallel (or not!) computations of the supremum of the difference b/w decontaminated cdf:
        l <- NULL
        S <-
          foreach::foreach(l = 1:ncol(subsample2_index), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
            subsample2_index[ ,l] <- c(1:length(samples[[k]]))[-subsample1_index[ ,l]]
            ## Comparison between the two populations :
            XY <- IBM_estimProp(sample1 = samples[[k]][subsample1_index[ ,l]], sample2 = samples[[k]][subsample2_index[ ,l]], known.prop = NULL,
                                comp.dist = couples.expr[[k]], comp.param = couples.param[[k]], with.correction = FALSE, n.integ = 1000)
            F_i1 <- decontaminated_cdf(sample1 = samples[[k]][subsample1_index[ ,l]], comp.dist = comp.dist[model.list[[k]]],
                                       comp.param = comp.param[model.list[[k]]], estim.p = XY[["prop.estim"]][1])
            F_i2 <- decontaminated_cdf(sample1 = samples[[k]][subsample2_index[ ,l]], comp.dist = comp.dist[model.list[[k]]],
                                       comp.param = comp.param[model.list[[k]]], estim.p = XY[["prop.estim"]][1])
            ## Assessment of the difference of interest at each specified x-value:
            max( sqrt(length(samples[[k]])/2) * abs(F_i1(x_val) - F_i2(x_val)) )
          }
        S_ii[k, which(sapply(X = S, FUN = class) != "numeric")] <- NA
        S_ii[k, which(sapply(X = S, FUN = class) == "numeric")] <- unlist(S[which(sapply(X = S, FUN = class) == "numeric")])
      }
      ## Extraction of the optimal gamma :
      gamma_opt <- max( apply(S_ii, 1, max, na.rm = TRUE) / log(sapply(X = samples, FUN = length)/2) )

      ##------ Calibration of tuning parameter C --------##
      couples.expr <- couples.param <- vector(mode = "list", length = length(samples))
      ## Split each population into three subsets (thus leading to be under the null):
      U_k <- vector(mode = "list", length = length(samples))
      ## Parallel (or not!) computations of the embedded test statistics for each original sample:
      U_k <-
        foreach::foreach (i = 1:length(samples), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
          couples.expr[[i]] <- append(comp.dist[model.list[[i]]], comp.dist[model.list[[i]]])
          couples.param[[i]] <- append(comp.param[model.list[[i]]], comp.param[model.list[[i]]])
          names(couples.expr[[i]]) <- names(couples.param[[i]]) <- c("f1","g1","f2","g2")
          ## Artificially create three subsamples from one given sample, thus under the null H0:
          subsample1_index <- sample(x = 1:length(samples[[i]]), size = floor(length(samples[[i]])/3), replace = FALSE, prob = NULL)
          subsample2_index <- sample(x = c(1:length(samples[[i]]))[-subsample1_index], size = length(subsample1_index), replace = FALSE, prob = NULL)
          subsample3_index <- c(1:length(samples[[i]]))[-sort(c(subsample1_index, subsample2_index))]
          ## Compute the test statistics for each combination of subsets :
          XY <- IBM_estimProp(sample1 = samples[[i]][subsample1_index], sample2 = samples[[i]][subsample2_index], known.prop = NULL,
                              comp.dist = couples.expr[[i]], comp.param = couples.param[[i]], with.correction = FALSE, n.integ = 1000)
          if (is.numeric(XY[["prop.estim"]])) {
            T_12 <- (length(samples[[i]])/3) *
              IBM_empirical_contrast(par = XY[["prop.estim"]], fixed.p.X = XY[["p.X.fixed"]], sample1 = samples[[i]][subsample1_index],
                                     sample2 = samples[[i]][subsample2_index], G = XY[["integ.supp"]],
                                     comp.dist = couples.expr[[i]], comp.param = couples.param[[i]])
          } else T_12 <- NA
          XZ <- IBM_estimProp(sample1 = samples[[i]][subsample1_index], sample2 = samples[[i]][subsample3_index], known.prop = NULL,
                              comp.dist = couples.expr[[i]], comp.param = couples.param[[i]], with.correction = FALSE, n.integ = 1000)
          if (is.numeric(XZ[["prop.estim"]])) {
            T_13 <- (length(samples[[i]])/3) *
              IBM_empirical_contrast(par = XZ[["prop.estim"]], fixed.p.X = XZ[["p.X.fixed"]], sample1 = samples[[i]][subsample1_index],
                                     sample2 = samples[[i]][subsample3_index], G = XZ[["integ.supp"]],
                                     comp.dist = couples.expr[[i]], comp.param = couples.param[[i]])
          } else T_13 <- NA
          YZ <- IBM_estimProp(sample1 = samples[[i]][subsample2_index], sample2 = samples[[i]][subsample3_index], known.prop = NULL,
                              comp.dist = couples.expr[[i]], comp.param = couples.param[[i]], with.correction = FALSE, n.integ = 1000)
          if (is.numeric(YZ[["prop.estim"]])) {
            T_23 <- (length(samples[[i]])/3) *
              IBM_empirical_contrast(par = YZ[["prop.estim"]], fixed.p.X = YZ[["p.X.fixed"]], sample1 = samples[[i]][subsample2_index],
                                     sample2 = samples[[i]][subsample3_index], G = YZ[["integ.supp"]],
                                     comp.dist = couples.expr[[i]], comp.param = couples.param[[i]])
          } else T_23 <- NA
          cumsum(x = c(T_12, T_13, T_23))
        }
      ## Remove elements of the list where an error was stored:
      U_k[which(sapply(X = U_k, FUN = function(x) any(is.na(x))) == TRUE)] <- NULL
      U_k[which(sapply(X = U_k, FUN = class) != "numeric")] <- NULL
      if (length(U_k) == 0) stop("The calibration of tuning parameters in the k-sample test procedure has failed. Please try again.")
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
    S_ij <- couples.expr <- couples.param <- vector(mode = "list", length = nrow(couples.list))
    for (k in 1:nrow(couples.list)) {
      couples.expr[[k]] <- comp.dist[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
      couples.param[[k]] <- comp.param[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
      names(couples.expr[[k]]) <- names(couples.param[[k]]) <- c("f1","g1","f2","g2")
    }
    ## Compute the supremum between decontaminated cdfs over the different couples of populations under study:
    S_ij <-
      foreach::foreach (k = 1:nrow(couples.list), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
        ## Comparison between the two populations :
        XY <- IBM_estimProp(sample1 = samples[[couples.list[k, ][1]]], sample2 = samples[[couples.list[k, ][2]]], known.prop = NULL,
                            comp.dist = couples.expr[[k]], comp.param = couples.param[[k]], with.correction = F, n.integ = 1000)
        emp.contr <- minimal_size * IBM_empirical_contrast(XY[["prop.estim"]], fixed.p.X = XY[["p.X.fixed"]],
                                                           sample1 = samples[[couples.list[k, ][1]]], sample2 = samples[[couples.list[k, ][2]]],
                                                           G = XY[["integ.supp"]], comp.dist = couples.expr[[k]], comp.param = couples.param[[k]])
        x_val <- seq(from = min(samples[[couples.list[k, ][1]]], samples[[couples.list[k, ][2]]]),
                     to = max(samples[[couples.list[k, ][1]]], samples[[couples.list[k, ][2]]]), length.out = 1000)
        F_1 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][1]]], comp.dist = couples.expr[[k]][1:2],
                                  comp.param = couples.param[[k]][1:2], estim.p = XY[["prop.estim"]][1])
        F_2 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][2]]], comp.dist = couples.expr[[k]][3:4],
                                  comp.param = couples.param[[k]][3:4], estim.p = XY[["prop.estim"]][2])
        ## Assessment of the supremum:
        list( sup = max( sqrt(minimal_size) * abs(F_1(x_val) - F_2(x_val)) ), contrast = emp.contr)
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
    if (tune.penalty) {
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
      U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 80, sample1 = samples[[which_row]], sample2 = samples[[which_col]], min_size = minimal_size,
                                     comp.dist = couples.expr[[which((couples.list[ ,1] == which_row) & (couples.list[ ,2] == which_col))]],
                                     comp.param = couples.param[[which((couples.list[ ,1] == which_row) & (couples.list[ ,2] == which_col))]],
                                     parallel = parallel, n_cpu = n_cpu)
      sim_U <- U[["U_sim"]]
    } else {
      sim_U <- sim_U
    }
    if (all(is.na(sim_U))) {
      return( list(rejection_rule = TRUE, p_val = NA, stat_name = NA, stat_value = NA, stat_rank = NA,
                   penalized_stat = NA, ordered_contrasts = NA, quantiles = NA, stat_terms = NA,
                   contr_mat = contrast.matrix) )
    } else {
      q_H <- stats::quantile(sim_U, conf.level)
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

    return( list(confidence_level = conf.level, rejection_rule = final_test, p_value = p_value, stat_name = finalStat_name,
                 test.stat = finalStat_value, stat_rank = selected.index, penalized_stat = penalized.stats, penalty_H0 = penalty_rule,
                 gamma = gamma_opt, tuning_constant = cst_selected, ordered_contrasts = sorted_contrasts, quantiles = unique(q_H),
                 sim_U = sim_U, stat_terms = stat.terms_final, contr_mat = contrast.matrix) )
  }

}
