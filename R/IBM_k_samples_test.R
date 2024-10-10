#' Equality test of K unknown component distributions
#'
#' Equality test of the unknown component distributions coming from K (K > 1) admixture models, based on the Inversion - Best
#' Matching (IBM) approach. Recall that we have K populations following admixture models, each one with probability
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
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024b}{admix}
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
#            F_i1 <- decontaminated_cdf(sample1 = samples[[k]][subsample1_index[ ,l]], comp.dist = comp.dist[model.list[[k]]],
#                                       comp.param = comp.param[model.list[[k]]], estim.p = XY[["p.X.fixed"]])
            F_i1 <- decontaminated_cdf(sample1 = samples[[k]][subsample1_index[ ,l]], comp.dist = comp.dist[model.list[[k]]],
                                       comp.param = comp.param[model.list[[k]]], estim.p = XY[["prop.estim"]])
            F_i2 <- decontaminated_cdf(sample1 = samples[[k]][subsample2_index[ ,l]], comp.dist = comp.dist[model.list[[k]]],
                                       comp.param = comp.param[model.list[[k]]], estim.p = XY[["prop.estim"]])
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
        if (length(XY[["prop.estim"]]) == 2) {
          F_1 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][1]]], comp.dist = couples.expr[[k]][1:2],
                                    comp.param = couples.param[[k]][1:2], estim.p = XY[["prop.estim"]][1])
          F_2 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][2]]], comp.dist = couples.expr[[k]][3:4],
                                    comp.param = couples.param[[k]][3:4], estim.p = XY[["prop.estim"]][2])
#          f1 <- decontaminated_density(sample1 = samples[[couples.list[k, ][1]]], comp.dist = couples.expr[[k]][1:2],
#                                       comp.param = couples.param[[k]][1:2], estim.p = XY[["prop.estim"]][1])
#          f2 <- decontaminated_density(sample1 = samples[[couples.list[k, ][2]]], comp.dist = couples.expr[[k]][3:4],
#                                       comp.param = couples.param[[k]][3:4], estim.p = XY[["prop.estim"]][2])
#          plot(f1, x_val = 0:120, type = "l")
#          plot(f2, x_val = 0:120, col = "red", type = "l", add_plot = TRUE)
        } else {
          F_1 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][1]]], comp.dist = couples.expr[[k]][1:2],
                                    comp.param = couples.param[[k]][1:2], estim.p = XY[["p.X.fixed"]])
          F_2 <- decontaminated_cdf(sample1 = samples[[couples.list[k, ][2]]], comp.dist = couples.expr[[k]][3:4],
                                    comp.param = couples.param[[k]][3:4], estim.p = XY[["prop.estim"]])
        }
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


#' Test equality of 2 unknown component distributions using IBM
#'
#' Two-sample equality test of the unknown component distributions coming from two admixture models, using Inversion - Best Matching
#' (IBM) approach. Recall that we have two admixture models with respective probability density functions (pdf)
#' l1 = p1 f1 + (1-p1) g1 and l2 = p2 f2 + (1-p2) g2, where g1 and g2 are known pdf and l1 and l2 are observed.
#' Perform the following hypothesis test: H0 : f1 = f2 versus H1 : f1 differs from f2.
#'
#' @param samples A list of the two observed samples, where each sample follows the mixture distribution given by l = p*f + (1-p)*g,
#'                with f and p unknown and g known.
#' @param known.p (default to NULL) Numeric vector with two elements, the known (true) mixture weights.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param sim_U Random draws of the inner convergence part of the contrast as defined in the IBM approach (see 'Details' below).
#' @param n_sim_tab Number of simulated gaussian processes used in the tabulation of the inner convergence distribution in the IBM approach.
#' @param min_size (default to NULL, only used with 'ICV' testing method in the k-sample case, otherwise useless) Minimal size among all samples (needed
#'                  to take into account the correction factor for the variance-covariance assessment).
#' @param conf.level The confidence level of the 2-samples test, i.e. the quantile level to which the test statistic is compared.
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when parallelizing.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return A list of five elements, containing : 1) the test statistic value; 2) the rejection decision; 3) the p-value of the
#'         test, 4) the estimated weights of the unknown component for each of the two admixture models, 5) the simulated distribution
#'         of the inner convergence regime (useful to perform the test when comparing to the extreme quantile of this distribution).
#'
#' @examples
#' \donttest{
#' ####### Under the null hypothesis H0 :
#' ## Simulate data:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm")
#' list.param <- list(f1 = list(mean = 1, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 1, sd = 1), g2 = list(mean = 3, sd = 1.2))
#' X.sim <- rsimmix(n= 1100, unknownComp_weight=0.85, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' Y.sim <- rsimmix(n= 1200, unknownComp_weight=0.75, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 3, sd = 1.2))
#' IBM_2samples_test(samples = list(X.sim, Y.sim), known.p = NULL, comp.dist = list.comp,
#'                   comp.param = list.param, sim_U = NULL, n_sim_tab = 6, min_size = NULL,
#'                   conf.level = 0.95, parallel = FALSE, n_cpu = 2)
#'
#' ####### Under the alternative H1 :
#' ## Simulate data:
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm")
#' list.param <- list(f1 = list(mean = 1, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 2, sd = 1), g2 = list(mean = 3, sd = 1.2))
#' X.sim <- rsimmix(n= 1100, unknownComp_weight=0.85, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' Y.sim <- rsimmix(n= 1200, unknownComp_weight=0.75, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 3, sd = 1.2))
#' IBM_2samples_test(samples = list(X.sim, Y.sim), known.p = NULL, comp.dist = list.comp,
#'                   comp.param = list.param, sim_U = NULL, n_sim_tab = 6, min_size = NULL,
#'                   conf.level = 0.95, parallel = FALSE, n_cpu = 2)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_2samples_test <- function(samples, known.p = NULL, comp.dist = NULL, comp.param = NULL, sim_U = NULL,
                              n_sim_tab = 50, min_size = NULL, conf.level = 0.95, parallel = FALSE, n_cpu = 2)
{
  ## Estimate the proportions of the mixtures:
  estim <- try(suppressWarnings(IBM_estimProp(sample1 = samples[[1]], sample2 = samples[[2]], known.prop = known.p, comp.dist = comp.dist,
                                              comp.param = comp.param, with.correction = FALSE, n.integ = 1000)),
               silent = TRUE)
  count_error <- 0
  while ((inherits(x = estim, what = "try-error", which = FALSE)) & (count_error < 3)) {
    estim <- NULL
    estim <- try(suppressWarnings(IBM_estimProp(sample1 = samples[[1]], sample2 = samples[[2]], known.prop = known.p, comp.dist = comp.dist,
                                                comp.param = comp.param, with.correction = FALSE, n.integ = 1000)),
                 silent = TRUE)
    count_error <- count_error + 1
  }
  if (inherits(x = estim, what = "try-error", which = FALSE)) {
    estim.weights <- 100
    contrast_val <- NA
  } else {
    #  if (length(estim[["prop.estim"]] == 2)) {
    estim.weights <- estim[["prop.estim"]]
    #  } else {
    #    estim.weights <- c(estim[["p.X.fixed"]], estim[["prop.estim"]])
    #  }
    if (is.null(min_size)) {
      sample.size <- min(length(samples[[1]]), length(samples[[2]]))
    } else {
      sample.size <- min_size
    }
    contrast_val <- sample.size * IBM_empirical_contrast(par = estim.weights, fixed.p.X = estim[["p.X.fixed"]], sample1 = samples[[1]],
                                                         sample2 = samples[[2]], G = estim[["integ.supp"]], comp.dist = comp.dist, comp.param = comp.param)
    ## Identify boolean 'green light' criterion to known whether we need to perform the test with stochastic integral tabulation:
    # green_light <- IBM_greenLight_criterion(estim.obj = estim, sample1 = samples[[1]], sample2 = samples[[2]], comp.dist = comp.dist,
    #                                        comp.param = comp.param, min_size = min_size, alpha = 0.05)
    # if (green_light) {
    #   if (is.null(sim_U)) {
    #     tab_dist <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 50, sample1 = samples[[1]], sample2 = samples[[2]], min_size = min_size,
    #                                          comp.dist = comp.dist, comp.param = comp.param, parallel = parallel, n_cpu = n_cpu)
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
    reject.decision <- TRUE
    p_value <- 1e-16
    sim_U <- NA
  } else {
    if (is.null(sim_U)) {
      tab_dist <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 60, sample1 = samples[[1]], sample2 = samples[[2]], min_size = min_size,
                                            comp.dist = comp.dist, comp.param = comp.param, parallel = parallel, n_cpu = n_cpu)
      sim_U <- tab_dist$U_sim
      contrast_val <- tab_dist$contrast_value
    }
    reject.decision <- contrast_val > stats::quantile(sim_U, conf.level)
    CDF_U <- stats::ecdf(sim_U)
    p_value <- 1 - CDF_U(contrast_val)
  }

  return( list(confidence_level = conf.level, rejection_rule = reject.decision, p_value = p_value,
               test.stat = contrast_val, weights = estim.weights, sim_U = sim_U) )
}


#' Green-light criterion for equality tests in admixture models
#'
#' Indicates whether one needs to perform the full version of the statistical test of equality between unknown component
#' distributions in two samples following admixture models (it is sometimes obvious that they differ), based on the IBM approach.
#' See more in 'Details' below.
#'
#' @param estim.obj Object obtained from the estimation of the component weights related to the proportions of the unknown component
#'                  in each of the two admixture models studied.
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param min_size (optional, NULL by default) In the k-sample case, useful to provide the minimal size among all samples
#'                 (needed to take into account the correction factor for variance-covariance assessment). Otherwise, useless.
#' @param alpha Confidence level at which the criterion is assessed (used to compute the confidence bands of the estimators
#'              of the unknown component weights).
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return A boolean indicating whether it is useful or useless to tabulate the contrast distribution in order to answer
#'         the testing problem (f1 = f2).
#'
#' @examples
#' \donttest{
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=550, unknownComp_weight=0.7, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=450, unknownComp_weight=0.8, comp.dist = list(list.comp$f2,list.comp$g2),
#'                                                    comp.param=list(list.param$f2,list.param$g2))
#' ## Estimate the unknown component weights in the two admixture models in real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
#'                    f2 = NULL, g2 = list(mean = 5, sd = 2))
#' estim <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], known.prop = NULL,
#'                        comp.dist = list.comp, comp.param = list.param,
#'                        with.correction = FALSE, n.integ = 1000)
#' IBM_greenLight_criterion(estim.obj = estim, sample1 = sample1[['mixt.data']],
#'                         sample2 = sample2[['mixt.data']], comp.dist = list.comp,
#'                         comp.param = list.param, min_size = NULL, alpha = 0.05)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

IBM_greenLight_criterion <- function(estim.obj, sample1, sample2, comp.dist = NULL, comp.param = NULL, min_size = NULL, alpha = 0.05)
{
  if (is.null(min_size)) {
    min_sample_size <- min(length(sample1), length(sample2))
  } else {
    min_sample_size <- min_size
  }

  length.support <- length(estim.obj[["integ.supp"]])
  z <- estim.obj[["integ.supp"]][round(floor(length.support/2))]
  varCov_estim <- IBM_estimVarCov_gaussVect(x = z, y = z, estim.obj = estim.obj, fixed.p1 = estim.obj[["p.X.fixed"]], known.p = NULL,
                                            sample1 = sample1, sample2 = sample2, min_size = min_size,
                                            comp.dist = comp.dist, comp.param = comp.param)

  if (length(estim.obj[["prop.estim"]]) == 2) {

    inf_bound.p1 <- estim.obj[["prop.estim"]][1] - sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p1 <- estim.obj[["prop.estim"]][1] + sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p1 <- c(inf_bound.p1, sup_bound.p1)
    inf_bound.p2 <- estim.obj[["prop.estim"]][2] - sqrt(varCov_estim[2,2]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p2 <- estim.obj[["prop.estim"]][2] + sqrt(varCov_estim[2,2]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p2 <- c(inf_bound.p2, sup_bound.p2)
    green_light_crit <- max(conf_interval.p1[1],conf_interval.p2[1]) <= 1

  } else {

    inf_bound.p <- estim.obj[["prop.estim"]][1] - sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    sup_bound.p <- estim.obj[["prop.estim"]][1] + sqrt(varCov_estim[1,1]/min_sample_size) * stats::qnorm(p=(1-alpha/4), mean=0, sd=1)
    conf_interval.p2 <- c(inf_bound.p, sup_bound.p)
    conf_interval.p1 <- NULL
    green_light_crit <- conf_interval.p2[1] <= 1

  }

  return( list(green_light = green_light_crit, conf_interval_p1 = conf_interval.p1, conf_interval_p2 = conf_interval.p2) )
}


#' Simulated distribution of the contrast using IBM
#'
#' Tabulate the distribution related to the inner convergence part of the contrast, by simulating trajectories of Gaussian
#' processes and applying some transformations. Useful to perform the test hypothesis, by retrieving the (1-alpha)-quantile
#' of interest. See 'Details' below and the cited paper therein for further information.
#'
#' @param n.sim Number of trajectories of simulated gaussian processes (number of random draws for tabulation).
#' @param n.varCovMat Number of time points on which gaussian processes are simulated.
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param min_size (default to NULL) In the k-sample case, useful to provide the minimal size among all samples. Otherwise, useless.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used for computations when parallelizing.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024a}{admix}
#'
#' @return A list with four elements, containing: 1) random draws of the quantity 'sample size times the empirical contrast',
#'         as defined in the IBM approach (see 'Details' above); 2) the estimated unknown component weights for the two admixture
#'         models; 3) the value of the quantity 'sample size times the empirical contrast'; 4) the sequence of points in the
#'         support that were used to evaluate the variance-covariance matrix of empirical processes.
#'
#' @examples
#' \donttest{
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 1, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 1, sd = 1), g2 = list(mean = 3, sd = 1.2))
#' X.sim <- rsimmix(n=1000, unknownComp_weight=0.7, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' Y.sim <- rsimmix(n=1200, unknownComp_weight=0.6, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' ## Tabulate 1st term of stochastic integral (inner convergence) in a real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 3, sd = 1.2))
#' U <- IBM_tabul_stochasticInteg(n.sim = 2, n.varCovMat = 20, sample1 = X.sim, sample2 = Y.sim,
#'                                min_size = NULL, comp.dist = list.comp, comp.param = list.param,
#'                                parallel = FALSE, n_cpu = 2)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

IBM_tabul_stochasticInteg <- function(n.sim = 200, n.varCovMat = 100, sample1 = NULL, sample2 = NULL, min_size = NULL,
                                      comp.dist = NULL, comp.param = NULL, parallel = FALSE, n_cpu = 2)
{
  i <- NULL
  if (parallel) {
    `%fun%` <- foreach::`%dopar%`
    doParallel::registerDoParallel(cores = n_cpu)
  } else {
    `%fun%` <- foreach::`%do%`
  }

  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
    stop("Known components must be specified in the two admixture models.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
  }

  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
  #  comp_tab <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_tab <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_tab)) assign(x = names(comp_tab)[i], value = comp_tab[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_tab)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                                                                                        paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_tab)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Differentiates when the 4 components were specified (simulations) and real-life cases (only known components are given).
  ## First, real-life case:
  if (all(unlist(sapply(comp.param, is.na)[c(1,3)]))) {
    if (any(exp.comp.dist == "stepfun")) {
      G1.fun <- eval(parse(text = expr[2]))
      G2.fun <- eval(parse(text = expr[4]))
      G1 <- function(z) G1.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      G1 <- function(z) { eval(parse(text = expr[2])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  } else {
    ## case of simulated data:
    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      G1.fun <- eval(parse(text = expr[2]))
      F2.fun <- eval(parse(text = expr[3]))
      G2.fun <- eval(parse(text = expr[4]))
      F1 <- function(z) F1.fun(z)
      G1 <- function(z) G1.fun(z)
      F2 <- function(z) F2.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
      F2 <- function(z) { eval(parse(text = expr[3])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  }
  ## Empirical cumulative distribution function from the two observed samples:
  L1 <- stats::ecdf(sample1)
  L2 <- stats::ecdf(sample2)

  ## Estimate the weights of the unknown component distributions in first and second samples:
  estim <- IBM_estimProp(sample1 = sample1, sample2 = sample2, known.prop = NULL, comp.dist = comp.dist, comp.param = comp.param,
                         with.correction = F, n.integ = 1000)
  if (is.null(min_size)) {
    sample.size <- min(length(sample1), length(sample2))
  } else {
    sample.size <- min_size
  }
  contrast_val <- sample.size *
    IBM_empirical_contrast(par = estim[["prop.estim"]], fixed.p.X = estim[["p.X.fixed"]], sample1 = sample1,
                           sample2 = sample2, G = estim[["integ.supp"]], comp.dist = comp.dist, comp.param = comp.param)
  ## Integration support:
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    densite.G <- stats::density(estim[["integ.supp"]], bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integ <- c(min(estim[["integ.supp"]]), max(estim[["integ.supp"]]))
    t_seq <- seq(from = supp.integ[1], to = supp.integ[2], length.out = n.varCovMat)
  } else {
    ## Case of multinomial distribution :
    supp.integ <- estim[["integ.supp"]]
    t_seq <- sort(unique(supp.integ))
  }

  ## Compute the normalization matrix M(.) at each point, to be used further when determining the full simulated trajectories:
  normalization_factors <-
    foreach::foreach (i = 1:length(t_seq), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
      IBM_normalization_term(t_seq[i], estim, estim[["p.X.fixed"]], NULL, sample1, sample2, min_size, comp.dist, comp.param)
    }

  ## Estimate the variance-covariance functions from the empirical processes:
  #  cov_mat_L1 <- sapply(t_seq, function(s1) {
  #                     sapply(t_seq, function(s2) {
  #                       estimVarCov_empProcess(x = s1, y = s2, obs.data = sample1) })
  #                     })
  cov_mat_L1 <- estimVarCov_empProcess_Rcpp(t = t_seq, obs_data = sample1)
  #  cov_mat_L2 <- sapply(t_seq, function(s1) {
  #                    sapply(t_seq, function(s2) {
  #                      estimVarCov_empProcess(x = s1, y = s2, obs.data = sample2) })
  #  })
  cov_mat_L2 <- estimVarCov_empProcess_Rcpp(t = t_seq, obs_data = sample2)

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)
  if (G1equalG2) {
    psi1 <- function(z) 2*( ((2-estim[["prop.estim"]])/estim[["prop.estim"]]^3) * G1(z) -
                              (2/estim[["prop.estim"]]^3)*L2(z) + (1/(estim[["prop.estim"]]^2*estim[["p.X.fixed"]]))*L1(z) -
                              ((1-estim[["p.X.fixed"]])/(estim[["prop.estim"]]^2*estim[["p.X.fixed"]])) * G1(z) )
    psi2 <- function(z) 2*( (1/(estim[["prop.estim"]]^2*estim[["p.X.fixed"]])) * (L2(z) - G1(z)) )
  } else {
    psi1.1 <- function(z) 2*( ((2-estim[["prop.estim"]][1])/estim[["prop.estim"]][1]^3) * G1(z) -
                                (2/estim[["prop.estim"]][1]^3)*L1(z) + (1/(estim[["prop.estim"]][1]^2*estim[["prop.estim"]][2]))*L2(z) -
                                ((1-estim[["prop.estim"]][2])/(estim[["prop.estim"]][1]^2*estim[["prop.estim"]][2])) * G2(z) )
    psi1.2 <- function(z) 2*( (1/(estim[["prop.estim"]][1]^2*estim[["prop.estim"]][2])) * (L1(z) - G1(z)) )
    psi2.1 <- function(z) 2*( ((2-estim[["prop.estim"]][2])/estim[["prop.estim"]][2]^3) * G2(z) -
                                (2/estim[["prop.estim"]][2]^3)*L2(z) + (1/(estim[["prop.estim"]][2]^2*estim[["prop.estim"]][1]))*L1(z) -
                                ((1-estim[["prop.estim"]][1])/(estim[["prop.estim"]][2]^2*estim[["prop.estim"]][1])) * G1(z) )
    psi2.2 <- function(z) 2*( (1/(estim[["prop.estim"]][2]^2*estim[["prop.estim"]][1])) * (L2(z) - G2(z)) )
  }

  U_sim <-
    foreach::foreach (i = 1:n.sim, .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
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
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi2(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi1 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi1(z) * B2_traj(z) * densite.G.dataPoint
        }
        Z <- matrix(NA, nrow = 4, ncol = length(t_seq))
        if (support == "continuous") {
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
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi1.1(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi2.2 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi2.2(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi2.1 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi2.1(z) * B2_traj(z) * densite.G.dataPoint
        }
        integrand.phi1.2 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi1.2(z) * B2_traj(z) * densite.G.dataPoint
        }
        Z <- matrix(NA, nrow = 6, ncol = length(t_seq))
        if (support == "continuous") {
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
        if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
        } else {
          densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
        }
        D_function(z)^2 * densite.G.dataPoint
      }

      if (support == "continuous") {
        U <- stats::integrate(integrand, lower = supp.integ[1], upper = supp.integ[2], subdivisions = 10000L, rel.tol = 1e-04)$value
      } else {
        U <- sum( unlist(sapply(supp.integ, integrand)) )
      }
      U
    }

  indexes.toRemove <- which( (substr(U_sim, start = 1, stop = 5) == "Error") | (substr(U_sim, start = 1, stop = 5) == "list(") )
  if (length(indexes.toRemove) != 0) U_sim <- U_sim[-indexes.toRemove]

  return( list(U_sim = as.numeric(U_sim), estimator = estim, contrast_value = contrast_val, integ.points = t_seq) )

}


#' Covariance matrix of some empirical process
#'
#' Estimate the variance-covariance matrix of some given empirical process, based on the Donsker correlation.
#' Compute Donsker correlation between two time points (x,y) for some given empirical process with R code
#' (another implementation in C++ is also available to speed up this computation).
#'
#' @param x First time point considered for the computation of the correlation given the empirical process.
#' @param y Second time point considered for the computation of the correlation given the same empirical process.
#' @param obs.data Sample that permits to estimate the cumulative distribution function (cdf).
#' @param known.p NULL by default (only useful to compute the exact Donsker correlation). The component weight dedicated to
#'                the unknown mixture component if available (in case of simulation studies)
#' @param comp.dist NULL by default (only useful to compute the exact Donsker correlation). Otherwise, a list with two elements
#'                  corresponding to component distributions (specified with R native names for these distributions) involved
#'                  in the admixture model. All elements must be specified, for instance list(f='norm', g='norm').
#' @param comp.param NULL by default (only useful to compute the exact Donsker correlation). Otherwise, a list with two elements
#'                   corresponding to the parameters of the component distributions, each element being a list itself. The names
#'                   used in this list must correspond to the native R argument names for these distributions.
#'                   All elements must be specified, for instance list(f=NULL, g=list(mean=0,sd=1)).
#'
#' @return The estimated variance-covariance matrix.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm')
#' list.param <- list(f1 = list(mean = 12, sd = 0.4),
#'                    g1 = list(mean = 16, sd = 0.7))
#' obs.data <- rsimmix(n=2500, unknownComp_weight=0.5, comp.dist=list.comp, comp.param= list.param)
#' ## Compute the variance-covariance matrix of the corresponding empirical process:
#' t <- seq(from = min(obs.data$mixt.data), to = max(obs.data$mixt.data), length = 50)
#' S2 <- sapply(t, function(s1) {
#'                 sapply(t, function(s2) {
#'                      estimVarCov_empProcess(x = s1, y = s2, obs.data = obs.data$mixt.data) })
#'                 })
#' lattice::wireframe(S2)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

estimVarCov_empProcess <- function(x, y, obs.data, known.p = NULL, comp.dist = NULL, comp.param = NULL)
{
  if (is.null(obs.data)) {
    stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
    if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) | is.null(known.p)) {
      stop("All parameters of the admixture model must be specified to compute the exact Donsker correlation.")
    }
    ## Extracts the information on component distributions:
    exp.comp.dist <- paste0("p", comp.dist)
    if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
    #    comp.ro <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
    comp.ro <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
    for (i in 1:length(comp.ro)) assign(x = names(comp.ro)[i], value = comp.ro[[i]])
    ## Creates the expression involved in future assessments of the CDF:
    make.expr.step <- function(i) paste(names(comp.ro)[i], "(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                                                                                                          paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
    make.expr <- function(i) paste(names(comp.ro)[i], "(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
    expr <- vector(mode = "character", length = length(exp.comp.dist))
    expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
    expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
    expr <- unlist(expr)

    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      F1 <- function(z) F1.fun(z)
      G1.fun <- eval(parse(text = expr[2]))
      G1 <- function(z) G1.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
    }

    L.CDF <- function(z) { known.p * F1(z) + (1-known.p) * G1(z) }

  } else {

    L.CDF <- stats::ecdf(obs.data)
  }

  ## Use the Donsker correlation formula:
  res <- L.CDF(min(x,y)) * (1 - L.CDF(max(x,y)))

  return(res)
}
