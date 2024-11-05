#' Clustering of K populations following admixture models
#'
#' Create clusters on the unknown components related to the K populations following admixture models. Based on the K-sample test
#' using Inversion - Best Matching (IBM) approach, see 'Details' below for further information.
#'
#' @param samples A list of the K (K>1) samples to be studied, all following admixture distributions.
#' @param admixMod A list of objects of class 'admix_model', containing useful information about distributions and parameters.
#' @param conf_level (default to 0.95) The confidence level of the k-sample tests used in the clustering procedure.
#' @param n_sim_tab (default to 100) Number of simulated Gaussian processes when tabulating the inner convergence distribution
#'                  in the IBM approach.
#' @param tune_penalty (default to FALSE) A boolean that allows to choose between a classical penalty term or an optimized penalty (embedding
#'                     some tuning parameters, automatically optimized). Optimized penalty is particularly useful for low or unbalanced sample sizes
#'                     to detect alternatives to the null hypothesis (H0).
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when paralleling computations.
#' @param tabul_dist (default to NULL) Only useful for comparisons of detected clusters at different confidence levels. A list of
#'                    the tabulated distributions of the stochastic integral used in the k-sample test, each element for each
#'                    cluster previously detected.
#' @param echo (default to FALSE) Display the remaining computation time.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024b}{admix}
#'
#' @return An object of class 'admix_cluster', containing 12 attributes: 1) the number of samples under study; 2) the sizes of samples;
#'         3) the information about mixture components in each sample (distributions and parameters); 4) the number of detected clusters;
#'         5) the list of p-values for each k-sample test at the origin of detected clusters; 6) the cluster affiliation for each sample;
#'         7) the confidence level of statistical tests; 8) which samples in which cluster; 9) the size of clusters; 10) the estimated
#'         weights of the unknown component distributions inside each cluster (remind that estimated weights are consistent only if
#'         unknown components are tested to be identical, which is the case inside clusters); 11) the matrix of pairwise discrepancies
#'         across all samples; 12) the list of tabulated distributions used for statistical tests involved in building the clusters.
#'
#' @examples
#' \donttest{
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 2600, weight = 0.8,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 16, "scale" = 1/4),
#'                                         list("rate" = 1/3.5)))
#' mixt2 <- twoComp_mixt(n = 3000, weight = 0.7,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 14, "scale" = 1/2),
#'                                         list("rate" = 1/5)))
#' mixt3 <- twoComp_mixt(n = 3500, weight = 0.6,
#'                       comp.dist = list("gamma", "gamma"),
#'                       comp.param = list(list("shape" = 16, "scale" = 1/4),
#'                                         list("shape" = 12, "scale" = 1/2)))
#' mixt4 <- twoComp_mixt(n = 4800, weight = 0.5,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 14, "scale" = 1/2),
#'                                         list("rate" = 1/7)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' data3 <- getmixtData(mixt3)
#' data4 <- getmixtData(mixt4)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
#'                          knownComp_param = mixt3$comp.param[[2]])
#' admixMod4 <- admix_model(knownComp_dist = mixt4$comp.dist[[2]],
#'                          knownComp_param = mixt4$comp.param[[2]])
#' ## Clustering procedure:
#' admix_cluster(samples = list(data1, data2, data3, data4),
#'               admixMod = list(admixMod1, admixMod2, admixMod3, admixMod4),
#'               conf_level = 0.95, n_sim_tab = 30, tune_penalty = TRUE,
#'               parallel = FALSE, n_cpu = 2, tabul_dist = NULL, echo = FALSE)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_cluster <- function(samples, admixMod, conf_level = 0.95, n_sim_tab = 30, tune_penalty = FALSE,
                          parallel = FALSE, n_cpu = 2, tabul_dist = NULL, echo = FALSE)
{
  old_options_warn <- base::options()$warn
  base::options(warn = -1)

  if (length(sapply(samples, length)) == 1) return("One single sample, no clusters to be found.")

  ## Control whether parallel computations were asked for or not:
  if (parallel) {
    `%fun%` <- doRNG::`%dorng%`
    doParallel::registerDoParallel(cores = n_cpu)
  } else {
    `%fun%` <- foreach::`%do%`
  }

  ## Get the minimal size among all sample sizes, useful for future tabulation (adjustment of variance-covariance):
  minimal_size <- min(sapply(X = samples, FUN = length))
  index_minimal_size <- which.min(sapply(X = samples, FUN = length))
  ## K*(K-1)/2 combinations of populations under study:
  couples.list <- NULL
  for (i in 1:(length(samples)-1)) {
    for (j in (i+1):length(samples)) {
      couples.list <- rbind(couples.list,c(i,j))
    }
  }
  couples.expr <- couples.param <- vector(mode = "list", length = nrow(couples.list))
  empirical.contr <-
  foreach::foreach (k = 1:nrow(couples.list), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
    ## Comparison between the two populations :
    XY <- estim_IBM(samples = samples[as.numeric(couples.list[k, ])],
                    admixMod = admixMod[as.numeric(couples.list[k, ])], n.integ = 1000)
    while (!exists("XY")) {
      XY <- estim_IBM(samples = samples[as.numeric(couples.list[k, ])],
                      admixMod = admixMod[as.numeric(couples.list[k, ])], n.integ = 1000)
    }
    minimal_size * IBM_empirical_contrast(XY$estimated_mixing_weights, samples = samples[as.numeric(couples.list[k, ])],
                                          admixMod = admixMod[as.numeric(couples.list[k, ])], G = XY$integ.supp, fixed.p.X = XY$p.X.fixed)
  }

  ## Manage cases when the optimization algorithm could not find a solution:
  if (length(which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE)) != 0) {
    for (k in 1:length(which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE))) {
      empirical.contr[[which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE)[k]]] <- NA
    }
  }

  ## Give a simple and useful representation of results for each couple:
  H0.test <- contrast.matrix <- discrepancy.id <- discrepancy.rank <- matrix(NA, nrow = length(samples), ncol = length(samples))
  ## Store results of pairwise tests for equality between unknown components (fasten the future research of potential clusters):
  weights.list <- matrix(data = NA, nrow = nrow(couples.list), ncol = 2)
  for (k in 1:nrow(couples.list)) {
    contrast.matrix[couples.list[k, ][1], couples.list[k, ][2]] <- empirical.contr[[k]]
    discrepancy.id[couples.list[k, ][1], couples.list[k, ][2]] <- paste(couples.list[k, ][1], couples.list[k, ][2], sep = "-")
    discrepancy.rank[couples.list[k, ][1], couples.list[k, ][2]] <- rank(unlist(empirical.contr))[k]
    ## Pairwise testing: test H0 between the two considered populations.
    pairwise_H0test <- NULL
    pairwise_H0test <- IBM_2samples_test(samples = samples[as.numeric(couples.list[k, ])],
                                         admixMod = admixMod[as.numeric(couples.list[k, ])],
                                         conf_level = conf_level, parallel = parallel, n_cpu = n_cpu,
                                         n_sim_tab = n_sim_tab)
    weights.list[k, ] <- pairwise_H0test$estimated_mixing_weights
    H0.test[couples.list[k, ][1], couples.list[k, ][2]] <- pairwise_H0test$reject_decision
  }

  ## Start algorithm to create the clusters :
  clusters <- tab_distrib <- vector(mode = "list")
  which_row <- unlist(apply(discrepancy.rank, 2, function(x) which(x == 1)))
  which_col <- unlist(apply(discrepancy.rank, 1, function(x) which(x == 1)))
  neighboors_index_init <- which.min(contrast.matrix)
  first.group <- discrepancy.id[which_row, which_col]
  if (is.null(tabul_dist)) {
    U <- IBM_tabul_stochasticInteg(samples = list(samples[[which_row]], samples[[which_col]]),
                                   admixMod = list(admixMod[[which_row]], admixMod[[which_col]]),
                                   min_size = minimal_size, parallel = parallel, n_cpu = n_cpu,
                                   n.varCovMat = 80, n_sim_tab = n_sim_tab)
    Usim <- tab_distrib[[1]] <- U$U_sim
    CDF_U <- stats::ecdf(U$U_sim)
    q_H <- stats::quantile(U$U_sim, conf_level)
  } else {
    Usim <- tabul_dist[[1]]
    tab_distrib[[1]] <- tabul_dist[[1]]
    CDF_U <- stats::ecdf(tabul_dist[[1]])
    q_H <- stats::quantile(tabul_dist[[1]], conf_level)
  }
  test.H0 <- contrast.matrix[which_row, which_col] > q_H
  ## Numerical vector storing the p-values associated to each test: each cluster is closed once the p-value lies below the H0-rejection
  ## threshold (1-conf_level). This enables to see whether we were close to close the cluster or not each time we add a new population.
  p_value <- numeric(length = 0L)
  if (!test.H0) {
    clusters[[1]] <- alreadyGrouped_samples <- as.character(c(which_row, which_col))
    neighboors <- c(discrepancy.id[which_row, ], discrepancy.id[which_col, ], discrepancy.id[ ,which_row], discrepancy.id[ ,which_col])
    new.n_clust <- n_clust <- 1
    ## Look for a another member to integrate the newly built cluster:
    indexesSamples_to_consider <- c(0,0,0)
    p_value <- c(p_value, 1 - CDF_U(contrast.matrix[which_row, which_col]))
    CDF_U <- NULL
    if (echo) {
      prog_bar <- utils::txtProgressBar(min = 0, max = length(samples), style = 3, width = 50, char = "=")
      utils::setTxtProgressBar(prog_bar, 2/length(samples))
    }
  } else {
    clusters[[1]] <- as.character(which_row)
    clusters[[2]] <- as.character(which_col)
    alreadyGrouped_samples <- c(as.character(which_row), as.character(which_col))
    neighboors <- c(discrepancy.id[which_col, ], discrepancy.id[ ,which_col])
    n_clust <- 1
    new.n_clust <- 2
    ## Look for a second member to integrate the newly second built cluster:
    indexesSamples_to_consider <- c(0,0)
    if (echo) {
      prog_bar <- utils::txtProgressBar(min = 0, max = length(samples), style = 3,  width = 50, char = "=")
      utils::setTxtProgressBar(prog_bar, 1/length(samples))
    }
  }

  while (length(alreadyGrouped_samples) < length(samples)) {
    ## Detect which are the neighboors of the current population under study:
    neighboors_index <- sort(unique(match(x = neighboors[-which(is.na(neighboors))], table = c(discrepancy.id))))
    ## Remove already studied couples:
    neighboors_index <- neighboors_index[-which(!is.na(match(x = neighboors_index, table = neighboors_index_init)))]
    nearest_neighboor_index <- neighboors_index[which.min(contrast.matrix[neighboors_index])]
    neighboors_index_init <- c(neighboors_index_init, nearest_neighboor_index)
    first_sample <- strsplit(x = discrepancy.id[nearest_neighboor_index], split="-")[[1]][1]
    second_sample <- strsplit(x = discrepancy.id[nearest_neighboor_index], split="-")[[1]][2]
    if (all(c(first_sample,second_sample) %in% alreadyGrouped_samples)) {
      #print("Already affiliated to one existing cluster")
      which_row <- which_col <- as.numeric(clusters[[length(clusters)]])
      neighboors <- c(discrepancy.id[which_row, ], discrepancy.id[ ,which_col])
    } else {
      indexesSamples_to_consider_new <- sort(as.numeric(unique(c(first_sample, second_sample, clusters[[length(clusters)]]))))
      ## Tabulate the new inner convergence distribution:
      if (new.n_clust == (n_clust+1)) {
        n_clust <- n_clust + 1
        if ( is.null(tabul_dist) | (new.n_clust > length(tabul_dist)) ) {
          U <- IBM_tabul_stochasticInteg(samples = list(samples[[as.numeric(first_sample)]], samples[[as.numeric(second_sample)]]),
                                         admixMod = list(admixMod[[as.numeric(first_sample)]], admixMod[[as.numeric(second_sample)]]),
                                         min_size = minimal_size, parallel = parallel, n_cpu = n_cpu, n.varCovMat = 80, n_sim_tab = n_sim_tab)
          ## Store the simulated tabulated distribution:
          Usim <- U$U_sim
          tab_distrib <- append(tab_distrib, list(Usim))
          q_H <- stats::quantile(U$U_sim, conf_level)
          CDF_U <- stats::ecdf(U$U_sim)
          p_val <- 1 - CDF_U(contrast.matrix[as.numeric(first_sample), as.numeric(second_sample)])
        } else {
          Usim <- tabul_dist[[new.n_clust]]
          q_H <- stats::quantile(tabul_dist[[new.n_clust]], conf_level)
          CDF_U <- stats::ecdf(tabul_dist[[new.n_clust]])
          p_val <- 1 - CDF_U(contrast.matrix[as.numeric(first_sample), as.numeric(second_sample)])
        }
        CDF_U <- NULL
      }
      if (length(setdiff(indexesSamples_to_consider_new, indexesSamples_to_consider)) > 0) {
        ## check whether new population could be affected to the existing cluster by looking at results from pairwise equality tests:
        which_to_test <- indexesSamples_to_consider_new[indexesSamples_to_consider_new %in% as.numeric(clusters[[length(clusters)]]) == FALSE]
        couples_to_test <- apply(X = t(expand.grid(list(as.numeric(clusters[[length(clusters)]]), which_to_test))), 1, sort)
        if (!is.matrix(couples_to_test)) {
          ## case of only one pair of two populations:
          couples_to_test <- sort(couples_to_test)
          if (H0.test[couples_to_test[1],couples_to_test[2]] == TRUE) {
            k_sample_decision <- TRUE
          } else {
            k_sample_decision <- FALSE
            k_sample_pval <- p_val
          }
        } else {
          couples_to_test <- t(apply(X = couples_to_test, 1, sort))
          ## Case of list of two populations:
          if (all(H0.test[couples_to_test]) == TRUE) {
            k_sample_decision <- TRUE
          } else {
            ## K-sample test :
            comp_indices <- sort( c(2*indexesSamples_to_consider_new-1, 2*indexesSamples_to_consider_new) )
            k_sample_test <- IBM_k_samples_test(samples = samples[indexesSamples_to_consider_new],
                                                admixMod = admixMod[indexesSamples_to_consider_new],
                                                conf_level = conf_level, parallel = parallel, n_cpu = n_cpu,
                                                sim_U = Usim, n_sim_tab = n_sim_tab, tune_penalty = tune_penalty)
            tab_distrib <- append(tab_distrib, list(k_sample_test$tabulated_dist))
            k_sample_decision <- k_sample_test$reject_decision
            k_sample_pval <- k_sample_test$p_value
          }
        }
      }
      if (!k_sample_decision) {
        alreadyGrouped_samples <- unique(c(alreadyGrouped_samples, first_sample, second_sample))
        clusters[[length(clusters)]] <- unique(append(clusters[[length(clusters)]], c(first_sample, second_sample) ))
        which_row <- which_col <- as.numeric(clusters[[length(clusters)]])
        neighboors <- c(discrepancy.id[which_row, ], discrepancy.id[ ,which_col])
        new.n_clust <- n_clust
        p_value <- c(p_value, k_sample_pval)
      } else {
        alreadyGrouped_samples <- unique(c(alreadyGrouped_samples, first_sample, second_sample))
        clusters <- list(clusters,
                         as.character(indexesSamples_to_consider_new[which((indexesSamples_to_consider_new %in% as.numeric(clusters[[length(clusters)]])) == F)]))
        which_row <- which_col <- as.numeric(clusters[[length(clusters)]])
        neighboors <- c(discrepancy.id[which_row, ], discrepancy.id[ ,which_col])
        new.n_clust <- n_clust + 1
        indexesSamples_to_consider_new <- as.numeric(clusters[[length(clusters)]])
        p_value <- c(p_value, 2e-16)
      }
    }
    indexesSamples_to_consider <- indexesSamples_to_consider_new
    if (echo) utils::setTxtProgressBar(prog_bar, length(alreadyGrouped_samples))
  } # End of While

  if (echo) base::close(prog_bar)

  ## Define function that copies upper triangle to lower triangle to make the matrix become symmetric:
  f <- function(mat) {
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    mat
  }
  symmetric_dist_mat <- f(contrast.matrix)
  diag(symmetric_dist_mat) <- 0

  n_clust_final <- length(rapply(clusters, function(x) length(x)))
  clusters_affiliation <- matrix(NA, nrow = 2, ncol = length(samples))
  clusters_affiliation[1, ] <- as.numeric(rapply(clusters, function(x) utils::head(x,length(samples))))
  clusters_affiliation[2, ] <- rep(1:n_clust_final, times = rapply(clusters, function(x) length(x)))
  clusters_affiliation <- clusters_affiliation[ ,order(clusters_affiliation[1, ])]
  rownames(clusters_affiliation) <- c("Id_sample","Id_cluster")

  clusters_components <- clusters_weights <- clusters_couples <- vector(mode = "list", length = n_clust_final)
  for (k in 1:n_clust_final) { clusters_components[[k]] <- which(clusters_affiliation["Id_cluster", ] == k) }
  order_clust <- order(sapply(X = clusters_components, FUN = min))
  for (l in 1:n_clust_final) {
    res <- vector(mode = "numeric", length = nrow(couples.list))
    for (m in 1:nrow(couples.list)) {
      res[m] <- all(couples.list[m, ] %in% clusters_components[[order_clust[l]]]) == TRUE
    }
    clusters_couples[[order_clust[l]]] <- couples.list[which(res == TRUE), ]
    clusters_weights[[order_clust[l]]] <- weights.list[which(res == TRUE), ]
  }

  obj <- list(
    n_populations = length(samples),
    population_sizes = sapply(X = samples, FUN = length),
    admixture_models = admixMod,
    n_clust = n_clust_final,
    pval_clust = round(p_value, 3),
    clusters = clusters_affiliation,
    confidence_level = conf_level,
    clust_pop = clusters_components,
    clust_sizes = sapply(X = clusters_components, length),
    clust_weights = clusters_weights,
    discrepancy_matrix = symmetric_dist_mat,
    tab_distributions = unique(tab_distrib)
    )
  class(obj) <- "admix_cluster"
  obj$call <- match.call()

  on.exit(base::options(warn = old_options_warn))
  return(obj)
}


#' Print method for object of class 'admix_cluster'
#'
#' Print the main results when clustering the unknown component distributions coming from various
#' admixture samples, i.e. the obtained clusters.
#'
#' @param x An object of class 'admix_cluster' (see ?admix_clustering).
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_cluster <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of detected clusters across the samples provided: ", x$n_clust, ".\n", sep = "")
  cat("\nList of samples involved in each built cluster (in numeric format, i.e. inside c()) :\n",
      paste("  - Cluster #", 1:length(x$clust_pop), ": vector of populations ", x$clust_pop, collapse="\n", sep = ""), sep="")
  cat("\n")
}


#' Summary method for object of class 'admix_cluster'
#'
#' Summarizes the results obtained when clustering the unknown component distributions coming from various
#' admixture samples.
#'
#' @param object An object of class 'admix_cluster' (see ?admix_clustering).
#' @param ... further arguments passed to or from other methods.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.admix_cluster <- function(object, ...)
{
  cat("Call:\n")
  print(object$call)
  cat("\n--------- About samples ---------\n")
  cat("The number of populations/samples under study is ", object$n_populations, ".\n", sep = "")
  cat(paste("Size of sample ", 1:object$n_populations, ": ", object$population_sizes, sep = ""), sep = "\n")
  cat("\n-------- About contamination (admixture) models -------")
  cat("\n")
  for (k in 1:object$n_populations) {
    cat("-> Distribution and parameters of the known component \n for admixture model #", k, ": ", sep="")
    cat(paste(sapply(object$admixture_models[[k]], "[[", "known")[1:2], collapse = " - "))
    cat("\n")
  }
  cat("\n--------- About clustering ---------\n")
  cat("The level of the underlying k-sample testing procedure is set to ", (1-object$confidence_level)*100, "%.", sep = "")
  cat("\nNumber of detected clusters across the samples provided: ", object$n_clust, ".", sep = "")
  cat("\nThe p-values of the k-sample tests (showing when to close the clusters (i.e. p-value < ", (1-object$confidence_level), ") equal: ",
      paste(object$pval_clust, collapse=", "), ".", sep="")
  cat("\n\nList of samples involved in each built cluster (in numeric format, i.e. inside c()) :\n",
      paste("  - Cluster #", 1:length(object$clust_pop), ": vector of populations ", object$clust_pop, collapse="\n", sep = ""), sep="")
  weights.list <- vector(mode = "list", length = length(object$clust_weights))
  for (i in 1:length(object$clust_weights)) {
    if (object$clust_sizes[i] > 2) {
      weights.list[[i]] <- c(round(object$clust_weights[[i]][cumsum(c(1,(object$clust_sizes[i]-1):2)),1], 2),
                             round(object$clust_weights[[i]][nrow(object$clust_weights[[i]]),2], 2))
    } else {
      weights.list[[i]] <- c(round(object$clust_weights[[i]], 2))
    }
  }
  cat("\n\nList of estimated weights for the unknown component distributions in each detected cluster
      (in the same order as listed samples in each detected clusters) :\n",
      paste("  - estimated weights of the unknown component distributions for cluster ", 1:length(object$clust_pop), ": ", weights.list, collapse="\n"), sep="")
  cat("\n\nMatrix providing the discrepancies between populations, used in the clustering procedure:\n")
  print(object$discrepancy_matrix)
}
