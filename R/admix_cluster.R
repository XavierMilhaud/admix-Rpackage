#' Cluster K populations following admixture models
#'
#' Create clusters on the unknown components related to the K populations following admixture models. Based on the K-sample test
#' using Inversion - Best Matching (IBM) approach, see 'Details' below for further information.
#'
#' @param samples A list of the K (K>1) samples to be studied, all following admixture distributions.
#' @param admixMod A list of objects of class \link[admix]{admix_model}, containing useful information about distributions and parameters.
#' @param conf_level (default to 0.95) The confidence level of the k-sample tests used in the clustering procedure.
#' @param tune_penalty (default to TRUE) A boolean that allows to choose between a classical penalty term or an optimized penalty (embedding
#'                     some tuning parameters, automatically optimized). Optimized penalty is particularly useful for low/mid-sized samples,
#'                     or unbalanced sample sizes to detect alternatives to the null hypothesis (H0). It is recommended to use it.
#' @param tabul_dist (default to NULL) Only useful for comparisons of detected clusters at different confidence levels. A list of
#'                    the tabulated distributions of the stochastic integral used in the k-sample test, each element for each
#'                    cluster previously detected.
#' @param echo (default to TRUE) Display the remaining computation time.
#' @param ... Optional arguments to \link[admix]{IBM_k_samples_test}; namely 'n_sim_tab', 'parallel' and 'n_cpu'. These are crucial
#'            to speed-up the building of clusters.
#'
#' @references
#' \insertRef{MilhaudPommeretSalhiVandekerkhove2024b}{admix}
#'
#' @return An object of class \link[admix]{admix_cluster}, containing 12 attributes: 1) the number of samples under study; 2) the sizes of samples;
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
#'               conf_level = 0.95, tune_penalty = TRUE, n_sim_tab = 30)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_cluster <- function(samples, admixMod, conf_level = 0.95, tune_penalty = TRUE,
                          tabul_dist = NULL, echo = TRUE, ...)
{
  if (!all(sapply(X = admixMod, FUN = inherits, what = "admix_model")))
    stop("Argument 'admixMod' is not correctly specified. See ?admix_model.")

  old_options_warn <- base::options()$warn
  on.exit(base::options(warn = old_options_warn))
  base::options(warn = -1)

  if (length(sapply(samples, length)) == 1) return("One single sample, no clusters to be found.")
  ## Get the minimal size among all sample sizes, useful for future tabulation (adjustment of variance-covariance):
  minimal_size <- min(sapply(X = samples, FUN = length))

  ## K*(K-1)/2 combinations of populations under study:
  couples.list <- NULL
  for (i in 1:(length(samples)-1)) {
    for (j in (i+1):length(samples)) {
      couples.list <- rbind(couples.list,c(i,j))
    }
  }

  weights.list <- matrix(data = NA, nrow = nrow(couples.list), ncol = 2)
  empirical.contr <- vector(mode = "list", length = nrow(couples.list))
  for (k in 1:nrow(couples.list)) {
    ## Comparison between the two populations :
    XY <- estim_IBM(samples = samples[as.numeric(couples.list[k, ])],
                    admixMod = admixMod[as.numeric(couples.list[k, ])], n.integ = 1000)
    while (!exists("XY")) {
      XY <- estim_IBM(samples = samples[as.numeric(couples.list[k, ])],
                      admixMod = admixMod[as.numeric(couples.list[k, ])], n.integ = 1000)
    }
    weights.list[k, ] <- XY$estimated_mixing_weights
    empirical.contr[[k]] <- minimal_size * IBM_empirical_contrast(XY$estimated_mixing_weights, samples = samples[as.numeric(couples.list[k, ])],
                                                                  admixMod = admixMod[as.numeric(couples.list[k, ])], G = XY$integ.supp, fixed.p.X = XY$p.X.fixed)
  }
  ## Manage cases when the optimization algorithm could not find a solution:
  if (length(which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE)) != 0) {
    for (k in 1:length(which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE))) {
      empirical.contr[[which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE)[k]]] <- NA
    }
  }

  ## Give a simple and useful representation of results for each couple:
  contrast.matrix <- discrepancy.id <- discrepancy.rank <- matrix(NA, nrow = length(samples), ncol = length(samples))
  for (k in 1:nrow(couples.list)) {
    contrast.matrix[couples.list[k, ][1], couples.list[k, ][2]] <- empirical.contr[[k]]
    discrepancy.id[couples.list[k, ][1], couples.list[k, ][2]] <- paste(couples.list[k, ][1], couples.list[k, ][2], sep = "-")
    discrepancy.rank[couples.list[k, ][1], couples.list[k, ][2]] <- rank(unlist(empirical.contr))[k]
  }

  ## Start algorithm to create the clusters :
  clusters <- tab_distrib <- vector(mode = "list")
  closest_couple <- couples.list[which.min(empirical.contr), ]
  tab_distrib <- NULL
  if (is.null(tabul_dist)) {
    twoSample_test <- IBM_k_samples_test(samples = samples[closest_couple], admixMod = admixMod[closest_couple],
                                         conf_level = conf_level, sim_U = NULL, tune_penalty = FALSE, ...)
    tab_distrib <- twoSample_test$tabulated_dist
  } else {
    twoSample_test <- IBM_k_samples_test(samples = samples[closest_couple], admixMod = admixMod[closest_couple],
                                         conf_level = conf_level, sim_U = tabul_dist[[1]], tune_penalty = FALSE, ...)
    tab_distrib <- tabul_dist[[1]]
  }

  Usim <- twoSample_test$tabulated_dist
  if (!all(is.na(Usim))) {
    CDF_U <- stats::ecdf(Usim)
    p_value <- 1 - CDF_U(empirical.contr[[which.min(empirical.contr)]])
  } else {
    Usim <- NULL
    p_value <- 1e-16
  }

  if (!twoSample_test$reject_decision) {
    CDF_U <- NULL
    clusters[[1]] <- as.character(closest_couple)
    ## Detect which are the neighboors of the current cluster under study:
    index_neighboors_couples <- c(1:nrow(couples.list))[-c(which.min(empirical.contr),
                                  which(rowSums(t(apply(couples.list, 1, is.element, set=closest_couple))) == 0))]
    neighboors_couples <- couples.list[index_neighboors_couples, ]
    n_clust <- 1
  } else {
    clusters[[1]] <- as.character(closest_couple[1])
    clusters[[2]] <- as.character(closest_couple[2])
    ## Detect which are the neighboors of the current population under study:
    index_neighboors_couples <- c(1:nrow(couples.list))[-c(which.min(empirical.contr),
                                  which(rowSums(t(apply(couples.list, 1, is.element, set=closest_couple[2]))) == 0))]
    neighboors_couples <- couples.list[index_neighboors_couples, ]
    n_clust <- 2
  }

  if (echo) {
    prog_bar <- utils::txtProgressBar(min = 0, max = length(samples), style = 3, width = 50, char = "=")
    utils::setTxtProgressBar(prog_bar, 2/length(samples))
  }
  alreadyGrouped_samples <- as.character(closest_couple)

  while (length(alreadyGrouped_samples) < length(samples)) {

    nearest_neighboor_index <- which.min(empirical.contr[index_neighboors_couples])
    first_sample <- as.character(ifelse(is.null(dim(neighboors_couples)), neighboors_couples[1],
                                        neighboors_couples[nearest_neighboor_index,1]))
    second_sample <- as.character(ifelse(is.null(dim(neighboors_couples)), neighboors_couples[2],
                                         neighboors_couples[nearest_neighboor_index,2]))
    new_sample <- setdiff(x = c(first_sample,second_sample), y = alreadyGrouped_samples)
    alreadyGrouped_samples <- c(alreadyGrouped_samples, new_sample)

    ## k-sample test (can be 2-sample test):
    index_samples <- unique(sort(as.numeric(c(clusters[[length(clusters)]], first_sample, second_sample))))
    if (is.null(tabul_dist)) {
      k_sample_test <- IBM_k_samples_test(samples = samples[index_samples], admixMod = admixMod[index_samples],
                                          conf_level = conf_level, sim_U = Usim, tune_penalty = tune_penalty, ...)
      tab_distrib <- append(tab_distrib, list(k_sample_test$tabulated_dist))
    } else {
      k_sample_test <- IBM_k_samples_test(samples = samples[index_samples], admixMod = admixMod[index_samples],
                                          conf_level = conf_level, sim_U = tabul_dist[[n_clust]], tune_penalty = tune_penalty, ...)
      tab_distrib <- append(tab_distrib, tabul_dist[[n_clust]])
    }
    p_value <- c(p_value, k_sample_test$p_value)

    if (k_sample_test$reject_decision) {
      n_clust <- n_clust + 1
      clusters <- append(clusters, new_sample)
      ## update tabulated distribution to unknown one:
      Usim <- NULL
    } else {
      clusters[[length(clusters)]] <- c(clusters[[length(clusters)]], new_sample)
    }

    ## define new neighboors to study:
    index_neighboors_couples <- c(1:nrow(couples.list))[-unique(sort(c(
      c(1:nrow(couples.list))[-which(lapply(apply(couples.list, 1, setdiff, y = as.numeric(alreadyGrouped_samples)), length) > 0)], # already in existing clusters
      which(apply(couples.list, 1, setequal, y = as.numeric(c(first_sample, second_sample)))), # last couple studied
      which(rowSums(t(apply(couples.list, 1, is.element, set = as.numeric(clusters[[length(clusters)]]) ))) == 0) # no link with studied sample
    )))]
    neighboors_couples <- couples.list[index_neighboors_couples, ]

    if (echo) utils::setTxtProgressBar(prog_bar, length(alreadyGrouped_samples))

  } # End of While

  if (echo) base::close(prog_bar)

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
    discrepancy_matrix = contrast.matrix,
    tab_distributions = lapply(X = unique(tab_distrib[-which(lapply(tab_distrib, length) == 1)]), FUN = sort)
  )
  class(obj) <- "admix_cluster"
  obj$call <- match.call()
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
  cat("\nNumber of detected clusters: ", x$n_clust, ".\n", sep = "")
  cat("List of samples involved in each built cluster:\n",
      gsub("\\)", "", gsub("c\\(", "", paste("  - Cluster #", 1:length(x$clust_pop), ": samples ",
                                             x$clust_pop, collapse="\n", sep = ""))))
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
  cat("\n------ About samples ------\n")
  cat("The number of populations/samples under study is ", object$n_populations, ".\n", sep = "")
  cat(paste("Size of sample ", 1:object$n_populations, ": ", object$population_sizes, sep = ""), sep = "\n")
  cat("\n----- About contamination (admixture) models -----")
  cat("\n")
  for (k in 1:object$n_populations) {
    cat("-> Distribution and parameters of the known component \n for admixture model #", k, ": ", sep="")
    cat(paste(sapply(object$admixture_models[[k]], "[[", "known")[1:2], collapse = " - "))
    cat("\n")
  }
  cat("\n------ About clustering ------\n")
  cat("Test level of the underlying k-sample testing procedure: ", (1-object$confidence_level)*100, "%.", sep = "")
  cat("\nNumber of detected clusters across the samples provided: ", object$n_clust, ".", sep = "")
  cat("\np-values of the k-sample tests (showing when to close the clusters (i.e. p-value < ", (1-object$confidence_level), ") equal: ",
      paste(object$pval_clust, collapse=", "), ".", sep="")
  cat("\n\nList of samples involved in each built cluster:\n",
      gsub("\\)", "", gsub("c\\(", "", paste("  - Cluster #", 1:length(object$clust_pop), ": samples ",
                                             object$clust_pop, collapse="\n", sep = ""))))
  weights.list <- vector(mode = "list", length = length(object$clust_weights))
  for (i in 1:length(object$clust_weights)) {
    if (object$clust_sizes[i] > 2) {
      weights.list[[i]] <- c(round(object$clust_weights[[i]][cumsum(c(1,(object$clust_sizes[i]-1):2)),1], 2),
                             round(object$clust_weights[[i]][nrow(object$clust_weights[[i]]),2], 2))
    } else {
      weights.list[[i]] <- c(round(object$clust_weights[[i]], 2))
    }
  }
  cat("\n\nList of estimated weights for the unknown distributions in each detected cluster
      (in the same order as listed samples in each detected clusters) :\n",
      gsub("\\)", "", gsub("c\\(", "", paste("- Estimated weights of the unknown distributions for cluster ",
                                             1:length(object$clust_pop), ": ", weights.list, collapse="\n", sep=""))))
  cat("\n\nMatrix of discrepancies between samples (used for clustering):\n")
  print(object$discrepancy_matrix)
}

