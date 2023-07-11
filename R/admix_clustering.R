#' Clustering of K populations following admixture models
#'
#' Create clusters on the unknown components related to the K populations following admixture models. Based on the K-sample test
#' using Inversion - Best Matching (IBM) approach, see 'Details' below for further information.
#'
#' @param samples A list of the K observed samples to be clustered, all following admixture distributions.
#' @param n_sim_tab Number of simulated gaussian processes used in the tabulation of the inner convergence distribution in the IBM approach.
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
#' @param tabul.dist Only useful for comparisons of detected clusters at different confidence levels. Is a list of the tabulated distributions
#'                   of the stochastic integral for each cluster previously detected.
#' @param conf.level The confidence level of the K-sample test used in the clustering procedure.
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when parallelizing.
#'
#' @details See the paper at the following HAL weblink: https://hal.science/hal-04129130
#'
#' @return A list with eight elements: 1) the number of populations under consideration; 2) the number of detected clusters;
#'         3) the list of p-values for each test performed; 4) the cluster affiliation for each population; 5) the chosen confidence
#'         level of statistical tests; 6) the cluster components; 7) the estimated weights of the unknown component distributions inside
#'         each cluster (remind that estimated weights are consistent only under the null); 8) the matrix of pairwise discrepancies
#'         among all populations.
#'
#' @examples
#' \donttest{
#' ## Simulate data (chosen parameters indicate 2 clusters (populations (1,3), and (2,4)):
#' list.comp <- list(f1 = "gamma", g1 = "exp",
#'                   f2 = "gamma", g2 = "exp",
#'                   f3 = "gamma", g3 = "gamma",
#'                   f4 = "gamma", g4 = "exp")
#' list.param <- list(f1 = list(shape = 16, rate = 4), g1 = list(rate = 1/3.5),
#'                    f2 = list(shape = 14, rate = 2), g2 = list(rate = 1/5),
#'                    f3 = list(shape = 16, rate = 4), g3 = list(shape = 12, rate = 2),
#'                    f4 = list(shape = 14, rate = 2), g4 = list(rate = 1/7))
#' A.sim <- rsimmix(n=2600, unknownComp_weight=0.8, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' B.sim <- rsimmix(n=3000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' C.sim <- rsimmix(n=3500, unknownComp_weight=0.6, comp.dist = list(list.comp$f3,list.comp$g3),
#'                  comp.param = list(list.param$f3, list.param$g3))$mixt.data
#' D.sim <- rsimmix(n=4800, unknownComp_weight=0.5, comp.dist = list(list.comp$f4,list.comp$g4),
#'                  comp.param = list(list.param$f4, list.param$g4))$mixt.data
#' ## Look for the clusters:
#' list.comp <- list(f1 = NULL, g1 = "exp",
#'                   f2 = NULL, g2 = "exp",
#'                   f3 = NULL, g3 = "gamma",
#'                   f4 = NULL, g4 = "exp")
#' list.param <- list(f1 = NULL, g1 = list(rate = 1/3.5),
#'                    f2 = NULL, g2 = list(rate = 1/5),
#'                    f3 = NULL, g3 = list(shape = 12, rate = 2),
#'                    f4 = NULL, g4 = list(rate = 1/7))
#' clusters <- admix_clustering(samples = list(A.sim,B.sim,C.sim,D.sim), n_sim_tab = 8,
#'                              comp.dist=list.comp, comp.param=list.param, conf.level = 0.95,
#'                              parallel = FALSE, n_cpu = 2)
#' clusters
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_clustering <- function(samples = NULL, n_sim_tab = 100, comp.dist = NULL, comp.param = NULL,
                             tabul.dist = NULL, conf.level = 0.95, parallel = FALSE, n_cpu = 2)
{
  ## Control whether parallel computations were asked for or not:
  if (parallel) {
    `%fun%` <- foreach::`%dopar%`
    doParallel::registerDoParallel(cores = n_cpu)
  } else {
    `%fun%` <- foreach::`%do%`
  }

  stopifnot( length(comp.dist) == (2*length(samples)) )
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    if ( (!all(sapply(comp.param, is.null)[seq.int(from = 2, to = length(comp.dist), by = 2)] == FALSE)) |
         (!all(sapply(comp.param, is.null)[seq.int(from = 1, to = length(comp.dist), by = 2)] == TRUE)) ) {
      stop("Component distributions and/or parameters must have been badly specified in the admixture models.")
    }
  }

  ## Get the minimal size among all sample sizes, useful for future tabulation (adjustment of variance-covariance):
  minimal_size <- min(sapply(X = samples, FUN = length))
  index_minimal_size <- which.min(sapply(X = samples, FUN = length))

  ## Look for all possible couples on which the discrepancy will be computed :
  model.list <- lapply(X = seq.int(from = 1, to = length(comp.dist), by = 2), FUN = seq.int, length.out = 2)
  ## K*(K-1)/2 combinations of populations under study:
  couples.list <- NULL
  for (i in 1:(length(samples)-1)) { for (j in (i+1):length(samples)) { couples.list <- rbind(couples.list,c(i,j)) } }

  couples.expr <- couples.param <- vector(mode = "list", length = nrow(couples.list))
  empirical.contr <-
  foreach::foreach (k = 1:nrow(couples.list), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
    couples.expr[[k]] <- comp.dist[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
    couples.param[[k]] <- comp.param[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
    names(couples.expr[[k]]) <- names(couples.param[[k]]) <- c("f1","g1","f2","g2")
    ## Comparison between the two populations :
    XY <- IBM_estimProp(sample1 = samples[[couples.list[k, ][1]]], sample2 = samples[[couples.list[k, ][2]]], known.prop = NULL,
                        comp.dist = couples.expr[[k]], comp.param = couples.param[[k]], with.correction = FALSE, n.integ = 1000)
    while (!exists("XY")) {
      XY <- IBM_estimProp(sample1 = samples[[couples.list[k, ][1]]], sample2 = samples[[couples.list[k, ][2]]], known.prop = NULL,
                          comp.dist = couples.expr[[k]], comp.param = couples.param[[k]], with.correction = FALSE, n.integ = 1000)
    }

    minimal_size * IBM_empirical_contrast(XY[["prop.estim"]], fixed.p.X = XY[["p.X.fixed"]], sample1 = samples[[couples.list[k, ][1]]],
                                          sample2 = samples[[couples.list[k, ][2]]], G = XY[["integ.supp"]],
                                          comp.dist = couples.expr[[k]], comp.param = couples.param[[k]])
  }
  ## Manage cases when the optimization algorithm could not find a solution:
  if (length(which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE)) != 0) {
    for (k in 1:length(which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE))) {
      empirical.contr[[which(lapply(X = empirical.contr, FUN = is.numeric) == FALSE)[k]]] <- NA
    }
  }

  ## Give a simple and useful representation of results for each couple:
  H0.test <- contrast.matrix <- discrepancy.id <- discrepancy.rank <- matrix(NA, nrow = length(model.list), ncol = length(model.list))
  ## Store results of pairwise tests for equality between unknown components (fasten the future research of potential clusters):
  couples.expr <- couples.param <- vector(mode = "list", length = nrow(couples.list))
  weights.list <- matrix(data = NA, nrow = nrow(couples.list), ncol = 2)
#  tabulated_dist <-  <- vector(mode = "list", length = nrow(couples.list))
  for (k in 1:nrow(couples.list)) {
    couples.expr[[k]] <- comp.dist[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
    couples.param[[k]] <- comp.param[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
    names(couples.expr[[k]]) <- names(couples.param[[k]]) <- c("f1","g1","f2","g2")
    contrast.matrix[couples.list[k, ][1], couples.list[k, ][2]] <- empirical.contr[[k]]
    discrepancy.id[couples.list[k, ][1], couples.list[k, ][2]] <- paste(couples.list[k, ][1], couples.list[k, ][2], sep = "-")
    discrepancy.rank[couples.list[k, ][1], couples.list[k, ][2]] <- rank(unlist(empirical.contr))[k]
    ## Pairwise testing: test H0 between the two considered populations.
    pairwise_H0test <- NULL
    pairwise_H0test <- IBM_2samples_test(samples = list(samples[[couples.list[k, ][1]]], samples[[couples.list[k, ][2]]]),
                                         known.p = NULL, comp.dist = couples.expr[[k]], comp.param = couples.param[[k]],
                                         sim_U = NULL, n_sim_tab = 30, min_size = NULL, conf.level = conf.level,
                                         parallel = parallel, n_cpu = n_cpu)
    weights.list[k, ] <- pairwise_H0test[["weights"]]
    H0.test[couples.list[k, ][1], couples.list[k, ][2]] <- pairwise_H0test$rejection_rule
#    tabulated_dist[[k]] <- pairwise_H0test$tabulated_dist
  }

  ## Start algorithm to create the clusters :
  clusters <- tab_distrib <- vector(mode = "list")
  which_row <- unlist(apply(discrepancy.rank, 2, function(x) which(x == 1)))
  which_col <- unlist(apply(discrepancy.rank, 1, function(x) which(x == 1)))
  neighboors_index_init <- which.min(contrast.matrix)
  first.group <- discrepancy.id[which_row, which_col]
  if (is.null(tabul.dist)) {
    U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 80, sample1 = samples[[which_row]], sample2 = samples[[which_col]],
                                   min_size = minimal_size, comp.dist = couples.expr[[which((couples.list[ ,1] == which_row) & (couples.list[ ,2] == which_col))]],
                                   comp.param = couples.param[[which((couples.list[ ,1] == which_row) & (couples.list[ ,2] == which_col))]],
                                   parallel = parallel, n_cpu = n_cpu)
    tab_distrib[[1]] <- U[["U_sim"]]
    Usim <- U[["U_sim"]]
    CDF_U <- stats::ecdf(U[["U_sim"]])
    q_H <- stats::quantile(U[["U_sim"]], conf.level)
  } else {
    Usim <- tabul.dist[[1]]
    tab_distrib[[1]] <- tabul.dist[[1]]
    CDF_U <- stats::ecdf(tabul.dist[[1]])
    q_H <- stats::quantile(tabul.dist[[1]], conf.level)
  }
  test.H0 <- contrast.matrix[which_row, which_col] > q_H
  ## Numerical vector storing the p-values associated to each test: each cluster is closed once the p-value lies below the H0-rejection
  ## threshold (1-conf.level). This enables to see whether we were close to terminate the cluster or not each time we add a new population.
  p_value <- numeric(length = 0L)
  if (!test.H0) {
    clusters[[1]] <- alreadyGrouped_samples <- as.character(c(which_row, which_col))
    neighboors <- c(discrepancy.id[which_row, ], discrepancy.id[which_col, ], discrepancy.id[ ,which_row], discrepancy.id[ ,which_col])
    new.n_clust <- n_clust <- 1
    ## Look for a another member to integrate the newly built cluster:
    indexesSamples_to_consider <- c(0,0,0)
    #CDF_U <- stats::ecdf(U[["U_sim"]])
    p_value <- c(p_value, 1 - CDF_U(contrast.matrix[which_row, which_col]))
    CDF_U <- NULL
    prog_bar <- utils::txtProgressBar(min = 0, max = length(samples), style = 3, width = 50, char = "=")
    utils::setTxtProgressBar(prog_bar, 2/length(samples))
  } else {
    clusters[[1]] <- as.character(which_row)
    clusters[[2]] <- as.character(which_col)
    alreadyGrouped_samples <- c(as.character(which_row), as.character(which_col))
    neighboors <- c(discrepancy.id[which_col, ], discrepancy.id[ ,which_col])
    n_clust <- 1
    new.n_clust <- 2
    ## Look for a second member to integrate the newly second built cluster:
    indexesSamples_to_consider <- c(0,0)
    prog_bar <- utils::txtProgressBar(min = 0, max = length(samples), style = 3,  width = 50, char = "=")
    utils::setTxtProgressBar(prog_bar, 1/length(samples))
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
        if ( is.null(tabul.dist) | (new.n_clust > length(tabul.dist)) ) {
          U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 80, sample1 = samples[[as.numeric(first_sample)]],
                                         sample2 = samples[[as.numeric(second_sample)]], min_size = minimal_size,
                                         comp.dist = couples.expr[[which((couples.list[ ,1] == as.numeric(first_sample)) & (couples.list[ ,2] == as.numeric(second_sample)))]],
                                         comp.param = couples.param[[which((couples.list[ ,1] == as.numeric(first_sample)) & (couples.list[ ,2] == as.numeric(second_sample)))]],
                                         parallel = parallel, n_cpu = n_cpu)
          ## Store the simulated tabulated distribution:
          Usim <- U[["U_sim"]]
          tab_distrib <- append(tab_distrib, list(Usim))
          q_H <- stats::quantile(U[["U_sim"]], conf.level)
          CDF_U <- stats::ecdf(U[["U_sim"]])
          p_val <- 1 - CDF_U(contrast.matrix[as.numeric(first_sample), as.numeric(second_sample)])
          CDF_U <- NULL
        } else {
          Usim <- tabul.dist[[new.n_clust]]
          q_H <- stats::quantile(tabul.dist[[new.n_clust]], conf.level)
          CDF_U <- stats::ecdf(tabul.dist[[new.n_clust]])
          p_val <- 1 - CDF_U(contrast.matrix[as.numeric(first_sample), as.numeric(second_sample)])
          CDF_U <- NULL
        }
      }
#      if (any(indexesSamples_to_consider_new != indexesSamples_to_consider)) {
      if (length(setdiff(indexesSamples_to_consider_new, indexesSamples_to_consider)) > 0) {
        ## First check whether the new population could eventually be affected to the existing cluster by looking at results from pairwise equality tests:
        which_to_test <- indexesSamples_to_consider_new[indexesSamples_to_consider_new %in% as.numeric(clusters[[length(clusters)]]) == FALSE]
        #couples_to_test <- apply(apply(X = expand.grid(list(clusters[[length(clusters)]], which_to_test)), 1, as.numeric), 1, sort)
        couples_to_test <- apply(X = t(expand.grid(list(as.numeric(clusters[[length(clusters)]]), which_to_test))), 1, sort)
        if (!is.matrix(couples_to_test)) {
          ## case of only one pair of two populations:
          couples_to_test <- sort(couples_to_test)
          if (H0.test[couples_to_test[1],couples_to_test[2]] == TRUE) { k_sample_decision <- TRUE
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
            #k_sample_test <- IBM_k_samples_test(samples = samples[indexesSamples_to_consider_new], sim_U = U[["U_sim"]],
            #                                    n_sim_tab = n_sim_tab, min_size = minimal_size, comp.dist = comp.dist[comp_indices],
            #                                    comp.param = comp.param[comp_indices], conf.level = conf.level,
            #                                    tune.penalty = TRUE, parallel = parallel, n_cpu = n_cpu)
            k_sample_test <- IBM_k_samples_test(samples = samples[indexesSamples_to_consider_new], sim_U = Usim, n_sim_tab = n_sim_tab,
                                                comp.dist = comp.dist[comp_indices], comp.param = comp.param[comp_indices],
                                                conf.level = conf.level, tune.penalty = TRUE, parallel = parallel, n_cpu = n_cpu)
            tab_distrib <- append(tab_distrib, list(k_sample_test[["sim_U"]]))
            k_sample_decision <- k_sample_test$rejection_rule
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
    utils::setTxtProgressBar(prog_bar, length(alreadyGrouped_samples))
  } # End of While
  base::close(prog_bar)

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

  obj <- list(n_popu = length(samples),
              n_clust = n_clust_final,
              pval_clust = round(p_value, 3),
              clusters = clusters_affiliation,
              confidence_level = conf.level,
              clust_pop = clusters_components,
              clust_sizes = sapply(X = clusters_components, length),
              clust_weights = clusters_weights,
              discrepancy_matrix = symmetric_dist_mat,
              tab_distributions = unique(tab_distrib))
  class(obj) <- "admix_cluster"
  obj$call <- match.call()

  return(obj)
  #return( list(confidence_level = conf.level, n_clust = n_clust_final, clusters = clusters_affiliation, discrepancy_matrix = symmetric_dist_mat,
  #             clust_pop = clusters_components, clust_couples = clusters_couples, clust_weights = clusters_weights) )
}
