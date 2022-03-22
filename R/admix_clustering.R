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
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when parallelizing.
#'
#' @details See the paper at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03201760
#'
#' @return A list with three elements: 1) the identified clusters; 2) the cluster affiliation; 3) the discrepancy matrix.
#'
#' @examples
#' \donttest{
#' ## Simulate data (chosen parameters indicate 2 clusters (populations (1,3), (2,4))!):
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
#'                              comp.dist=list.comp, comp.param=list.param, parallel=FALSE, n_cpu=2)
#' clusters$clustering
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

admix_clustering <- function(samples = NULL, n_sim_tab = 100, comp.dist = NULL, comp.param = NULL, parallel = FALSE, n_cpu = 2)
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
  #q_H <- numeric(length = nrow(couples.list))
  #for (k in 1:nrow(couples.list)) {
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
#    empirical.contr[k] <- minimal_size * IBM_empirical_contrast(XY[["prop.estim"]], fixed.p.X = XY[["p.X.fixed"]], sample1 = samples[[couples.list[k, ][1]]],
#                                                                sample2 = samples[[couples.list[k, ][2]]], G = XY[["integ.supp"]],
#                                                                comp.dist = couples.expr[[k]], comp.param = couples.param[[k]])
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
  contrast.matrix <- discrepancy.id <- discrepancy.rank <- matrix(NA, nrow = length(model.list), ncol = length(model.list))
  couples.expr <- couples.param <- vector(mode = "list", length = nrow(couples.list))
  for (k in 1:nrow(couples.list)) {
    couples.expr[[k]] <- comp.dist[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
    couples.param[[k]] <- comp.param[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
    names(couples.expr[[k]]) <- names(couples.param[[k]]) <- c("f1","g1","f2","g2")
#    contrast.matrix[couples.list[k, ][1], couples.list[k, ][2]] <- empirical.contr[k]
    contrast.matrix[couples.list[k, ][1], couples.list[k, ][2]] <- empirical.contr[[k]]
    discrepancy.id[couples.list[k, ][1], couples.list[k, ][2]] <- paste(couples.list[k, ][1], couples.list[k, ][2], sep = "-")
#    discrepancy.rank[couples.list[k, ][1], couples.list[k, ][2]] <- rank(empirical.contr)[k]
    discrepancy.rank[couples.list[k, ][1], couples.list[k, ][2]] <- rank(unlist(empirical.contr))[k]
  }

  ## Start algorithm to create the clusters :
  clusters <- vector(mode = "list")
  which_row <- unlist(apply(discrepancy.rank, 2, function(x) which(x == 1)))
  which_col <- unlist(apply(discrepancy.rank, 1, function(x) which(x == 1)))
  neighboors_index_init <- which.min(contrast.matrix)
  first.group <- discrepancy.id[which_row, which_col]
  U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 80, sample1 = samples[[which_row]], sample2 = samples[[which_col]],
                                 min_size = minimal_size, comp.dist = couples.expr[[which((couples.list[ ,1] == which_row) & (couples.list[ ,2] == which_col))]],
                                 comp.param = couples.param[[which((couples.list[ ,1] == which_row) & (couples.list[ ,2] == which_col))]],
                                 parallel = parallel, n_cpu = n_cpu)
  q_H <- stats::quantile(U[["U_sim"]], 0.95)
  test.H0 <- contrast.matrix[which_row, which_col] > q_H

  if (!test.H0) {
    clusters[[1]] <- alreadyGrouped_samples <- as.character(c(which_row, which_col))
    neighboors <- c(discrepancy.id[which_row, ], discrepancy.id[which_col, ], discrepancy.id[ ,which_row], discrepancy.id[ ,which_col])
    new.n_clust <- n_clust <- 1
    ## Look for a another member to integrate the newly built cluster:
    indexesSamples_to_consider <- c(0,0,0)
  } else {
    clusters[[1]] <- as.character(which_row)
    clusters[[2]] <- as.character(which_col)
    alreadyGrouped_samples <- c(as.character(which_row), as.character(which_col))
    neighboors <- c(discrepancy.id[which_col, ], discrepancy.id[ ,which_col])
    n_clust <- 1
    new.n_clust <- 2
    ## Look for a second member to integrate the newly second built cluster:
    indexesSamples_to_consider <- c(0,0)
  }

  while (length(alreadyGrouped_samples) < length(samples)) {
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
        U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 80, sample1 = samples[[as.numeric(first_sample)]],
                sample2 = samples[[as.numeric(second_sample)]], min_size = minimal_size,
                comp.dist = couples.expr[[which((couples.list[ ,1] == as.numeric(first_sample)) & (couples.list[ ,2] == as.numeric(second_sample)))]],
                comp.param = couples.param[[which((couples.list[ ,1] == as.numeric(first_sample)) & (couples.list[ ,2] == as.numeric(second_sample)))]],
                parallel = parallel, n_cpu = n_cpu)
        q_H <- stats::quantile(U[["U_sim"]], 0.95)
      }
      if (any(indexesSamples_to_consider_new != indexesSamples_to_consider)) {
        ## K-sample test :
        comp_indices <- sort( c(2*indexesSamples_to_consider_new-1, 2*indexesSamples_to_consider_new) )
        k_sample_test <- IBM_k_samples_test(samples = samples[indexesSamples_to_consider_new], sim_U = U[["U_sim"]],
                                            n_sim_tab = n_sim_tab, min_size = minimal_size, comp.dist = comp.dist[comp_indices],
                                            comp.param = comp.param[comp_indices], parallel = parallel, n_cpu = n_cpu)
      }
      if (!k_sample_test$rejection_rule) {
        alreadyGrouped_samples <- unique(c(alreadyGrouped_samples, first_sample, second_sample))
        clusters[[length(clusters)]] <- unique(append(clusters[[length(clusters)]], c(first_sample, second_sample) ))
        which_row <- which_col <- as.numeric(clusters[[length(clusters)]])
        neighboors <- c(discrepancy.id[which_row, ], discrepancy.id[ ,which_col])
        new.n_clust <- n_clust
      } else {
        alreadyGrouped_samples <- unique(c(alreadyGrouped_samples, first_sample, second_sample))
        clusters <- list(clusters,
                         as.character(indexesSamples_to_consider_new[which((indexesSamples_to_consider_new %in% as.numeric(clusters[[length(clusters)]])) == F)]))
        which_row <- which_col <- as.numeric(clusters[[length(clusters)]])
        neighboors <- c(discrepancy.id[which_row, ], discrepancy.id[ ,which_col])
        new.n_clust <- n_clust + 1
        indexesSamples_to_consider_new <- as.numeric(clusters[[length(clusters)]])
      }
    }
#    print(clusters)
    indexesSamples_to_consider <- indexesSamples_to_consider_new

  } # End of While

  ## Define function that copies upper triangle to lower triangle to make the lmatrix become symmetric:
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

  return( list(n_clust = n_clust_final, clusters = clusters_affiliation, discrepancy_matrix = symmetric_dist_mat) )
}
