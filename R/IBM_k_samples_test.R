#' Equality test of unknown component distributions in K admixture models, with IBM approach
#'
#' Test hypothesis on the unknown component of K (K > 1) admixture models using Inversion - Best Matching method.
#' K-samples test of the unknown component distribution in admixture models using Inversion - Best Matching
#' (IBM) method. Recall that we have K populations following admixture models, each one with probability
#' density functions (pdf) l_k = p_k*f_k + (1-p_k)*g_k, where g_k is the known pdf and l_k corresponds to the
#' observed sample. Perform the following hypothesis test:
#'    H0 : f_1 = ... = f_K  against  H1 : f_i differs from f_j (i diff j, and i,j in 1,...,K).
#'
#' @param samples A list of the samples to be studied, all following admixture distributions.
#' @param sim_U Random draws of the inner convergence part of the contrast as defined in the IBM approach (see 'Details' below).
#' @param n_sim_tab Number of simulated gaussian processes when tabulating the inner convergence distribution in the IBM approach.
#' @param min_size useful to provide the minimal size among all samples (needed to take into account the correction factor for the
#'                 variance-covariance assessment). Otherwise, useless.
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
#' @param parallel (default to FALSE) Boolean indicating whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used when parallelizing.
#'
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03201760
#'
#' @return A list of ten elements, containing: 1) the rejection decision; 2) the p-value of the test; 3) the terms involved in
#'         the test statistic; 4) the test statistic value; 5) the selected rank (number of terms involved in the test statistic);
#'         6) the value of the penalized test statistic; 7) the sorted contrast values; 8) the 95th-quantile of the contrast
#'         distribution; 9) the final terms of the statistic; and 10) the contrast matrix.
#'
#' @examples
#' \donttest{
#' ####### Under the alternative hypothesis H1 (with K=3 populations):
#' list.comp <- list(f1 = "norm", g1 = "norm",
#'                   f2 = "norm", g2 = "norm",
#'                   f3 = "norm", g3 = "norm")
#' list.param <- list(f1 = list(mean = 0, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 2, sd = 1), g2 = list(mean = 4, sd = 1.1),
#'                    f3 = list(mean = 0, sd = 1), g3 = list(mean = 3, sd = 0.8))
#' sim1 <- rsimmix(n = 5000, unknownComp_weight = 0.8, comp.dist = list(list.comp$f1,list.comp$g1),
#'                 comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' sim2 <- rsimmix(n= 5300, unknownComp_weight = 0.6, comp.dist = list(list.comp$f2,list.comp$g2),
#'                 comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' sim3 <- rsimmix(n = 5100, unknownComp_weight = 0.7, comp.dist = list(list.comp$f3,list.comp$g3),
#'                 comp.param = list(list.param$f3, list.param$g3))$mixt.data
#' ## Perform the 3-samples test in a real-life setting:
#' list.comp <- list(f1 = NULL, g1 = "norm",
#'                   f2 = NULL, g2 = "norm",
#'                   f3 = NULL, g3 = "norm")
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 4, sd = 1.1),
#'                    f3 = NULL, g3 = list(mean = 3, sd = 0.8))
#' obj <- IBM_k_samples_test(samples=list(sim1,sim2,sim3), sim_U=NULL, n_sim_tab=4, min_size=NULL,
#'                           comp.dist=list.comp, comp.param=list.param, parallel=FALSE, n_cpu=2)
#' obj$rejection_rule
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_k_samples_test <- function(samples = NULL, sim_U = NULL, n_sim_tab = 100, min_size = NULL, comp.dist = NULL,
                               comp.param = NULL, parallel = FALSE, n_cpu = 2)
{
  stopifnot( length(comp.dist) == (2*length(samples)) )
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    if ( (!all(sapply(comp.param, is.null)[seq.int(from = 2, to = length(comp.dist), by = 2)] == FALSE)) |
         (!all(sapply(comp.param, is.null)[seq.int(from = 1, to = length(comp.dist), by = 2)] == TRUE)) ) {
      stop("Component distributions and/or parameters must have been badly specified in the admixture models.")
    }
  }

  ## Get the minimal size among all sample sizes, useful for future tabulation (adjustment of variance-covariance):
  if (is.null(min_size)) {
    minimal_size <- min(sapply(X = samples, FUN = length))
  } else {
    minimal_size <- min_size
  }

  ## Look for all possible couples on which the discrepancy will be computed :
  model.list <- lapply(X = seq.int(from = 1, to = length(comp.dist), by = 2), FUN = seq.int, length.out = 2)
  couples.list <- NULL
  for (i in 1:(length(samples)-1)) {
    for (j in (i+1):length(samples)) { couples.list <- rbind(couples.list,c(i,j)) }
  }

  couples.expr <- couples.param <- vector(mode = "list", length = nrow(couples.list))
  empirical.contr <- q_H <- numeric(length = nrow(couples.list))

  for (k in 1:nrow(couples.list)) {
    couples.expr[[k]] <- comp.dist[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
    couples.param[[k]] <- comp.param[c(model.list[[couples.list[k, ][1]]], model.list[[couples.list[k, ][2]]])]
    names(couples.expr[[k]]) <- names(couples.param[[k]]) <- c("f1","g1","f2","g2")
    ## Comparison between the two populations :
    XY <- IBM_estimProp(sample1 = samples[[couples.list[k, ][1]]], sample2 = samples[[couples.list[k, ][2]]], known.prop = NULL,
                        comp.dist = couples.expr[[k]], comp.param = couples.param[[k]], with.correction = F, n.integ = 1000)
    sol.XY <- XY[["prop.estim"]]
    empirical.contr[k] <- minimal_size * IBM_empirical_contrast(sol.XY, fixed.p.X = XY[["p.X.fixed"]], sample1 = samples[[couples.list[k, ][1]]],
                                                            sample2 = samples[[couples.list[k, ][2]]], G = XY[["integ.supp"]],
                                                            comp.dist = couples.expr[[k]], comp.param = couples.param[[k]])
  }

  ## Give a simple and useful representation of results for each couple:
  contrast.matrix <- discrepancy.id <- discrepancy.rank <- matrix(NA, nrow = length(samples), ncol = length(samples))
  for (k in 1:nrow(couples.list)) {
    contrast.matrix[couples.list[k, ][1], couples.list[k, ][2]] <- empirical.contr[k]
    discrepancy.id[couples.list[k, ][1], couples.list[k, ][2]] <- paste(couples.list[k, ][1], couples.list[k, ][2], sep = "-")
    discrepancy.rank[couples.list[k, ][1], couples.list[k, ][2]] <- rank(empirical.contr)[k]
  }

  if (is.null(sim_U)) {
    which_row <- unlist(apply(discrepancy.rank, 2, function(x) which(x == 1)))
    which_col <- unlist(apply(discrepancy.rank, 1, function(x) which(x == 1)))
    U <- IBM_tabul_stochasticInteg(n.sim = n_sim_tab, n.varCovMat = 100, sample1 = samples[[which_row]], sample2 = samples[[which_col]], min_size = minimal_size,
                    comp.dist = couples.expr[[which((couples.list[ ,1] == which_row) & (couples.list[ ,2] == which_col))]],
                    comp.param = couples.param[[which((couples.list[ ,1] == which_row) & (couples.list[ ,2] == which_col))]],
                    parallel = parallel, n_cpu = n_cpu)
    sim_U <- U[["U_sim"]]
  }
  if (all(is.na(sim_U))) {
    return( list(rejection_rule = TRUE, p_val = NA, stat_name = NA, stat_value = NA, stat_rank = NA,
                 penalized_stat = NA, ordered_contrasts = NA, quantiles = NA, stat_terms = NA,
                 contr_mat = contrast.matrix) )
  } else {
    q_H <- rep(stats::quantile(sim_U, 0.95), times = length(q_H))
  }

  ## Sort the discrepancies and corresponding information:
  sorted_contrasts <- sort(empirical.contr)
  penalized.stats <- numeric(length = (length(sorted_contrasts)+1))
  stat.terms <- stat.terms_final <- character(length = nrow(couples.list))
  for (k in 1:length(empirical.contr)) {
    stat.terms[k] <- paste("d_", discrepancy.id[which(discrepancy.rank == k)], sep = "")
    stat.terms_final[k] <- paste(stat.terms[1:k], collapse = " + ", sep = "")
    ## k*quantile*log(n) ?
    penalized.stats[k+1] <- penalized.stats[k] + sorted_contrasts[k] - minimal_size^0.25
  }
  penalized.stats <- round(penalized.stats[-1], 1)

  selected.index <- which.max(penalized.stats)
  finalStat_name <- stat.terms_final[selected.index]
  finalStat_value <- cumsum(sorted_contrasts)[selected.index]
  final_test <- finalStat_value > q_H[1]
  CDF_U <- stats::ecdf(sim_U)
  p_value <- 1 - CDF_U(finalStat_value)

  return( list(rejection_rule = final_test, p_val = p_value, stat_name = finalStat_name, stat_value = finalStat_value,
               stat_rank = selected.index, penalized_stat = penalized.stats, ordered_contrasts = sorted_contrasts,
               quantiles = unique(q_H), stat_terms = stat.terms_final, contr_mat = contrast.matrix) )
}
