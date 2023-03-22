#' Results of the clustering algorithm performed over the K populations following admixture models.
#'
#' Print the detected clusters among the populations under study. This method also prints the number of clusters,
#' the estimated weights of the unknown component distributions inside each cluster, and the discrepancy matrix.
#' The latter represents some kind of distance between the populations.
#'
#' @param x An object of class 'admix_cluster' (see ?admix_clustering).
#' @param ... further arguments passed to or from other methods.
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
#'                              comp.dist=list.comp, comp.param=list.param, conf.level = 0.95,
#'                              parallel=FALSE, n_cpu=2)
#' print(clusters)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.admix_cluster <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nThe number of populations/samples under study is ", x$n_popu, ".", sep = "")
  cat("\nThe level of the underlying k-sample testing procedure is set to ", (1-x$confidence_level)*100, "%.", sep = "")
  cat("\n\nThe number of detected clusters in these populations equals ", x$n_clust, ".", sep = "")
  cat("\n\nThe list of clusters with populations belonging to them (in numeric format, i.e. inside c()) :\n",
      paste("  - Cluster id. ", 1:length(x$clust_pop), ": ", x$clust_pop, collapse="\n", sep = ""))
  cat("\n\nThe list of estimated weights for the unknown component distributions in each detected cluster
      (in the same format and order as listed populations for clusters just above) :\n",
      paste("  - estimated weights of the unknown component distributions for cluster ", 1:length(x$clust_pop), ": ", x$clust_weights, collapse="\n"))
  cat("\n\nThe matrix giving the distances between populations, used in the clustering procedure through the k-sample tests:\n")
  print(x$discrepancy_matrix)
}
