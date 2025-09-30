
#' Extractor for known component(s) in admixture model(s)
#'
#' Get the known component of the admixture model considered for estimation,
#' test, or clustering.
#'
#' @param x An object of class \code{admix_estim}, \code{gaussianity_test}, \code{orthobasis_test},
#'          \code{IBM_test}, or \code{admix_cluster}.
#'
#' @return A list providing information on the known component (distribution, parameters).
#'
#' @details
#' This is a generic extractor, providing with the same information whatever the object class.
#'
#' @examples
#' ## Simulate a two-component Gaussian mixture:
#' mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Estimate the unknown quantities:
#' x <- admix_estim(samples = list(data1), admixMod = list(admixMod1), est_method = "BVdk")
#' ## Extract the information about the known component:
#' get_known_component(x)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_known_component <- function(x)
{
  UseMethod("get_known_component",x)
}

get_known_component.admix_estim <- function(x)
{
  if (!inherits(x, "admix_estim")) stop("This method must be used with objects of class 'admix_estim'.")
  admix_mod <- vector(mode = "list")
  for (i in 1:length(x$estim_objects)) admix_mod[[i]] <- x$estim_objects[[i]]$admixture_models
  admix_mod
}
get_known_component.gaussianity_test <- function(x)
{
  if (!inherits(x, "gaussianity_test")) stop("This method must be used with objects of class 'gaussianity_test'.")
  x$admixture_models
}
get_known_component.orthobasis_test <- function(x)
{
  if (!inherits(x, "orthobasis_test")) stop("This method must be used with objects of class 'orthobasis_test'.")
  x$admixture_models
}
get_known_component.IBM_test <- function(x)
{
  if (!inherits(x, "IBM_test")) stop("This method must be used with objects of class 'IBM_test'.")
  x$admixture_models
}
get_known_component.admix_cluster <- function(x)
{
  if (!inherits(x, "admix_cluster")) stop("This method must be used with objects of class 'admix_cluster'.")
  x$admixture_models
}

#' Extractor for simulated data from two-component mixture
#'
#' Get the mixture data generated from method \code{twoComp_mixt()}.
#'
#' @param x An object of class \code{twoComp_mixt}.
#'
#' @return A numeric vector of the simulated data.
#'
#' @examples
#' sim.X <- twoComp_mixt(n = 20, weight = 0.5,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean"=3, "sd"=0.5),
#'                                         list("mean"=0, "sd"=1)))
#' get_mixture_data(sim.X)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_mixture_data <- function(x)
{
  UseMethod("get_mixture_data",x)
}
get_mixture_data.twoComp_mixt <- function(x)
{
  if (!inherits(x, "twoComp_mixt")) stop("This method must be used with objects of class 'twoComp_mixt'.")
  x$mixt.data
}

#' Extractor for estimated mixing weights
#'
#' Extracts the estimated mixing weights from fitted objects of class
#' \code{admix_estim}, \code{gaussianity_test} and \code{orthobasis_test}.
#'
#' @param x An object of class \code{admix_estim}, \code{gaussianity_test} or
#'        \code{orthobasis_test}.
#'
#' @return A numeric vector of estimated mixing weight(s).
#'
#' @details
#' This is a generic extractor. The exact behavior depends on the class
#' of the input object:
#' \itemize{
#'   \item \code{admix_estim}: returns the estimated mixture proportions.
#'   \item \code{gaussianity_test}, \code{orthobasis_test}: returns weights derived from hypothesis testing results.
#' }
#'
#' @examples
#' ## Simulate a two-component Gaussian mixture:
#' mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Estimate the unknown quantities:
#' x <- admix_estim(samples = list(data1), admixMod = list(admixMod1), est_method = "BVdk")
#' ## Extract the information about the known component:
#' get_mixing_weights(x)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_mixing_weights <- function(x)
{
  UseMethod("get_mixing_weights", x)
}

get_mixing_weights.admix_estim <- function(x)
{
  if (!inherits(x, "admix_estim")) stop("This method must be used with objects of class 'admix_estim'.")
  weights <- sapply(X = x$estim_objects, "[[", 'estimated_mixing_weights')
  unlist(Filter(Negate(is.null), weights))
}
get_mixing_weights.gaussianity_test <- function(x)
{
  if (!inherits(x, "gaussianity_test")) stop("This method must be used with objects of class 'gaussianity_test'.")
  as.numeric(x$estimate[base::grepl(pattern = "Weight", x = names(x$estimate))])
}
get_mixing_weights.orthobasis_test <- function(x)
{
  if (!inherits(x, "orthobasis_test")) stop("This method must be used with objects of class 'orthobasis_test'.")
  as.numeric(x$weights[base::grepl(pattern = "Weight", x = names(x$weights))])
}


#' Extractor for the test decision
#'
#' Provide the decision of the statistical test: reject or
#' do not reject the null hypothesis.
#'
#' @param x An object of class \code{gaussianity_test}, \code{orthobasis_test} or \code{IBM_test}.
#'
#' @return A boolean giving the result of the test, TRUE if the null
#'         hypothesis is rejected, otherwise FALSE.
#'
#' @examples
#' mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 350, weight = 0.85,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = -1, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' data2 <- get_mixture_data(mixt2)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' x <- admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'                 conf_level = 0.95, test_method = "poly", ask_poly_param = FALSE, support = "Real")
#' reject_nullHyp(x)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
reject_nullHyp <- function(x)
{
  UseMethod("reject_nullHyp",x)
}
reject_nullHyp.gaussianity_test <- function(x)
{
  if (!inherits(x, "gaussianity_test")) stop("This method must be used with objects of class 'gaussianity_test'.")
  as.logical(x$reject_decision)
}
reject_nullHyp.orthobasis_test <- function(x)
{
  if (!inherits(x, "orthobasis_test")) stop("This method must be used with objects of class 'orthobasis_test'.")
  as.logical(x$reject_decision)
}
reject_nullHyp.IBM_test <- function(x)
{
  if (!inherits(x, "IBM_test")) stop("This method must be used with objects of class 'IBM_test'.")
  as.logical(x$reject_decision)
}


#' Extractor for the selected rank in the test statistic
#'
#' Provide the selected rank of the test statistic (connected to the expansion order
#' of the densities in the orthonormal polynomial basis if method 'poly' was chosen;
#' or to the number of terms, i.e. discrepancies between couples of samples, included
#' in the test statistic with method 'icv').
#'
#' @param x An object of class \code{gaussianity_test}, \code{orthobasis_test} or \code{IBM_test}.
#'
#' @return An integer corresponding to the selected rank in the test statistics,
#'         i.e. how many terms were kept in the test statistic.
#'
#' @examples
#' mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 350, weight = 0.85,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = -1, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' data2 <- get_mixture_data(mixt2)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' x <- admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'                 conf_level = 0.95, test_method = "poly", ask_poly_param = FALSE, support = "Real")
#' which_rank(x)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
which_rank <- function(x)
{
  UseMethod("which_rank",x)
}
which_rank.gaussianity_test <- function(x)
{
  if (!inherits(x, "gaussianity_test")) stop("This method must be used with objects of class 'gaussianity_test'.")
  as.numeric(x$selected_rank)
}
which_rank.orthobasis_test <- function(x)
{
  if (!inherits(x, "orthobasis_test")) stop("This method must be used with objects of class 'orthobasis_test'.")
  as.numeric(x$selected_rank)
}
which_rank.IBM_test <- function(x)
{
  if (!inherits(x, "IBM_test")) stop("This method must be used with objects of class 'IBM_test'.")
  if (!is.na(x$selected_rank)) { as.numeric(x$selected_rank)
  } else { print("Two-sample test, hence no rank to be selected!") }
}


#' Extractor for tabulated distribution in the k-sample test
#'
#' Provide (the list of) tabulated distribution(s) that allow to define
#' the extreme quantile(s) against which the test statistic(s) is compared.
#'
#' @param x An object of class \code{IBM_test} or \code{admix_cluster}.
#'
#' @return A numeric vector containing the simulated values of the tabulated
#'         distribution, sorted in increasing order.
#'
#' @examples
#' mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 350, weight = 0.85,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = -1, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' data2 <- get_mixture_data(mixt2)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' x <- admix_test(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'                 conf_level = 0.95, test_method = "icv", n_sim_tab = 10)
#' get_tabulated_dist(x)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_tabulated_dist <- function(x)
{
  UseMethod("get_tabulated_dist",x)
}
get_tabulated_dist.IBM_test <- function(x)
{
  if (!inherits(x, "IBM_test")) stop("This method must be used with objects of class 'IBM_test'.")
  sort(x$tabulated_dist)
}
get_tabulated_dist.admix_cluster <- function(x)
{
  if (!inherits(x, "admix_cluster")) stop("This method must be used with objects of class 'admix_cluster'.")
  x$tab_distributions
}


#' Extractor for pairwise discrepancy rankings
#'
#' Provide the matrix storing the ranks of discrepancies using Inversion-Best Matching
#' approach between all couples among the K (K>2) samples under study.
#'
#' @param x An object of class \code{IBM_test}.
#'
#' @return A matrix of ranks, from the closest couple (rank 1) in terms of discrepancy
#'         measure to the most different one.
#'
#' @examples
#' mixt1 <- twoComp_mixt(n = 380, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 350, weight = 0.85,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = -1, "sd" = 1)))
#' mixt3 <- twoComp_mixt(n = 500, weight = 0.6,
#'                       comp.dist = list("gamma", "gamma"),
#'                       comp.param = list(list("shape" = 16, "scale" = 1/4),
#'                                         list("shape" = 12, "scale" = 1/2)))
#' data1 <- get_mixture_data(mixt1)
#' data2 <- get_mixture_data(mixt2)
#' data3 <- get_mixture_data(mixt3)
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
#'                          knownComp_param = mixt3$comp.param[[2]])
#' x <- admix_test(samples = list(data1,data2,data3),
#'                 admixMod = list(admixMod1,admixMod2,admixMod3),
#'                 conf_level = 0.95, test_method = "icv", n_sim_tab = 10)
#' get_discrepancy_rank(x)
#' get_discrepancy_matrix(x)
#' get_statistic_components(x)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_discrepancy_rank <- function(x)
{
  UseMethod("get_discrepancy_rank",x)
}
get_discrepancy_rank.IBM_test <- function(x)
{
  if (!inherits(x, "IBM_test")) stop("This method must be used with objects of class 'IBM_test'")
  if (!all(is.na(x$discrepancy_rank))) { x$discrepancy_rank
  } else { print("Two-sample test, hence no discrepancy ranks to be stored.") }
}


#' Extractor for discrepancies b/w unknown components
#'
#' Provide the matrix storing pairwise discrepancies b/w unknown components
#' in admixture models, using Inversion-Best Matching approach.
#'
#' @param x An object of class \code{IBM_test} or \code{admix_cluster}.
#'
#' @return A matrix of pairwise discrepancies among the K (K>2) samples
#'         under study.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_discrepancy_matrix <- function(x)
{
  UseMethod("get_discrepancy_matrix",x)
}
get_discrepancy_matrix.IBM_test <- function(x)
{
  if (!inherits(x, "IBM_test")) stop("This method must be used with objects of class 'IBM_test'.")
  if (!all(is.na(x$discrepancy_matrix))) { x$discrepancy_matrix
  } else { print("Two-sample test, hence no matrix for discrepancies!") }
}
get_discrepancy_matrix.admix_cluster <- function(x)
{
  if (!inherits(x, "admix_cluster")) stop("This method must be used with objects of class 'admix_cluster'.")
  x$discrepancy_matrix
}


#' Extractor for components involved in test statistic
#'
#' Provide the terms (or discrepancies) that compose the test statistic
#' for the k-sample test.
#'
#' @param x An object of class \code{IBM_test}.
#'
#' @return The components finally included in the test statistic, i.e. the discrepancies
#'         of the couples that were aggregated in the built sequence of statistics.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_statistic_components <- function(x)
{
  UseMethod("get_statistic_components",x)
}
get_statistic_components.IBM_test <- function(x)
{
  if (!inherits(x, "IBM_test")) stop("This method must be used with objects of class 'IBM_test'.")
  if (!is.na(x$statistic_name)) { x$statistic_name
  } else { print("Two-sample test, hence only one component!") }
}


#' Extractor for members of clusters
#'
#' Extract the clusters that were discovered among K samples, where belonging to
#' one given cluster means having equal unknown component distributions.
#'
#' @param x An object of class \code{admix_cluster}.
#'
#' @return The samples included in each detected cluster.
#'
#' @examples
#' \donttest{
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1600, weight = 0.8,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 16, "scale" = 1/4),
#'                                         list("rate" = 1/3.5)))
#' mixt2 <- twoComp_mixt(n = 2000, weight = 0.7,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 14, "scale" = 1/2),
#'                                         list("rate" = 1/5)))
#' mixt3 <- twoComp_mixt(n = 2500, weight = 0.6,
#'                       comp.dist = list("gamma", "gamma"),
#'                       comp.param = list(list("shape" = 16, "scale" = 1/4),
#'                                         list("shape" = 12, "scale" = 1/2)))
#' mixt4 <- twoComp_mixt(n = 3800, weight = 0.5,
#'                       comp.dist = list("gamma", "exp"),
#'                       comp.param = list(list("shape" = 14, "scale" = 1/2),
#'                                         list("rate" = 1/7)))
#' data1 <- get_mixture_data(mixt1) ; data2 <- get_mixture_data(mixt2)
#' data3 <- get_mixture_data(mixt3) ; data4 <- get_mixture_data(mixt4)
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
#' x <- admix_cluster(samples = list(data1, data2, data3, data4),
#'               admixMod = list(admixMod1, admixMod2, admixMod3, admixMod4),
#'               conf_level = 0.95, tune_penalty = TRUE, n_sim_tab = 10)
#' get_cluster_members(x)
#' get_cluster_sizes(x)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_cluster_members <- function(x)
{
  UseMethod("get_cluster_members",x)
}
get_cluster_members.admix_cluster <- function(x)
{
  if (!inherits(x, "admix_cluster")) stop("This method must be used with objects of class 'admix_cluster'.")
  x$clusters
}


#' Extractor for cluster sizes
#'
#' Provide the number of samples in each cluster.
#'
#' @param x An object of class \code{admix_cluster}.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
get_cluster_sizes <- function(x)
{
  UseMethod("get_cluster_sizes",x)
}
get_cluster_sizes.admix_cluster <- function(x)
{
  if (!inherits(x, "admix_cluster")) stop("This method must be used with objects of class 'admix_cluster'")
  x$clust_sizes
}
