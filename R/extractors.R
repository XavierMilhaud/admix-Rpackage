
#' Extractor for objects of class 'admix_estim', 'htest' or 'admix_cluster'
#'
#' Get the known component(s) of admixture model(s) considered for estimation, test, or clustering.
#'
#' @param x An object of class 'admix_estim', 'htest', or 'admix_cluster'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_known_component <- function(x)
{
  if (!inherits(x, "admix_estim") & !inherits(x, "htest") & !inherits(x, "admix_cluster"))
    stop("This function must be used with objects of class 'admix_estim', 'htest' or 'admix_cluster'.")

  admix_mod <- vector(mode = "list")
  if (inherits(x, "admix_estim")) {
    for (i in 1:length(x$estim_objects))
      admix_mod[[i]] <- x$estim_objects[[i]]$admixture_models
  } else {
    admix_mod <- x$admixture_models
  }
  admix_mod
}


#' Extractor for objects of class 'admix_estim' or 'htest'
#'
#' Get the estimated unknown mixing proportion(s) related to the weight(s) of
#' the unknown component distribution(s) of the admixture model(s).
#'
#' @param x An object of class 'admix_estim' or 'htest'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_mixing_weights <- function(x)
{
  if (!inherits(x, "admix_estim") & !inherits(x, "htest"))
    stop("This function must be used with objects of class 'admix_estim' or 'htest'.")
  if (inherits(x, "admix_estim")) {
    weights <- sapply(X = x$estim_objects, "[[", 'estimated_mixing_weights')
    unlist(Filter(Negate(is.null), weights))
  } else {  # case of htest objects
    as.numeric(x$estimate[base::grepl(pattern = "Weight", x = names(x$estimate))])
  }
}


#' Extractor for object of class 'htest'
#'
#' Provide the decision of the statistical test: reject or
#' do not reject the null hypothesis.
#'
#' @param x An object of class 'htest'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

reject_nullHyp <- function(x)
{
  if (!inherits(x, "htest"))
    stop("This function must be used with objects of class 'admix_test'")
  as.logical(x$reject_decision)
}


#' Extractor for object of class 'htest'
#'
#' Provide the selected rank of the test statistic (connected to the expansion order
#' of the densities in the orthonormal polynomial basis if method 'poly' was chosen;
#' or to the number of terms, i.e. discrepancies between couples of samples, included
#' in the test statistic).
#'
#' @param x An object of class 'htest'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

which_rank <- function(x)
{
  if (!inherits(x, "htest"))
    stop("This function must be used with objects of class 'htest'")
  as.numeric(x$selected_rank)
}


#' Extractor for object of class 'IBM_test'
#'
#' Provide the tabulated distribution that allows to define the extreme
#' quantile against which the test statistic is compared.
#'
#' @param x An object of class 'IBM_test'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_tabulated_dist <- function(x)
{
  if (!inherits(x, "IBM_test"))
    stop("This function must be used with objects of class 'IBM_test'")
  sort(x$tabulated_dist)
}


#' Extractor for object of class 'IBM_test'
#'
#' Provide the matrix storing the ranks of discrepancies using Inversion-Best Matching
#' approach between all couples among the K samples under study.
#'
#' @param x An object of class 'IBM_test'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_discrepancy_rank <- function(x)
{
  if (!inherits(x, "IBM_test"))
    stop("This function must be used with objects of class 'IBM_test'")
  x$discrepancy_rank
}


#' Extractor for object of class 'IBM_test'
#'
#' Provide the matrix storing the identifiers of discrepancies using Inversion-Best Matching
#' approach between all couples among the K samples under study.
#'
#' @param x An object of class 'IBM_test'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_discrepancy_id <- function(x)
{
  if (!inherits(x, "IBM_test"))
    stop("This function must be used with objects of class 'IBM_test'")
  x$discrepancy_id
}


#' Extractor for object of class 'IBM_test'
#'
#' Provide the terms (or discrepancies) that compose the test statistic
#' for the k-sample test.
#'
#' @param x An object of class 'IBM_test'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_statistic_components <- function(x)
{
  if (!inherits(x, "IBM_test"))
    stop("This function must be used with objects of class 'IBM_test'")
  x$statistic_name
}


#' Extractor for object of class 'admix_cluster'
#'
#' Extract the clusters that were discovered among K samples, where belonging to
#' one given cluster means having equal unknown component distributions.
#'
#' @param x An object of class 'admix_cluster'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_cluster_members <- function(x)
{
  if (!inherits(x, "admix_cluster"))
    stop("This function must be used with objects of class 'admix_cluster'")
  x$clusters
}


#' Extractor for object of class 'admix_cluster'
#'
#' Provide the number of samples in each cluster.
#'
#' @param x An object of class 'admix_cluster'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_cluster_sizes <- function(x)
{
  if (!inherits(x, "admix_cluster"))
    stop("This function must be used with objects of class 'admix_cluster'")
  x$clust_sizes
}


#' Extractor for object of class 'IBM_test' or 'admix_cluster'
#'
#' Provide the matrix storing discrepancies using Inversion-Best Matching
#' approach between all couples among the K samples under study.
#'
#' @param x An object of class 'IBM_test' or 'admix_cluster'.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

get_discrepancy_matrix <- function(x)
{
  if (!inherits(x, "admix_cluster") & !inherits(x, "IBM_test"))
    stop("This function must be used with objects of class 'IBM_test' or 'admix_cluster'")
  x$discrepancy_matrix
}
