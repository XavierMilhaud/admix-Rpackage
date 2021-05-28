#' Test for equality of the known components between two admixture models
#'
#' Test if the known components coming from the two two-components admixture models are the same.
#'
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#'
#' @return A boolean (TRUE if the known components are the same, otherwise FALSE).
#'
#' @examples
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                    f2 = list(mean = 2, sd = 0.3), g2 = list(mean = 0, sd = 1))
#' is_equal_knownComp(comp.dist = list.comp, comp.param = list.param)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

is_equal_knownComp <- function(comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (comp.dist[[2]] == comp.dist[[4]]) {
    vect.par <- FALSE
    if ( any(sapply(comp.param[[2]], length) != 1) | any(sapply(comp.param[[4]], length) != 1) ) vect.par <- TRUE
    if (!vect.par) {
      param_matrix <- as.matrix(sapply(comp.param[[2]], "==", comp.param[[4]]))
      if ( (all(diag(param_matrix) == TRUE)) & (nrow(param_matrix) != 0) & (ncol(param_matrix) != 0) ) {
        G1equalG2 <- TRUE
      } else {
        G1equalG2 <- FALSE
      }
    } else {
      G1equalG2 <- all(unlist(comp.param[[2]][which(sapply(comp.param[[2]],length) !=1 )]) == unlist(comp.param[[4]][which(sapply(comp.param[[4]],length) !=1 )]))
    }
  } else {
    G1equalG2 <- FALSE
  }

  return(G1equalG2)
}
