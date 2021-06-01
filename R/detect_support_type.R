#' Detect the support of the random variables under study
#'
#' Given one or two sets of observations (two samples), the function provides with the most plausible type of support for the
#' underlying random variables to be studied. Basically, if less than 3 percent of the observations have different values,
#' we consider that the support is discrete. Otherwise, we consider it as a continuous support.
#'
#' @param sample1 The first sample of observations under study.
#' @param sample2 The second sample of observations under study.
#'
#' @return The type of support, either discrete or continuous.
#'
#' @examples
#' ## Simulate the two mixture samples:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
#'                     f2 = list(mean = 1, sd = 0.1), g2 = list(mean = 5, sd = 2))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                    comp.param=list(list.param$f1,list.param$g1))
#' sample2 <- rsimmix(n=2000, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
#'                    comp.param=list(list.param$f2,list.param$g2))
#' ## Test the type of support:
#' detect_support_type(sample1[['mixt.data']], sample2[['mixt.data']])
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

detect_support_type <- function(sample1, sample2 = NULL)
{
  if (is.null(sample2)) {
    if ((length(unique(sample1)) / length(sample1)) < 0.03) {
      support <- "discrete"
    }  else {
      support <- "continuous"
    }
  } else {
    if ( ((length(unique(sample1)) / length(sample1)) < 0.03) & ((length(unique(sample2)) / length(sample2)) < 0.03) ) {
      support <- "discrete"
    }  else {
      support <- "continuous"
    }
  }

  return(support)
}
