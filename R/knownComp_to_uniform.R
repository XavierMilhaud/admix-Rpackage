#' Transforms the known component of the admixture distribution to a Uniform distribution
#'
#' In admixture such that the probability density function (pdf) follows l = p*f + (1-p)*g, where p is the unknown
#' weight and f is the unknown component distribution: transforms g of the two-component mixture ditribution to a
#' Uniform distribution. Useful to use Patra and Sen estimator for the estimation of the unknown weight p.
#'
#' @param data Observations of the sample under study, following an admixture distribution.
#' @param comp.dist A list with two elements corresponding to component distributions (specified with R native names for these distributions) involved
#'                  in the admixture model. Unknown elements must be specified as 'NULL' objects, e.g. when 'f' is unknown: list(f=NULL, g='norm').
#' @param comp.param A list with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   Unknown elements must be specified as 'NULL' objects, e.g. if 'f' is unknown: list(f=NULL, g=list(mean=0,sd=1)).
#'
#' @return The transformed data, i.e. the transformed mixture distribution where the known component g now follows a
#'         Uniform(0,1) distribution.
#'
#' @examples
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm')
#' list.param <- list(f1 = list(mean = 3, sd = 0.5),
#'                    g1 = list(mean = 0, sd = 1))
#' sample1 <- rsimmix(n=1500, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
#'                                                    comp.param=list(list.param$f1,list.param$g1))
#' plot_admix(sim.X = sample1[['mixt.data']], support = 'continuous')
#' ## Transform the known component into a Uniform(0,1) distribution:
#' list.comp <- list(f1 = NULL, g1 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1))
#' transformed_data <- knownComp_to_uniform(data = sample1[['mixt.data']],
#'                                          comp.dist = list.comp, comp.param = list.param)
#' plot_admix(sim.X = transformed_data, support = 'continuous')
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

knownComp_to_uniform <- function(data, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")

  ## Extracts the information on component distributions for inversion (transformation to uniform distribution of the known component):
  comp.dist.inv <- paste0("p", comp.dist[[2]])
  comp.inv <- sapply(X = comp.dist.inv, FUN = get, pos = "package:stats", mode = "function")
  assign(x = names(comp.inv)[1], value = comp.inv[[1]])

  ## Check if arguments of R core functions were correctly specified:
  arg.names <- sapply(X = comp.inv, FUN = methods::formalArgs)
  n.arg.user <- length(comp.param[[2]])
  if (inherits(arg.names, what = "matrix")) {
    arg.names.supplied <- names(comp.param[[2]])
    if (is.character(arg.names.supplied)) {
      common.args <- match(x = arg.names.supplied, table = arg.names)
    } else {
      common.args <- apply(sapply(comp.param, names), 2, match, table = arg.names)
    }
    if (any(is.na(common.args))) stop("Parameters of the mixture components were not correctly specified")
  } else {
    common.args <- vector(mode = "list", length = length(comp.dist))
    for (i in 1:length(comp.dist)) {
      if (class(sapply(comp.param, names))[1] != "matrix") {
        common.args[[i]] <- sapply(sapply(comp.param, names)[[i]], match, arg.names[[i]])
      } else { common.args[[i]] <- match(sapply(comp.param, names)[, i], arg.names[[i]])  }
      if (any(is.na(common.args[[i]]))) stop("Parameters of the mixture components were not correctly specified")
    }
  }

  ## Inversion of the second component to get a Uniform distribution for the second component:
  ## Creates the expression allowing further to generate the right data:
  make.expr.inv <- function(z) paste(names(comp.inv)[1],"(q=", z, ",", paste(names(comp.param[[2]]), "=", comp.param[[2]], sep="", collapse=","), ")", sep="")
  expr.inv <- parse(text  = make.expr.inv(data))
  data.transformed <- sapply(expr.inv, eval)

  return(data.transformed)
}
