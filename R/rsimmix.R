#' Simulation of a two-component mixture model
#'
#' Simulate a two-component mixture model following the probability density function (pdf) l such that l = p*f + (1-p)*g,
#' with f and g mixture component distributions, and p the mixture weight.
#'
#' @param n Number of observations to be drawn.
#' @param unknownComp_weight Weight of the component distribution f (representing the unknown component in admixture models).
#' @param comp.dist A list with two elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the mixture model. These elements respectively refer to the two components f and g.
#'                  No unknown elements permitted. For instance, 'comp.dist' could be set equal to list(f = 'rnorm', g = 'norm').
#' @param comp.param A list with two elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   These elements respectively refer to the parameters of f and g distributions of the mixture model. No unknown elements permitted.
#'                   For instance, 'comp.param' could be set equal to list(f=list(mean=2,sd=0.3), g=list(mean=0,sd=1)).
#'
#' @return A list of three components. The first, named 'mixt.data', is the simulated sample from the specified mixture distribution.
#'         The second, named 'unknown.data', refers to the data simulated corresponding to the distribution f. The third, named 'known.data',
#'         corresponds to the observations affiliated to the known component g.
#'
#' @examples
#' sim.X <- rsimmix(n = 2000, unknownComp_weight = 0.7, comp.dist = list(f = 'norm', g = 'norm'),
#'                  comp.param = list(f = list(mean = 3, sd = 0.5), g = list(mean = 0, sd = 1)))
#' class(sim.X)
#' attributes(sim.X)
#' plot_admix(sim.X = sim.X$mixt.data, sim.Y = NULL, user.bounds = NULL, support = 'continuous')
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

rsimmix <- function(n = 1000, unknownComp_weight = 0.5,
                    comp.dist = list(f = 'norm', g = 'norm'),
                    comp.param = list(f = c(mean=0, sd=1), g = c(mean=2, sd=1)))
{
  ## Some arguments to check:
  if ( !((unknownComp_weight > 0) & (unknownComp_weight < 1)) ) stop("Weight of the unknown component is not appropriate and should belong to ]0,1[.")
  if (length(comp.dist) != 2) stop("Please provide two distributions, one for each of the mixture components.")
  if (comp.dist[[1]] == "multinom") {
    if (comp.dist[[2]] != "multinom") stop("Multinomial distribution component is not supported when mixed with other than another multinomial distribution")
  }

  ## Extracts the information on component distributions:
  comp.dist <- paste0("r", comp.dist)
  #comp_sim <- sapply(X = comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_sim <- sapply(X = comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_sim)) assign(x = names(comp_sim)[i], value = comp_sim[[i]])

  ## Check if arguments of R core functions were correctly specified:
  arg.names <- sapply(X = comp_sim, FUN = methods::formalArgs)
  n.arg.user <- sapply(X = comp.param, FUN = length)
  if (inherits(arg.names, what = "matrix")) {
    arg.names.supplied <- sapply(comp.param, names)
    if (is.character(arg.names.supplied)) {
      common.args <- match(x = arg.names.supplied, table = arg.names)
    } else {
      common.args <- apply(sapply(comp.param, names), 2, match, table = arg.names)
    }
    if (any(is.na(common.args))) stop("Parameters of the mixture components were not correctly specified")
  } else {
    common.args <- vector(mode = "list", length = length(comp.dist))
    for (i  in 1:length(comp.dist)) {
      if (class(sapply(comp.param, names))[1] != "matrix") {
        common.args[[i]] <- sapply(sapply(comp.param, names)[[i]], match, arg.names[[i]])
      } else { common.args[[i]] <- match(sapply(comp.param, names)[, i], arg.names[[i]])  }
      if (any(is.na(common.args[[i]]))) stop("Parameters of the mixture components were not correctly specified")
    }
  }

  ## Creates the expression allowing further to generate the right data:
  make.expr_sim <- function(i) paste(names(comp_sim)[i], "(1,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr_sim <- sapply(1:length(comp.dist), make.expr_sim)
  ## Generates the label for each observation:
  z <- sample(x = 2, size = n, replace = TRUE, prob = c(unknownComp_weight,1-unknownComp_weight))
  ## Generates the final mixture data:
  data.gen <- parse(text = expr_sim[z])
  res <- sapply(data.gen, eval)

  if (any(comp.dist == "rmultinom")) {
    res_tmp <- rowSums(res)
    res <- unlist( apply( as.data.frame(1:length(res_tmp)), 1, function(k) { rep(k, res_tmp[k]) } ) )
  }

  ## Catch the information whether the observation comes from the first or the second component of the mixture distribution:
  f.obs <- res[z == 1]           # observations corresponding to the unknown component f
  g.obs <- res[z == 2]           # observations corresponding to the known component g

  return( list(mixt.data = res, unknown.data = f.obs, known.data = g.obs) )
}


#' Simulation of a two-component mixture with one component following a two-component mixture
#'
#' simulate a two-component admixture model, where the first component is a mixture itself
#'
#' @param n is the number of observations to be drawn
#' @param a the shift of the mean for the two distributions that are embedded in the unknown component
#' @param m the mean (up to the shift a) of the unknown components
#' @param s the standard deviation of the unknown components
#' @param p the weight of the unknown component (itself a mixture).
#'
#' @return a list containing the data generated from a mixture of mixture distribution, the data where the known component density has
#'         been made uniform(0,1), and the known data (corresponding to the part of data generated from the known component density).
#'
#' @examples
#' sample1 <- rsimmix_mix(n = 3000, m = 5, s = 0.5, p = 0.3, a = 2)[['mixt.data']]
#' plot(stats::density(sample1))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'
rsimmix_mix <- function(n, m, s, p, a) {
  z <- stats::rbinom(n, 1, p)
  xx <- stats::rnorm(n, 0, 1)
  z2 <- stats::rbinom(n, 1, 0.5)
  x2 <- stats::rnorm(n, m, s)
  w <- (1-z2) * (x2+a) + z2 * (x2-a)
  simul <- (1-z) * xx + z * (w)
  data.transformees <- stats::pnorm(simul, mean = 0, sd = 1)
  return( list(mixt.data = simul, data.transform = data.transformees, known.data = xx) )
}
