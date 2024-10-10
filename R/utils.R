# .onLoad <- function(libname, pkgname)
# {
#   library.dynam("admix", pkgname, libname)
# }

admixStartupMessage <- function()
{
  msg <- c(paste0(
    "This is package admix, version ",
    utils::packageVersion("admix")),
    "\n-------------------------------\n",
    "Type 'citation(\"admix\")' for citing this R package in publications.",
    "\n-------------------------------\n",
    "This work was partly conducted within the Research Chair DIALog under the aegis of the Risk Foundation, an initiative by CNP Assurances.\n")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # unlock .admix variable allowing its modification
  #unlockBinding(".admix", asNamespace("admix"))
  # startup message
  msg <- admixStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'admix' version", packageVersion("admix"))
  base::packageStartupMessage(msg)
  base::invisible()
}

#' Detect the support of some random variables
#'
#' Given one or two sets of observations (two samples), the function provides with the most plausible type of support for the
#' underlying random variables to be studied. If less than 3 percents of the observations have different values,
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
#' @noRd

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


#' Simulation of a Gaussian process
#'
#' Simulate the trajectory of a Gaussian process, given a mean vector and a variance-covariance structure.
#'
#' @param mean_vec Vector (if multimensional) of means for the increments following gaussian distribution.
#' @param varCov_mat Corresponding variance-covariance structure.
#' @param from Initial time point at which the process is simulated.
#' @param to Last time point at which the process is simulated.
#' @param start Useful if the user wants to make the trajectory start from some given value.
#' @param nb.points Number of points at which the process is simulated.
#'
#' @return The trajectory of the Gaussian processes after simulating the multivariate Gaussian distributions with
#'         specified variance-covariance structure.
#'
#' @examples
#' list.comp <- list(f1 = "norm", g1 = "norm")
#' list.param <- list(f1 = list(mean = 12, sd = 0.4),
#'                    g1 = list(mean = 16, sd = 0.7))
#' sample1 <- rsimmix(n = 2000, unknownComp_weight = 0.5, comp.dist = list.comp,
#'                    comp.param = list.param)$mixt.data
#' ## First get the variance-covariance matrix of the empirical process (Donsker correlation):
#' cov_mat <- .Call('_admix_estimVarCov_empProcess_Rcpp', PACKAGE = 'admix',
#'                   seq(from = min(sample1), to = max(sample1), length.out = 100), sample1)
#' ## Plug it into the simulation of the gaussian process:
#' B1 <- sim_gaussianProcess(mean_vec=rep(0,nrow(cov_mat)), varCov_mat=cov_mat, from=min(sample1),
#'                           to = max(sample1), start = 0, nb.points = nrow(cov_mat))
#' plot(x = B1$dates, y = B1$traj1, type="l", xlim = c(min(sample1),max(sample1)), ylim = c(-1,1))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

sim_gaussianProcess <- function(mean_vec, varCov_mat, from = 0, to = 1, start = 0, nb.points = 10)
{
  t <- seq(from = from, to = to, length.out = nb.points)
  rownames(varCov_mat) <- paste("x=", round(t,2), sep = "")
  colnames(varCov_mat) <- paste("y=", round(t,2), sep = "")

  ## Simulation of the multivariate gaussian distributions, without trend :
  trajectory.L1 <- MASS::mvrnorm(n = 1, m = mean_vec, Sigma = varCov_mat)
  traj1 <- trajectory.L1 - trajectory.L1[1] + start

  return( list(dates = t, traj1 = traj1) )
}


#' Equality of known components in 2 admixture models
#'
#' Test if the known component distributions coming from 2 two-components admixtures are identical.
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
#' @noRd

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


#' Kernel estimation of some CDF
#'
#' Functions to perform the estimation of cumulative distribution function (CDF) by kernel estimators
#' (with a non-gaussian kernel).
#'
#' @param h window of the kernel estimation.
#' @param u the point at which the estimation is made.
#'
#' @return the estimated value of the cdf.
#'
#' @examples
#' kernel_cdf(0.4,0.5)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

kernel_cdf <- function(u, h)
{
  return( (0.5 + u/h + (u/h)^2 / 2) * ((u > -h) & (u <= 0)) +
            (0.5 + u/h - (u/h)^2 / 2) * ((0 < u) & (u < h)) + 1*(u >= h) )
}


#' Kernel estimation of some density function
#'
#' Functions to perform the estimation of probability density function (pdf) by kernel estimators (with a non-gaussian kernel).
#'
#' @param h window of the kernel estimation.
#' @param u the point at which the estimation is made.
#'
#' @return the estimated value of the pdf.
#'
#' @examples
#' kernel_density(0.4,0.5)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

kernel_density <- function(u, h)
{
  return( (1 - abs(u/h)) * ((u > -h) & (u < h)) / h )
}


#' Transforms the known component distribution to a Uniform
#'
#' In an admixture such that the probability density function (pdf) follows l = p*f + (1-p)*g, where p is the unknown
#' weight and f is the unknown component distribution: transforms the distribution g to a Uniform distribution.
#' Useful to use Patra and Sen estimator for the estimation of the unknown weight p.
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
#' plot_mixt_density(samples = list(sample1[['mixt.data']]), support = 'continuous')
#' ## Transform the known component into a Uniform(0,1) distribution:
#' list.comp <- list(f1 = NULL, g1 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1))
#' transformed_data <- knownComp_to_uniform(data = sample1[['mixt.data']],
#'                                          comp.dist = list.comp, comp.param = list.param)
#' plot_mixt_density(samples = list(transformed_data), support = 'continuous')
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

knownComp_to_uniform <- function(data, comp.dist, comp.param)
{
  stopifnot( (length(comp.dist) == 2) & (length(comp.param) == 2) )
  if (is.null(comp.dist[[2]]) | is.null(comp.param[[2]])) stop("Known component must be specified.")

  ## Extracts the information on component distributions for inversion (transformation to uniform distribution of the known component):
  comp.dist.inv <- paste0("p", comp.dist[[2]])
  #  comp.inv <- sapply(X = comp.dist.inv, FUN = get, pos = "package:stats", mode = "function")
  comp.inv <- sapply(X = comp.dist.inv, FUN = get, mode = "function")
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


#' Builds a polynomial orthonormal basis
#'
#' Builds a polynomial orthonormal basis, needed to decompose the probability density function (pdf) of the unknown component
#' from the admixture, depending on the support under consideration.
#'
#' @param support Support of the random variables implied in the two-component mixture distribution.
#' @param deg Degree up to which the basis is built.
#' @param x (NULL by default) Only used when support is 'Integer'. The point at which the polynomial value will be evaluated.
#' @param m (NULL by default) Only used when support is 'Integer'. Corresponds to the mean of the reference measure, i.e. Poisson(m).
#'
#' @return the orthonormal polynomial basis used to decompose the density of the unknown component of the mixture distribution.
#'
#' @examples
#' poly_orthonormal_basis(support = 'Real', deg = 10, x = NULL, m = NULL)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#' @noRd

poly_orthonormal_basis <- function(support = c("Real","Integer","Positive","Bounded.continuous","Bounded.discrete"), deg, x, m)
{
  if (support == "Real") {
    ## Hermite polynomials :
    poly_basis <- orthopolynom::hermite.he.polynomials(n = deg, normalized = FALSE)

  } else if (support == "Integer") {
    ## Charlier polynomials : 'm' is the mean of the reference measure, i.e. P(m).
    poly_basis <- function(x, m, deg) {
      CH <- CHA <- matrix(NA, nrow = length(x), ncol = deg)
      CH[ ,1] <- (m-x) / m
      CH[ ,2] <- ((m+1-x) * (m-x) / m - 1) / m
      for (k in 3:deg) { CH[ ,k] <- ((m+k-1-x) * CH[ ,k-1] - (k-1) * CH[ ,k-2]) / m }
      for (k in 1:deg) { CHA[ ,k] <- CH[ ,k] / (sqrt(factorial(k)*m^{-k})) }
      return(CHA)
    }

  } else if (support == "Positive") {
    ## Laguerre polynomials :
    poly_basis <- orthopolynom::laguerre.polynomials(n = deg, normalized = FALSE)

  } else if (support == "Bounded.continuous") {
    ## Legendre polynomials :
    poly_basis <- orthopolynom::legendre.polynomials(n = deg, normalized = FALSE)

  } else stop("Please give a correct argument for the support of the distributions!")

  return(poly_basis)
}
