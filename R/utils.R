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

#' Detect the type of support of some random variables
#'
#' Given one or two sets of observations (samples), the function provides with the most plausible type of support for the
#' underlying random variables to be studied. If less than 3 percents of the observations have different values,
#' we consider that the support is discrete. Otherwise, we consider it as a continuous support.
#'
#' @param sample1 The first sample of observations under study.
#' @param sample2 The second sample of observations under study.
#'
#' @return The type of support, either discrete or continuous.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1500, weight = 0.5,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' mixt2 <- twoComp_mixt(n = 2000, weight = 0.7,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 5, "sd" = 2)))
#' data2 <- getmixtData(mixt2)
#' ## Test the type of support:
#' detect_support_type(data1, data2)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

detect_support_type <- function(sample1, sample2 = NULL)
{
  if (is.null(sample2)) {
    if ((length(unique(sample1)) / length(sample1)) < 0.03) { support <- "Discrete"
    } else { support <- "Continuous" }
  } else {
    if ( ((length(unique(sample1)) / length(sample1)) < 0.03) &
         ((length(unique(sample2)) / length(sample2)) < 0.03) ) { support <- "Discrete"
    } else { support <- "Continuous" }
  }
  return(support)
}


#' Simulation of a Gaussian process
#'
#' Simulate the trajectory of a Gaussian process, given a mean vector and a variance-covariance structure.
#'
#' @param mean_vec Vector (if multimensional) of means for the increments following Gaussian distribution.
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
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 450, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## First get the variance-covariance matrix of the empirical process (Donsker correlation):
#' cov_mat <- .Call('_admix_estimVarCov_empProcess_Rcpp', PACKAGE = 'admix',
#'                   seq(from = min(data1), to = max(data1), length.out = 100), data1)
#' ## Plug it into the simulation of the gaussian process:
#' B1 <- sim_gaussianProcess(mean_vec=rep(0,nrow(cov_mat)), varCov_mat=cov_mat, from=min(data1),
#'                           to = max(data1), start = 0, nb.points = nrow(cov_mat))
#' plot(x = B1$dates, y = B1$traj1, type="l", xlim = c(min(data1),max(data1)), ylim = c(-1,1))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
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
#' @param admixMod1 An object of class 'admix_model' for the first admixture model.
#' @param admixMod2 An object of class 'admix_model' for the second admixture model.
#'
#' @return A boolean (TRUE if the known components are the same, otherwise FALSE).
#'
#' @examples
#' admixMod1 <- admix_model(knownComp_dist = "norm",
#'                          knownComp_param = list("mean"=0, "sd"=1))
#' admixMod2 <- admix_model(knownComp_dist = "norm",
#'                          knownComp_param = list("mean"=0, "sd"=1))
#' is_equal_knownComp(admixMod1, admixMod2)
#'
#' admixMod1 <- admix_model(knownComp_dist = "multinom",
#'                          knownComp_param = list("size"=1, "prob"=c(0.2,0.5,0.3)))
#' admixMod2 <- admix_model(knownComp_dist = "multinom",
#'                          knownComp_param = list("size"=1, "prob"=c(0.2,0.4,0.4)))
#' is_equal_knownComp(admixMod1, admixMod2)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

is_equal_knownComp <- function(admixMod1, admixMod2)
{
  if (!inherits(x = admixMod1, what = "admix_model"))
    stop("Argument 'admixMod1' is not correctly specified. See ?admix_model.")
  if (!inherits(x = admixMod2, what = "admix_model"))
    stop("Argument 'admixMod2' is not correctly specified. See ?admix_model.")

  if ( (admixMod1$comp.dist$known == admixMod2$comp.dist$known) &
       all(names(admixMod1$comp.param$known) == names(admixMod2$comp.param$known)) &
       all(unlist(admixMod1$comp.param$known) == unlist(admixMod2$comp.param$known)) ) {
    G1equalG2 <- TRUE
  } else G1equalG2 <- FALSE
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
#' @noRd

kernel_density <- function(u, h)
{
  return( (1 - abs(u/h)) * ((u > -h) & (u < h)) / h )
}


#' Transforms the known component distribution to a Uniform distribution
#'
#' In an admixture such that the probability density function (pdf) follows l = p*f + (1-p)*g, where p is the unknown
#' weight and f is the unknown component distribution: transforms the known distribution g to a Uniform distribution.
#' Useful when dealing with the Patra and Sen estimator (for the estimation of the unknown weight p).
#'
#' @param data Observations of the sample under study, following an admixture distribution.
#' @param admixMod An object of class 'admix_model', containing useful information about the known components and their parameter(s).
#'
#' @return The transformed data, i.e. the transformed mixture distribution where the known component g now follows a
#'         Uniform(0,1) distribution.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 1500, weight = 0.5,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' plot(mixt1)
#' data1 <- getmixtData(mixt1)
#'
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Transform the known component into a Uniform(0,1) distribution:
#' transformed_data <- knownComp_to_uniform(data = data1, admixMod = admixMod1)
#' plot(density(transformed_data))
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

knownComp_to_uniform <- function(data, admixMod)
{
  if (!inherits(x = admixMod, what = "admix_model"))
    stop("Argument 'admixMod' is not correctly specified. See ?admix_model.")

  ## Extracts the information about component distributions for inversion
  ## (transformation to uniform distribution of the known component):
  comp.dist.inv <- paste0("p", admixMod$comp.dist$known)
  if (comp.dist.inv == "pmultinom") eff <- as.numeric(table(data))
  comp.inv <- sapply(X = comp.dist.inv, FUN = get, mode = "function")
  assign(x = names(comp.inv)[1], value = comp.inv[[1]])

  ## Creates the adequate expression:
  make.expr.multinom <- function(z) {
    paste(names(comp.inv)[1],"(q = eff, n = sum(eff), prob = ",
          paste("c(", paste(admixMod$comp.param$known$prob, collapse = ","), "), lower.tail=TRUE)", sep = ""), sep = "")
  }
  #make.expr.multinom <- function(z) {
  #  paste(names(comp.inv)[1],"(q=c(rep(0,", z-1, "), 1, rep(0,", length(admixMod$comp.param$known$prob)-z, ")), n = 1, ",
  #        paste("c(", paste(admixMod$comp.param$known$prob, collapse = ","), "), lower.tail=TRUE)", sep = ""), sep = "")
  #}
  make.expr.inv <- function(z) paste(names(comp.inv)[1],"(q=", z, ",", paste(names(admixMod$comp.param$known),
                                     "=", admixMod$comp.param$known, sep="", collapse=","), ")", sep="")
  if (comp.dist.inv == "pmultinom") {
    expr.inv <- parse(text = make.expr.multinom(data))
  } else {
    expr.inv <- parse(text = make.expr.inv(data))
  }

  ## Inversion of the second component to get a Uniform distribution for the second component:
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
#' poly_orthonormal_basis(support = "Real", deg = 10, x = NULL, m = NULL)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
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


#' Expansion coefficients of some given density in an orthonormal polynomial basis.
#'
#' Compute the coefficients when decomposing some density in a given orthonormal polynomial basis.
#'
#' @param data Observed sample from which the coefficients are calculated.
#' @param supp Support of the density considered.
#' @param degree Degree up to which the polynomial basis is built.
#' @param m (default to 3) Only used when support is 'Integer'. Corresponds to the mean of the reference measure, i.e. Poisson(m).
#' @param other (default to NULL) A list to precise bounds when the support is bounded, where the second and fourth elements give bounds.
#'
#' @return The list composed of 'degree' elements, each element being a numeric vector (with sample size) where each value represents
#'         the k-th order coefficient found when decomposing the density in the orthonormal polynomial basis.
#'
#' @examples
#' ## Simulate data:
#' sample1 <- rnorm(n = 7000, mean = 3, sd = 1)
#' ## Compute the expansion coefficients in the orthonormal polynomial basis:
#' coeff <- orthoBasis_coef(data = sample1, supp = "Real", degree = 3, m=NULL, other=NULL)
#' sapply(coeff, mean)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

orthoBasis_coef <- function(data, supp = c('Real','Integer','Positive','Bounded.continuous'),
                            degree, m = 3, other = NULL)
{
  ## Builds the orthonormal polynomial basis:
  if (supp == "Integer") {
    stop("Still to implement in 'orthoBasis_coef'")
  } else {
    poly_basis <- poly_orthonormal_basis(support = supp, deg = degree, x = NULL, m = NULL)
  }

  ## Store the coefficients to be computed :
  coef.list <- vector(mode = "list", length = degree)

  ## Depending on the support :
  if (supp == "Real") {
    ## Reference measure N(0,1)
    for (i in 1:degree) coef.list[[i]] <- orthopolynom::polynomial.values(poly_basis, data)[[i+1]] / sqrt(factorial(i))

  } else if (supp == "Integer") {
    ## Reference measure P(3)
    for (i in 1:degree) coef.list[[i]] <- poly_basis[ ,i]

  } else if (supp == "Positive") {
    ## Reference measure Exp(1)
    for (i in 1:degree) coef.list[[i]] <- orthopolynom::polynomial.values(poly_basis, data)[[i+1]]

  } else if (supp == "Bounded.continuous") {
    ## Reference measure Unif(a,b)
    if (is.null(other)) { bounds <- c(min(data), max(data))
    } else { bounds <- other[[2]] }
    for (i in 1:degree) coef.list[[i]] <- (orthopolynom::polynomial.values(poly_basis, (2*data-bounds[1]-bounds[2])/(bounds[2]-bounds[1]))[[i+1]]) / sqrt(2*i+1)

  } else stop("Change the support since the choosen one is not considered!")

  return(coef.poly = coef.list)
}


#' Variance-covariance matrix of the empirical process in an admixture model
#'
#' Estimate the variance-covariance matrix of some given empirical process, based on the Donsker correlation.
#' Compute Donsker correlation between two time points (x,y) for some given empirical process with R code
#' (another implementation in C++ is also available to speed up this computation).
#'
#' @param x First time point considered for the computation of the correlation given the empirical process.
#' @param y Second time point considered for the computation of the correlation given the same empirical process.
#' @param obs.data Sample that permits to estimate the cumulative distribution function (cdf).
#'
#' @return The estimated variance-covariance matrix.
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 2500, weight = 0.5,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 12, "sd" = 0.4),
#'                                         list("mean" = 16, "sd" = 0.7)))
#' data1 <- getmixtData(mixt1)
#' ## Compute the variance-covariance matrix of the corresponding empirical process:
#' t <- seq(from = min(data1), to = max(data1), length = 50)
#' S2 <- sapply(t, function(s1) {
#'                 sapply(t, function(s2) {
#'                      estimVarCov_empProcess(x = s1, y = s2, obs.data = data1) })
#'                 })
#' lattice::wireframe(S2)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @noRd

estimVarCov_empProcess <- function(x, y, obs.data)
{
  L.CDF <- stats::ecdf(obs.data)
  ## Use the Donsker correlation formula:
  res <- L.CDF(min(x,y)) * (1 - L.CDF(max(x,y)))
  return(res)
}
