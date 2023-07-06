#' Distribution of the contrast in the Inversion - Best Matching (IBM) method
#'
#' Tabulate the distribution related to the inner convergence part of the contrast, by simulating trajectories of gaussian
#' processes and applying some transformations. Useful to perform the test hypothesis then, by retrieving the (1-alpha)-quantile
#' of interest. See 'Details' below and the cited paper therein for further information.
#'
#' @param n.sim Number of trajectories of simulated gaussian processes (number of random draws for tabulation).
#' @param n.varCovMat Number of time points on which gaussian processes are simulated.
#' @param sample1 Observations of the first sample under study.
#' @param sample2 Observations of the second sample under study.
#' @param min_size (default to NULL) In the k-sample case, useful to provide the minimal size among all samples. Otherwise, useless.
#' @param comp.dist A list with four elements corresponding to the component distributions (specified with R native names for these distributions)
#'                  involved in the two admixture models. The two first elements refer to the unknown and known components of the 1st admixture model,
#'                  and the last two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                  For instance, 'comp.dist' could be specified as follows: list(f1=NULL, g1='norm', f2=NULL, g2='norm').
#' @param comp.param A list with four elements corresponding to the parameters of the component distributions, each element being a list
#'                   itself. The names used in this list must correspond to the native R argument names for these distributions.
#'                   The two first elements refer to the parameters of unknown and known components of the 1st admixture model, and the last
#'                   two ones to those of the second admixture model. If there are unknown elements, they must be specified as 'NULL' objects.
#'                   For instance, 'comp.param' could be specified as follows: : list(f1=NULL, g1=list(mean=0,sd=1), f2=NULL, g2=list(mean=3,sd=1.1)).
#' @param parallel (default to FALSE) Boolean to indicate whether parallel computations are performed (speed-up the tabulation).
#' @param n_cpu (default to 2) Number of cores used for computations when parallelizing.
#'
#' @details See the paper presenting the IBM approach at the following HAL weblink: https://hal.science/hal-03201760
#'
#' @return A list with four elements, containing: 1) random draws of the quantity 'sample size times the empirical contrast',
#'         as defined in the IBM approach (see 'Details' above); 2) the estimated unknown component weights for the two admixture
#'         models; 3) the value of the quantity 'sample size times the empirical contrast'; 4) the sequence of points in the
#'         support that were used to evaluate the variance-covariance matrix of empirical processes.
#'
#' @examples
#' \donttest{
#' ## Simulate data:
#' list.comp <- list(f1 = 'norm', g1 = 'norm',
#'                   f2 = 'norm', g2 = 'norm')
#' list.param <- list(f1 = list(mean = 1, sd = 1), g1 = list(mean = 2, sd = 0.7),
#'                    f2 = list(mean = 1, sd = 1), g2 = list(mean = 3, sd = 1.2))
#' X.sim <- rsimmix(n=1000, unknownComp_weight=0.7, comp.dist = list(list.comp$f1,list.comp$g1),
#'                  comp.param = list(list.param$f1, list.param$g1))$mixt.data
#' Y.sim <- rsimmix(n=1200, unknownComp_weight=0.6, comp.dist = list(list.comp$f2,list.comp$g2),
#'                  comp.param = list(list.param$f2, list.param$g2))$mixt.data
#' ## Tabulate 1st term of stochastic integral (inner convergence) in a real-life setting:
#' list.comp <- list(f1 = NULL, g1 = 'norm',
#'                   f2 = NULL, g2 = 'norm')
#' list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
#'                    f2 = NULL, g2 = list(mean = 3, sd = 1.2))
#' U <- IBM_tabul_stochasticInteg(n.sim = 2, n.varCovMat = 20, sample1 = X.sim, sample2 = Y.sim,
#'                                min_size = NULL, comp.dist = list.comp, comp.param = list.param,
#'                                parallel = FALSE, n_cpu = 2)
#' }
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

IBM_tabul_stochasticInteg <- function(n.sim = 200, n.varCovMat = 100, sample1 = NULL, sample2 = NULL, min_size = NULL,
                                      comp.dist = NULL, comp.param = NULL, parallel = FALSE, n_cpu = 2)
{
  i <- NULL
  if (parallel) {
    `%fun%` <- foreach::`%dopar%`
    doParallel::registerDoParallel(cores = n_cpu)
  } else {
    `%fun%` <- foreach::`%do%`
  }

  stopifnot( (length(comp.dist) == 4) & (length(comp.param) == 4) )
  if (is.null(comp.dist[[2]]) | is.null(comp.dist[[4]]) | is.null(comp.param[[2]]) | is.null(comp.param[[4]])) {
    stop("Known components must be specified in the two admixture models.")
  }
  if ( any(sapply(comp.dist, is.null)) | any(sapply(comp.param, is.null)) ) {
    comp.dist[which(sapply(comp.dist, is.null) == TRUE)] <- "norm"
    comp.param[which(sapply(comp.param, is.null) == TRUE)] <- NA
    if (!all(unlist(sapply(comp.param, is.na)[c(1,3)]))) stop("Mixture distributions/parameters were not correctly specified")
  }

  exp.comp.dist <- paste0("p", comp.dist)
  if (any(exp.comp.dist == "pmultinom")) { exp.comp.dist[which(exp.comp.dist == "pmultinom")] <- "stepfun" }
#  comp_tab <- sapply(X = exp.comp.dist, FUN = get, pos = "package:stats", mode = "function")
  comp_tab <- sapply(X = exp.comp.dist, FUN = get, mode = "function")
  for (i in 1:length(comp_tab)) assign(x = names(comp_tab)[i], value = comp_tab[[i]])
  ## Create the expression involved in future assessments of the CDF:
  make.expr.step <- function(i) paste(names(comp_tab)[i],"(x = 1:", length(comp.param[[i]][[2]]), paste(", y = ", paste("cumsum(c(0,",
                                          paste(comp.param[[i]][[2]], collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  make.expr <- function(i) paste(names(comp_tab)[i],"(z,", paste(names(comp.param[[i]]), "=", comp.param[[i]], sep = "", collapse = ","), ")", sep="")
  expr <- vector(mode = "character", length = length(exp.comp.dist))
  expr[which(exp.comp.dist == "stepfun")] <- sapply(which(exp.comp.dist == "stepfun"), make.expr.step)
  expr[which(expr == "")] <- sapply(which(expr == ""), make.expr)
  expr <- unlist(expr)

  ## Differentiates when the 4 components were specified (simulations) and real-life cases (only known components are given).
  ## First, real-life case:
  if (all(unlist(sapply(comp.param, is.na)[c(1,3)]))) {
    if (any(exp.comp.dist == "stepfun")) {
      G1.fun <- eval(parse(text = expr[2]))
      G2.fun <- eval(parse(text = expr[4]))
      G1 <- function(z) G1.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      G1 <- function(z) { eval(parse(text = expr[2])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  } else {
    ## case of simulated data:
    if (any(exp.comp.dist == "stepfun")) {
      F1.fun <- eval(parse(text = expr[1]))
      G1.fun <- eval(parse(text = expr[2]))
      F2.fun <- eval(parse(text = expr[3]))
      G2.fun <- eval(parse(text = expr[4]))
      F1 <- function(z) F1.fun(z)
      G1 <- function(z) G1.fun(z)
      F2 <- function(z) F2.fun(z)
      G2 <- function(z) G2.fun(z)
    } else {
      F1 <- function(z) { eval(parse(text = expr[1])) }
      G1 <- function(z) { eval(parse(text = expr[2])) }
      F2 <- function(z) { eval(parse(text = expr[3])) }
      G2 <- function(z) { eval(parse(text = expr[4])) }
    }
  }
  ## Empirical cumulative distribution function from the two observed samples:
  L1 <- stats::ecdf(sample1)
  L2 <- stats::ecdf(sample2)

  ## Estimate the weights of the unknown component distributions in first and second samples:
  estim <- IBM_estimProp(sample1 = sample1, sample2 = sample2, known.prop = NULL, comp.dist = comp.dist, comp.param = comp.param,
                         with.correction = F, n.integ = 1000)
  if (is.null(min_size)) {
    sample.size <- min(length(sample1), length(sample2))
  } else {
    sample.size <- min_size
  }
  contrast_val <- sample.size *
    IBM_empirical_contrast(par = estim[["prop.estim"]], fixed.p.X = estim[["p.X.fixed"]], sample1 = sample1,
                          sample2 = sample2, G = estim[["integ.supp"]], comp.dist = comp.dist, comp.param = comp.param)
  ## Integration support:
  support <- detect_support_type(sample1, sample2)
  if (support == "continuous") {
    densite.G <- stats::density(estim[["integ.supp"]], bw = "SJ", adjust = 0.5, kernel = "gaussian")
    supp.integ <- c(min(estim[["integ.supp"]]), max(estim[["integ.supp"]]))
    t_seq <- seq(from = supp.integ[1], to = supp.integ[2], length.out = n.varCovMat)
  } else {
    ## Case of multinomial distribution :
    supp.integ <- estim[["integ.supp"]]
    t_seq <- sort(unique(supp.integ))
  }

  ## Compute the normalization matrix M(.) at each point, to be used further when determining the full simulated trajectories:
  normalization_factors <-
    foreach::foreach (i = 1:length(t_seq), .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
      IBM_normalization_term(t_seq[i], estim, estim[["p.X.fixed"]], NULL, sample1, sample2, min_size, comp.dist, comp.param)
    }

  ## Estimate the variance-covariance functions from the empirical processes:
#  cov_mat_L1 <- sapply(t_seq, function(s1) {
#                     sapply(t_seq, function(s2) {
#                       estimVarCov_empProcess(x = s1, y = s2, obs.data = sample1) })
#                     })
  cov_mat_L1 <- estimVarCov_empProcess_Rcpp(t = t_seq, obs_data = sample1)
#  cov_mat_L2 <- sapply(t_seq, function(s1) {
#                    sapply(t_seq, function(s2) {
#                      estimVarCov_empProcess(x = s1, y = s2, obs.data = sample2) })
#  })
  cov_mat_L2 <- estimVarCov_empProcess_Rcpp(t = t_seq, obs_data = sample2)

  ##------- Differentiates the cases where G1 = G2 or not --------##
  G1equalG2 <- is_equal_knownComp(comp.dist, comp.param)
  if (G1equalG2) {
    psi1 <- function(z) 2*( ((2-estim[["prop.estim"]])/estim[["prop.estim"]]^3) * G1(z) -
                              (2/estim[["prop.estim"]]^3)*L2(z) + (1/(estim[["prop.estim"]]^2*estim[["p.X.fixed"]]))*L1(z) -
                              ((1-estim[["p.X.fixed"]])/(estim[["prop.estim"]]^2*estim[["p.X.fixed"]])) * G1(z) )
    psi2 <- function(z) 2*( (1/(estim[["prop.estim"]]^2*estim[["p.X.fixed"]])) * (L2(z) - G1(z)) )
  } else {
    psi1.1 <- function(z) 2*( ((2-estim[["prop.estim"]][1])/estim[["prop.estim"]][1]^3) * G1(z) -
                                (2/estim[["prop.estim"]][1]^3)*L1(z) + (1/(estim[["prop.estim"]][1]^2*estim[["prop.estim"]][2]))*L2(z) -
                                ((1-estim[["prop.estim"]][2])/(estim[["prop.estim"]][1]^2*estim[["prop.estim"]][2])) * G2(z) )
    psi1.2 <- function(z) 2*( (1/(estim[["prop.estim"]][1]^2*estim[["prop.estim"]][2])) * (L1(z) - G1(z)) )
    psi2.1 <- function(z) 2*( ((2-estim[["prop.estim"]][2])/estim[["prop.estim"]][2]^3) * G2(z) -
                                (2/estim[["prop.estim"]][2]^3)*L2(z) + (1/(estim[["prop.estim"]][2]^2*estim[["prop.estim"]][1]))*L1(z) -
                                ((1-estim[["prop.estim"]][1])/(estim[["prop.estim"]][2]^2*estim[["prop.estim"]][1])) * G1(z) )
    psi2.2 <- function(z) 2*( (1/(estim[["prop.estim"]][2]^2*estim[["prop.estim"]][1])) * (L2(z) - G2(z)) )
  }

  U_sim <-
    foreach::foreach (i = 1:n.sim, .inorder = TRUE, .errorhandling = 'pass', .export = ls(globalenv())) %fun% {
      ## Simulate random gaussian processes using the variance-covariance functions determined from empirical processes:
      B1_path <- sim_gaussianProcess(mean_vec = rep(0,nrow(cov_mat_L1)), varCov_mat = cov_mat_L1, from = supp.integ[1],
                                     to = supp.integ[length(supp.integ)], start = 0, nb.points = nrow(cov_mat_L1))
      B1_traj <- stats::approxfun(x = t_seq, y = B1_path$traj1, method = "linear")
      B2_path <- sim_gaussianProcess(mean_vec = rep(0,nrow(cov_mat_L2)), varCov_mat = cov_mat_L2, from = supp.integ[1],
                                     to = supp.integ[length(supp.integ)], start = 0, nb.points = nrow(cov_mat_L2))
      B2_traj <- stats::approxfun(x = t_seq, y = B2_path$traj1, method = "linear")

      ## Introdude vector 'Z' :
      if (G1equalG2) {

        integrand.phi2 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi2(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi1 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi1(z) * B2_traj(z) * densite.G.dataPoint
        }
        Z <- matrix(NA, nrow = 4, ncol = length(t_seq))
        if (support == "continuous") {
          Z[1, ] <- stats::integrate(Vectorize(integrand.phi2), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[2, ] <- B1_path$traj1
          Z[3, ] <- stats::integrate(Vectorize(integrand.phi1), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[4, ] <- B2_path$traj1
        } else {
          Z[1, ] <- sum( unlist(sapply(supp.integ, integrand.phi2)) )
          Z[2, ] <- B1_path$traj1
          Z[3, ] <- sum( unlist(sapply(supp.integ, integrand.phi1)) )
          Z[4, ] <- B2_path$traj1
        }
        estim.vector <- matrix(NA, nrow = 2, ncol = length(t_seq))

      } else {

        integrand.phi1.1 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi1.1(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi2.2 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi2.2(z) * B1_traj(z) * densite.G.dataPoint
        }
        integrand.phi2.1 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi2.1(z) * B2_traj(z) * densite.G.dataPoint
        }
        integrand.phi1.2 <- function(z) {
          if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
          } else {
            densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
          }
          psi1.2(z) * B2_traj(z) * densite.G.dataPoint
        }
        Z <- matrix(NA, nrow = 6, ncol = length(t_seq))
        if (support == "continuous") {
          Z[1, ] <- stats::integrate(Vectorize(integrand.phi1.1), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[2, ] <- stats::integrate(Vectorize(integrand.phi2.2), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[3, ] <- B1_path$traj1
          Z[4, ] <- stats::integrate(Vectorize(integrand.phi2.1), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[5, ] <- stats::integrate(Vectorize(integrand.phi1.2), lower=supp.integ[1], upper=supp.integ[2], subdivisions = 1000L, rel.tol = 1e-3)$value
          Z[6, ] <- B2_path$traj1
        } else {
          Z[1, ] <- sum( unlist(sapply(supp.integ, integrand.phi1.1)) )
          Z[2, ] <- sum( unlist(sapply(supp.integ, integrand.phi2.2)) )
          Z[3, ] <- B1_path$traj1
          Z[4, ] <- sum( unlist(sapply(supp.integ, integrand.phi2.1)) )
          Z[5, ] <- sum( unlist(sapply(supp.integ, integrand.phi1.2)) )
          Z[6, ] <- B2_path$traj1
        }
        estim.vector <- matrix(NA, nrow = 3, ncol = length(t_seq))
      }

      ## Get the trajectory of function 'sqrt(n)*(Dn(.)-D(.))' :
      for (j in 1:length(t_seq)) estim.vector[ ,j] <- normalization_factors[[j]] %*% Z[ ,j]
      D_function <- stats::approxfun(x = t_seq, y = estim.vector[nrow(estim.vector), ], method = "linear")

      integrand <- function(z) {
        if (support == "continuous") { densite.G.dataPoint <- stats::approx(densite.G$x, densite.G$y, xout = z)$y
        } else {
          densite.G.dataPoint <- 1 / length(table(c(sample1,sample2)))
        }
        D_function(z)^2 * densite.G.dataPoint
      }

      if (support == "continuous") {
        U <- stats::integrate(integrand, lower = supp.integ[1], upper = supp.integ[2], subdivisions = 10000L, rel.tol = 1e-04)$value
      } else {
        U <- sum( unlist(sapply(supp.integ, integrand)) )
      }
      U
    }

  indexes.toRemove <- which( (substr(U_sim, start = 1, stop = 5) == "Error") | (substr(U_sim, start = 1, stop = 5) == "list(") )
  if (length(indexes.toRemove) != 0) U_sim <- U_sim[-indexes.toRemove]

  return( list(U_sim = as.numeric(U_sim), estimator = estim, contrast_value = contrast_val, integ.points = t_seq) )

}
