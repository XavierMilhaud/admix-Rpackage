#' Probability density function of the unknown component
#'
#' Estimates the decontaminated probability density function (PDF) of the unknown component in
#' an admixture model, based on the inversion of the admixture density equation \eqn{l = p f + (1-p) g}.
#'
#' @param sample1 Numeric vector, sample under study.
#' @param admixMod An object of class \code{admix_model}, containing useful information about known distribution(s) and parameter(s).
#' @param estim.p Numeric. The estimated mixing proportion \eqn{\hat{p}} of the unknown component.
#'
#' @details
#' The decontaminated density \eqn{f} is computed as:
#' \deqn{f(x) = (1 / \hat{p}) [ \hat{l}(x) - (1 - \hat{p}) g(x) ]}
#' where:
#' \itemize{
#'   \item \eqn{\hat{l}(x)} is the empirical density of the sample,
#'   \item \eqn{g(x)} is the known component’s theoretical density,
#'   \item \eqn{\hat{p}} is the estimated mixture weight.
#' }
#' For continuous data, \eqn{\hat{l}(x)} is estimated using kernel density estimation.
#' For discrete data, it is approximated from normalized frequencies.
#'
#' @return An object of class \code{decontaminated_density} containing:
#' \describe{
#'   \item{data}{Original sample.}
#'   \item{support}{Type of support ("Continuous" or "Discrete").}
#'   \item{decontaminated_density_fun}{A function returning the estimated decontaminated density.}
#' }
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 400, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1), admixMod = list(admixMod1),
#'                    est_method = 'PS')
#' ## Determine the decontaminated version of the unknown density by inversion:
#' x <- decontaminated_density(sample1 = data1, admixMod = admixMod1,
#'                             estim.p = get_mixing_weights(est))
#' print(x)
#' summary(x)
#' plot(x)
#'
#' ####### Discrete support:
#' mixt1 <- twoComp_mixt(n = 4000, weight = 0.6,
#'                       comp.dist = list("pois", "pois"),
#'                       comp.param = list(list("lambda" = 3),
#'                                         list("lambda" = 2)))
#' mixt2 <- twoComp_mixt(n = 3000, weight = 0.8,
#'                       comp.dist = list("pois", "pois"),
#'                       comp.param = list(list("lambda" = 3),
#'                                         list("lambda" = 4)))
#' mixt3 <- twoComp_mixt(n = 1500, weight = 0.5,
#'                       comp.dist = list("pois", "pois"),
#'                       comp.param = list(list("lambda" = 7),
#'                                         list("lambda" = 1)))
#' data1 <- get_mixture_data(mixt1)
#' data2 <- get_mixture_data(mixt2)
#' data3 <- get_mixture_data(mixt3)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' admixMod3 <- admix_model(knownComp_dist = mixt3$comp.dist[[2]],
#'                          knownComp_param = mixt3$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1,data2),
#'                    admixMod = list(admixMod1,admixMod2), est_method = 'IBM')
#' est2 <- admix_estim(samples = list(data3), admixMod = list(admixMod3), est_method = 'PS')
#' ## Determine the decontaminated version of the unknown density by inversion:
#' x <- decontaminated_density(sample1 = data1, admixMod = admixMod1,
#'                             estim.p = get_mixing_weights(est)[1])
#' y <- decontaminated_density(sample1 = data2, admixMod = admixMod2,
#'                             estim.p = get_mixing_weights(est)[2])
#' z <- decontaminated_density(sample1 = data3, admixMod = admixMod3,
#'                             estim.p = get_mixing_weights(est2))
#' plot(x, offset = -0.2, bar_width = 0.2, col = "steelblue")
#' plot(y, add_plot = TRUE, offset = 0, bar_width = 0.2, col = "red")
#' plot(z, add_plot = TRUE, offset = 0.2, bar_width = 0.2, col = "orange")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export
#'

decontaminated_density <- function(sample1, admixMod, estim.p)
{
  ##--- Safety checks ---##
  if (!inherits(x = admixMod, what = "admix_model"))
    stop("Argument 'admixMod' must be of class 'admix_model'. See ?admix_model.")
  if (!is.numeric(sample1) || !is.numeric(estim.p))
    stop("'sample1' and 'estim.p' must be numeric.")
  if (estim.p <= 0 || estim.p > 1)
    stop("'estim.p' must be in (0,1].")

  ##--- Extract known component density ---##
  knownComp.dens <- paste0("d", admixMod$comp.dist$known)
  if (knownComp.dens == "dmultinom") knownComp.dens <- "approxfun"
  comp_dens_fun <- get(knownComp.dens, mode = "function")
  if (knownComp.dens == "approxfun") { # Discrete case
    expr <- paste0("approxfun(x = 1:", length(admixMod$comp.param$known$prob),
                   ", y = c(", paste(admixMod$comp.param$known$prob, collapse = ","), "))")
    g1.fun <- eval(parse(text = expr))
    g1 <- function(x) g1.fun(x)
  } else { # Continuous distribution
    expr <- paste0(knownComp.dens, "(x, ",
                   paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known,
                         collapse = ", "), ")")
    g1 <- function(x) eval(parse(text = expr))
  }

  ##--- Detect support type ---##
  support <- detect_support_type(sample1)

  ##--- Empirical density estimation ---##
  if (support == "Continuous") {
    l1_dens <- stats::density(sample1, n = 2048)
    l1_emp <- stats::approxfun(l1_dens$x, l1_dens$y, yleft = 0, yright = 0)
  } else {
    tbl <- table(sample1)
    l1_emp <- stats::approxfun(
      x = as.numeric(names(tbl)),
      y = as.numeric(tbl) / sum(tbl),
      method = "constant", yleft = 0, yright = 0
    )
  }
  ##--- Compute decontaminated density ---##
  f1_decontamin <- function(x) {
    val <- (1 / estim.p) * (l1_emp(x) - (1 - estim.p) * g1(x))
    val[val < 0] <- 0  # avoid negatives due to estimation noise
    return(val)
  }
  ##--- Construct output object ---##
  obj <- list(
    data = sample1,
    support = support,
    decontaminated_density_fun = f1_decontamin,
    call = match.call()
  )
  class(obj) <- "decontaminated_density"
  return(obj)
}

#' Print method for object of class \code{decontaminated_density}
#'
#' Print some overview of the decontaminated density function.
#'
#' @param x An object of class \code{decontaminated_density} (see ?decontaminated_density).
#' @param ... Arguments to be passed to generic method \code{plot}, such as graphical parameters (see ?par).
#'
#' @return The function related to the estimated decontaminated density.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.decontaminated_density <- function(x, ...)
{
  cat("Call:")
  print(x$call)
  cat("\n")
  cat("Statistics about the estimated decontaminated density function:\n")
  print(summary(x$decontaminated_density_fun(x$data)))
  cat("\n")
}


#' Summary method for object of class \code{decontaminated_density}
#'
#' Summarizes information about the estimated decontaminated density function.
#'
#' @param object An object of class \code{decontaminated_density} (see ?decontaminated_density).
#' @param ... Arguments to be passed to generic method \code{summary}.
#'
#' @return Classical statistical indicators about the decontaminated density.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

summary.decontaminated_density <- function(object, ...)
{
  cat("Call:")
  print(object$call)
  cat("\n")
  cat("Type of support: ", object$support, sep = "")
  cat("\n")
  if (object$support == "Continuous") {
    cat("Statistical indicators about the support:\n")
    print(summary(object$data))
  } else {
    cat("Count table:")
    print(table(object$data))
  }
  cat("\n")
  cat("Statistics about the estimated decontaminated density function:\n")
  print(summary(object$decontaminated_density_fun(object$data)))
  cat("\n")
}


#' Plot method for object of class \code{decontaminated_density}
#'
#' Plot the decontaminated density of the unknown component from some admixture model, after inversion of the admixture
#' cumulative distribution functions.
#'
#' @param x An object of class \code{decontaminated_density} (see ?decontaminated_density).
#' @param x_val Values at which to evaluate the decontaminated density.
#' @param add_plot Boolean, TRUE when a new plot is added to the existing one.
#' @param offset Numeric. Position of the bars relative to the labels on the x-axis.
#' @param bar_width Width of bars to be plotted.
#' @param ... Arguments to be passed to generic method \code{plot}, such as graphical parameters (see ?par).
#'
#' @details The decontaminated density is obtained by inverting the admixture density, given by l = p*f + (1-p)*g, to isolate the
#'          unknown component f after having estimated p and l.
#'
#' @return The plot of the decontaminated density if one sample is provided, or the comparison of decontaminated
#'         densities plotted on the same graph in the case of multiple samples.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot.decontaminated_density <- function(x, x_val = NULL, add_plot = FALSE, offset = 0, bar_width = 0.3, ...)
{
  support <- x$support
  if (!is.null(x_val)) {
    decontamin_dens_values <- x$decontaminated_density_fun(x_val)
  } else {
    if (support == "Discrete") {
      x_val <- seq(from = min(x$data), to = max(x$data))
    } else {
      x_val <- seq(from = min(x$data), to = max(x$data), length.out = 100)
    }
    decontamin_dens_values <- x$decontaminated_density_fun(x_val)
  }
  decontamin_dens_values[is.na(decontamin_dens_values)] <- 0
  decontamin_dens_values <- pmax(decontamin_dens_values, 0)

  if (!add_plot) {
    if (support == "Discrete") {
      ## initialise an empty plot
      plot(range(x_val), range(0, decontamin_dens_values * 1.1),
           type="n", xaxt="n", xlab="x", ylab="Density", ...)
      graphics::axis(1, at=x_val, labels=as.character(x_val))
    } else {
      plot(x = x_val, y = decontamin_dens_values, type="l", ...)
    }
  }

  if (support == "Discrete") {
    for (i in seq_along(x_val)) {
      graphics::rect(xleft = x_val[i] - bar_width/2 + offset,
                     xright = x_val[i] + bar_width/2 + offset,
                     ybottom = 0, ytop = decontamin_dens_values[i], ...)
    }
  } else {
    graphics::lines(x = x_val, y = decontamin_dens_values, ...)
  }
}


#' Estimates the decontaminated CDF of the unknown component in an admixture
#'
#' Estimates the decontaminated cumulative distribution function (CDF) of the unknown component in an admixture model, using inversion of the admixture
#' CDF. Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with g a known CDF and p and f unknown quantities.
#'
#' @param sample1 Observations of the sample under study.
#' @param admixMod An object of class \code{admix_model}, containing useful information about distributions and parameters.
#' @param estim.p The estimated weight of the unknown component distribution, related to the proportion of the unknown component
#'                  in the admixture model studied.
#'
#' @details The decontaminated CDF is obtained by inverting the admixture CDF, given by L = p*F + (1-p)*G, to isolate the
#'          unknown component F after having estimated p. This means that F = (1/hat(p)) * (hat(L)-(1-p)*G).
#'
#' @return The decontaminated CDF F of the admixture model, of class 'stepfun' (step function).
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 400, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 300, weight = 0.6,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 5, "sd" = 2)))
#' data1 <- get_mixture_data(mixt1)
#' data2 <- get_mixture_data(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1,data2),
#'                    admixMod = list(admixMod1,admixMod2), est_method = 'PS')
#' prop <- get_mixing_weights(est)
#' ## Determine the decontaminated version of the unknown CDF by inversion:
#' F1 <- decontaminated_cdf(sample1 = data1, admixMod = admixMod1, estim.p = prop[1])
#' F2 <- decontaminated_cdf(sample1 = data2, admixMod = admixMod2, estim.p = prop[2])
#' abs <- seq(from=-1, to=4, length.out=100)
#' plot(x=abs, y=F1(abs), xlim=c(2,4), ylim=c(0,1), type="l")
#' par(new = TRUE)
#' plot(x=abs, y=F2(abs), xlim=c(2,4), ylim=c(0,1), type="l", col="red")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @keywords internal
#' @export

decontaminated_cdf <- function(sample1, admixMod, estim.p)
{
  if (!inherits(x = admixMod, what = "admix_model"))
    stop("Argument 'admixMod' is not correctly specified. See ?admix_model.")

  ## Extract the information on component distributions:
  knownComp.CDF <- paste0("p", admixMod$comp.dist$known)
  if (any(knownComp.CDF == "pmultinom")) { knownComp.CDF[which(knownComp.CDF == "pmultinom")] <- "stepfun" }
  comp_cdf <- sapply(X = knownComp.CDF, FUN = get, mode = "function")
  for (i in 1:length(comp_cdf)) assign(x = names(comp_cdf)[i], value = comp_cdf[[i]])

  ## Create the expression involved in future assessments of the densities:
  if (knownComp.CDF == "stepfun") {
    expr <- paste(names(comp_cdf),"(x = 1:", length(admixMod$comp.param$known$prob), paste(", y = ", paste("cumsum(c(0,",
              paste(admixMod$comp.param$known$prob, collapse = ","), "))", sep = ""), ")", sep = ""), sep = "")
  } else {
    expr <- paste(names(comp_cdf),"(z,", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known,
                                                sep = "", collapse = ","), ")", sep="")
  }

  ## Differentiates the case of the multinomial distribution from other cases:
  if (any(knownComp.CDF == "stepfun")) {
    G1.fun <- eval(parse(text = expr))
    G1 <- function(z) G1.fun(z)
  } else {
    G1 <- function(z) { eval(parse(text = expr)) }
  }
  ## Empirical cumulative distribution function from the observed sample:
  L1 <- stats::ecdf(sample1)						# hat(L1)

  ##------- Defines the support --------##
  support <- detect_support_type(sample1)
  ## Range of observations:
  if (support == "Continuous") {
    x_values <- seq(from = min(sample1), to = max(sample1), length.out = 3000)
  } else {
    x_values <- unique(sort(c(unique(sample1))))
  }

  F1_decontamin <- sapply(X = x_values, FUN = function(x) (1/estim.p) * (L1(x) - (1-estim.p) * G1(x)))
  F1_decontamin[which(is.na(F1_decontamin))] <- 0
  F1_monotonous <- Iso::pava(y = c(0, F1_decontamin), decreasing = FALSE)
  F1_monotonous[which(F1_monotonous <= 0)] <- 0
  F1_monotonous[which(F1_monotonous >= 1)] <- 1
  F1_fun <- stats::stepfun(x = x_values, y = F1_monotonous)

  return(F1_fun)
}
