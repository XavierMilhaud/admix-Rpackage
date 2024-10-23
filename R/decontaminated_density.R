#' Estimates the decontaminated density of the unknown component in an admixture
#'
#' Estimate the decontaminated density of the unknown component in the admixture model under study, after inversion of the admixture
#' cumulative distribution function. Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with g a known CDF and p and f unknown quantities.
#'
#' @param sample1 Sample under study.
#' @param estim.p The estimated weight of the unknown component distribution, related to the proportion of the unknown component
#'                in the admixture model studied.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#'
#' @details The decontaminated density is obtained by inverting the admixture density, given by l = p*f + (1-p)*g, to isolate the
#'          unknown component f after having estimated p.
#'
#' @return An object of class 'decontaminated_density', containing 2 attributes: 1) the decontaminated density function;
#'         2) the type of support for the underlying distribution (either discrete or continuous, useful for plots).
#'
#' @examples
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 400, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = -2, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' data1 <- getmixtData(mixt1)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1), admixMod = list(admixMod1),
#'                    est.method = 'PS')
#' ## Determine the decontaminated version of the unknown density by inversion:
#' decontaminated_density(sample1 = data1, estim.p = est$estimated_mixing_weights[1],
#'                        admixMod = admixMod1)
#'
#' ####### Discrete support:
#' mixt1 <- twoComp_mixt(n = 5000, weight = 0.6,
#'                       comp.dist = list("pois", "pois"),
#'                       comp.param = list(list("lambda" = 3),
#'                                         list("lambda" = 2)))
#' mixt2 <- twoComp_mixt(n = 4000, weight = 0.8,
#'                       comp.dist = list("pois", "pois"),
#'                       comp.param = list(list("lambda" = 3),
#'                                         list("lambda" = 4)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1, data2),
#'                    admixMod = list(admixMod1, admixMod2), est.method = 'IBM')
#' ## Determine the decontaminated version of the unknown density by inversion:
#' decontaminated_density(sample1 = data1, estim.p = est$estimated_mixing_weights[1],
#'                        admixMod = admixMod1)
#'
#' ####### Finite discrete support:
#' mixt1 <- twoComp_mixt(n = 12000, weight = 0.6,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size" = 1, "prob" = c(0.3,0.4,0.3)),
#'                                         list("size" = 1, "prob" = c(0.6,0.3,0.1))))
#' mixt2 <- twoComp_mixt(n = 10000, weight = 0.8,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size" = 1, "prob" = c(0.3,0.4,0.3)),
#'                                         list("size" = 1, "prob" = c(0.2,0.6,0.2))))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1, data2),
#'                    admixMod = list(admixMod1, admixMod2), est.method = 'IBM')
#' ## Determine the decontaminated version of the unknown density by inversion:
#' decontaminated_density(sample1 = data1, estim.p = est$estimated_mixing_weights[1],
#'                        admixMod = admixMod1)
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

decontaminated_density <- function(sample1, estim.p, admixMod)
{
  ## Extract the information on component distributions:
  knownComp.dens <- paste0("d", admixMod$comp.dist$known)
  if (any(knownComp.dens == "dmultinom")) { knownComp.dens[which(knownComp.dens == "dmultinom")] <- "approxfun" }
  comp_dens <- sapply(X = knownComp.dens, FUN = get, mode = "function")
  for (i in 1:length(comp_dens)) assign(x = names(comp_dens)[i], value = comp_dens[[i]])

  ## Create the expression involved in future assessments of the densities:
  if (knownComp.dens == "approxfun") {
    expr <- paste(names(comp_dens),"(x = 1:", length(admixMod$comp.param$known$prob), paste(", y = c(",
                                paste(paste(admixMod$comp.param$known$prob, collapse = ","), ")", sep = ""), ")", sep = ""), sep = "")
  } else {
    expr <- paste(names(comp_dens),"(x,", paste(names(admixMod$comp.param$known), "=", admixMod$comp.param$known,
                                                sep = "", collapse = ","), ")", sep="")
  }

  ## Evaluate density of known components, differentiates the case of the multinomial distribution from others:
  if (any(knownComp.dens == "approxfun")) {
    g1.fun <- eval(parse(text = expr))
    g1 <- function(x) g1.fun(x)
  } else {
    g1 <- function(x) { eval(parse(text = expr)) }
  }

  ##------- Defines the support --------##
  support <- detect_support_type(sample1)

  ## Retrieves the empirical cumulative distribution functions:
  if (support == "Continuous") {
    l1_dens <- stats::density(sample1)
    l1_emp <- stats::approxfun(x = l1_dens$x, y = l1_dens$y)
  } else {
    l1_emp <- stats::approxfun(x = as.numeric(names(table(sample1))),
                               y = table(sample1) / sum(table(sample1)), method = "constant")
  }
  f1_decontamin <- function(x) (1/estim.p) * (l1_emp(x) - (1-estim.p) * g1(x))

  obj <- list(decontaminated_density_fun = f1_decontamin,
              support = support)
  class(obj) <- "decontaminated_density"
  obj$call <- match.call()

  return(obj)
}


#' Print method for object of class 'decontaminated_density'
#'
#' Print some overview of the decontaminated density function.
#'
#' @param x An object of class 'decontaminated_density' (see ?decontaminated_density).
#' @param ... Arguments to be passed to generic method 'plot', such as graphical parameters (see par).
#'
#' @return More important information about the decontaminated density.
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

print.decontaminated_density <- function(x, ...)
{
  cat("Call:")
  print(x$call)
  cat("\n")
  print(x$decontaminated_density_fun)
  cat("\n")
  cat("Type of support: ", x$support, sep = "")
}



#' Plot method for class 'decontaminated_density'
#'
#' Plot the decontaminated density of the unknown component from some admixture model, after inversion of the admixture
#' cumulative distribution functions.
#'
#' @param x An object of class 'decontaminated_density' (see ?decontaminated_density).
#' @param x_val (numeric) A vector of points at which to evaluate the probability mass/density function.
#' @param add_plot (default to FALSE) A boolean specifying if one plots the decontaminated density over an existing plot. Used for visual
#'                 comparison purpose.
#' @param ... Arguments to be passed to generic method 'plot', such as graphical parameters (see par).
#'
#' @details The decontaminated density is obtained by inverting the admixture density, given by l = p*f + (1-p)*g, to isolate the
#'          unknown component f after having estimated p and l.
#'
#' @return The plot of the decontaminated density.
#'
#' @examples
#' ####### Continuous support:
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 400, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 350, weight = 0.6,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 5, "sd" = 2)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'                    est.method = 'PS')
#' prop <- getmixingWeight(est)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' res1 <- decontaminated_density(sample1 = data1, estim.p = prop[1], admixMod = admixMod1)
#' res2 <- decontaminated_density(sample1 = data2, estim.p = prop[2], admixMod = admixMod2)
#' ## Use appropriate sequence of x values:
#' plot(x = res1, x_val = seq(from = 0, to = 6, length.out = 100), add_plot = FALSE)
#' plot(x = res2, col = "red", x_val = seq(from = 0, to = 6, length.out = 100), add_plot = TRUE)
#'
#' ####### Countable discrete support:
#' mixt1 <- twoComp_mixt(n = 4000, weight = 0.7,
#'                       comp.dist = list("pois", "pois"),
#'                       comp.param = list(list("lambda" = 3),
#'                                         list("lambda" = 2)))
#' mixt2 <- twoComp_mixt(n = 3500, weight = 0.85,
#'                       comp.dist = list("pois", "pois"),
#'                       comp.param = list(list("lambda" = 3),
#'                                         list("lambda" = 4)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'                    est.method = "IBM")
#' prop <- getmixingWeight(est)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' res1 <- decontaminated_density(sample1 = data1, estim.p = prop[1],
#'                                admixMod = admixMod1)
#' res2 <- decontaminated_density(sample1 = data2, estim.p = prop[2],
#'                                admixMod = admixMod2)
#' ## Use appropriate sequence of x values:
#' plot(x = res1, x_val = seq(from = 0, to = 15, by = 1), add_plot = FALSE)
#' plot(x = res2, x_val = seq(from = 0, to = 15, by = 1), add_plot = TRUE, col = "red")
#'
#' ####### Finite discrete support:
#' mixt1 <- twoComp_mixt(n = 4000, weight = 0.7,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size" = 1, "prob" = c(0.3,0.4,0.3)),
#'                                         list("size" = 1, "prob" = c(0.6,0.3,0.1))))
#' mixt2 <- twoComp_mixt(n = 3500, weight = 0.85,
#'                       comp.dist = list("multinom", "multinom"),
#'                       comp.param = list(list("size" = 1, "prob" = c(0.3,0.4,0.3)),
#'                                         list("size" = 1, "prob" = c(0.2,0.6,0.2))))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'                    est.method = "IBM")
#' prop <- getmixingWeight(est)
#' ## Determine the decontaminated version of the unknown density by inversion:
#' res1 <- decontaminated_density(sample1 = data1, estim.p = prop[1],
#'                                admixMod = admixMod1)
#' res2 <- decontaminated_density(sample1 = data2, estim.p = prop[2],
#'                                admixMod = admixMod2)
#' ## Use appropriate sequence of x values:
#' plot(x = res1, x_val = seq(from=1, to=3, by=1), add_plot = FALSE)
#' plot(x = res2, x_val = seq(from=1, to=3, by=1), add_plot = TRUE, col = "red")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

plot.decontaminated_density <- function(x, x_val, add_plot = FALSE, ...)
{
  support <- x$support
  decontamin_dens_values <- x$decontaminated_density_fun(x_val)
  decontamin_dens_values[which(is.na(decontamin_dens_values))] <- 0
  decontamin_dens_values <- pmax(decontamin_dens_values, 0)
  y.range <- c(min(decontamin_dens_values), max(decontamin_dens_values)*1.1)

  if (!add_plot) {
    if (support == "Discrete") {
      graphics::barplot(height = decontamin_dens_values, names = as.character(x_val), space = 0.1, ...)
    } else {
      plot(x = x_val, y = decontamin_dens_values, ...)
    }
  } else {
    old_par_new <- graphics::par()$new
    graphics::par(new = TRUE)
    if (support == "Discrete") {
      graphics::barplot(height = decontamin_dens_values, names = as.character(x_val), add = add_plot,
                        space = 0.1, ...)
    } else {
      graphics::lines(x = x_val, y = decontamin_dens_values, ...)
    }
    on.exit(graphics::par(new = old_par_new))
  }
}


#' Estimates the decontaminated CDF of the unknown component in an admixture
#'
#' Estimates the decontaminated cumulative distribution function (CDF) of the unknown component in an admixture model, using inversion of the admixture
#' CDF. Recall that an admixture model follows the cumulative distribution function (CDF) L, where
#' L = p*F + (1-p)*G, with g a known CDF and p and f unknown quantities.
#'
#' @param sample1 Observations of the sample under study.
#' @param estim.p The estimated weight of the unknown component distribution, related to the proportion of the unknown component
#'                  in the admixture model studied.
#' @param admixMod An object of class 'admix_model', containing useful information about distributions and parameters.
#'
#' @details The decontaminated CDF is obtained by inverting the admixture CDF, given by L = p*F + (1-p)*G, to isolate the
#'          unknown component F after having estimated p. This means that F = (1/hat(p)) * (hat(L)-(1-p)*G).
#'
#' @return The decontaminated CDF F of the admixture model, of class 'stepfun' (step function).
#'
#' @examples
#' ####### Continuous support:
#' ## Simulate mixture data:
#' mixt1 <- twoComp_mixt(n = 400, weight = 0.4,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 0, "sd" = 1)))
#' mixt2 <- twoComp_mixt(n = 300, weight = 0.6,
#'                       comp.dist = list("norm", "norm"),
#'                       comp.param = list(list("mean" = 3, "sd" = 0.5),
#'                                         list("mean" = 5, "sd" = 2)))
#' data1 <- getmixtData(mixt1)
#' data2 <- getmixtData(mixt2)
#' ## Define the admixture models:
#' admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]],
#'                          knownComp_param = mixt1$comp.param[[2]])
#' admixMod2 <- admix_model(knownComp_dist = mixt2$comp.dist[[2]],
#'                          knownComp_param = mixt2$comp.param[[2]])
#' ## Estimation:
#' est <- admix_estim(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2),
#'                    est.method = 'PS')
#' prop <- getmixingWeight(est)
#' ## Determine the decontaminated version of the unknown CDF by inversion:
#' F1 <- decontaminated_cdf(sample1 = data1, estim.p = prop[1], admixMod = admixMod1)
#' F2 <- decontaminated_cdf(sample1 = data2, estim.p = prop[2], admixMod = admixMod2)
#' abs <- seq(from=-1, to=4, length.out=100)
#' plot(x=abs, y=F1(abs), xlim=c(-1,4), ylim=c(0,1), type="l")
#' par(new = TRUE)
#' plot(x=abs, y=F2(abs), xlim=c(-1,4), ylim=c(0,1), type="l", col="red")
#'
#' @author Xavier Milhaud <xavier.milhaud.research@gmail.com>
#' @export

decontaminated_cdf <- function(sample1, estim.p, admixMod)
{
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
  F1_decontamin <- pmax(F1_decontamin, 0)
  F1_decontamin <- pmin(F1_decontamin, 1)

  F1_fun <- stats::stepfun(x = x_values, y = c(0,F1_decontamin))
  #F1_fun <- Iso::pava(y = c(0,F1_decontamin), decreasing = FALSE, stepfun = TRUE)
  #plot.stepfun(x = F1_fun, xlim = c(min(x_values[knots(F1_fun)-1]), max(x_values[knots(F1_fun)-1])))
  return(F1_fun)
}
