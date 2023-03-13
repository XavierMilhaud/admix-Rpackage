test_that("good estimation of variance-covariance matrix of the gaussian process", {
  testthat::skip_on_cran()
  library(foreach)
  library(doParallel)
  n1 <- 5000
  n2 <- 4800
  z <- 2
  list.comp <- list(f1 = "gamma", g1 = "exp",
                    f2 = "gamma", g2 = "exp")
  list.param <- list(f1 = list(shape = 2, scale = 2), g1 = list(rate = 1/3),
                     f2 = list(shape = 2, scale = 2), g2 = list(rate = 1/5))
  resultats <-
    foreach::foreach (k = 1:10, .inorder = TRUE, .errorhandling = "remove") %do% {
      X.sim <- rsimmix(n=n1, unknownComp_weight=0.4, comp.dist = list(list.comp$f1,list.comp$g1),
                       comp.param = list(list.param$f1, list.param$g1))$mixt.data
      Y.sim <- rsimmix(n=n2, unknownComp_weight=0.6, comp.dist = list(list.comp$f2,list.comp$g2),
                       comp.param = list(list.param$f2, list.param$g2))$mixt.data
      estim <- IBM_estimProp(sample1 = X.sim, sample2 = Y.sim, known.prop = c(0.4,0.6),
                             comp.dist = list.comp, comp.param = list.param,
                             with.correction = FALSE, n.integ = 1000)
      D_z <- IBM_gap(z = z, par = estim[["prop.estim"]], fixed.p1 = NULL, sample1 = X.sim,
                     sample2 = Y.sim, comp.dist = list.comp, comp.param = list.param)
      D_z.theo <- IBM_theoretical_gap(z = z, par = estim[["theo.prop.estim"]], known.p = c(0.4,0.6),
                                      comp.dist = list.comp, comp.param = list.param)
      list(estim[["prop.estim"]], estim[["theo.prop.estim"]], D_z, D_z.theo)
    }
  vect.p1 <- sapply(resultats, "[[", 1)[1, ] ; vect.p1.theo <- sapply(resultats, "[[", 2)[1, ]
  vect.p2 <- sapply(resultats, "[[", 1)[2, ] ; vect.p2.theo <- sapply(resultats, "[[", 2)[2, ]
  vect.D_z <- sapply(resultats, "[[", 3) ; vect.D_z.theo <- sapply(resultats, "[[", 4)
  empirical_var_gaussVectcomponents <- round(c(var(sqrt(min(n1,n2)) * (vect.p1 - vect.p1.theo)),
                                               var(sqrt(min(n1,n2)) * (vect.p2 - vect.p2.theo)),
                                               var(sqrt(min(n1,n2)) * (vect.D_z - vect.D_z.theo))))
  ## Then validate theoretical results by comparing them to empirical experience:
  X.sim <- rsimmix(n=n1, unknownComp_weight=0.4, comp.dist = list(list.comp$f1,list.comp$g1),
                   comp.param = list(list.param$f1, list.param$g1))$mixt.data
  Y.sim <- rsimmix(n=n2, unknownComp_weight=0.6, comp.dist = list(list.comp$f2,list.comp$g2),
                   comp.param = list(list.param$f2, list.param$g2))$mixt.data
  estim <- IBM_estimProp(sample1 = X.sim, sample2 = Y.sim, known.prop = c(0.4,0.6),
                         comp.dist = list.comp, comp.param = list.param,
                         with.correction = FALSE, n.integ = 1000)
#  varCov_estim <- IBM_estimVarCov_gaussVect(x = z, y = z, estim.obj = estim,
#                                            fixed.p1 = estim[["p.X.fixed"]], known.p = c(0.4,0.6),
#                                            sample1 = X.sim, sample2 = Y.sim, min_size = NULL,
#                                            comp.dist = list.comp, comp.param = list.param)
  expect_equal(0, 0)
})
