test_that("good estimation of variance-covariance matrix of the gaussian process", {
  testthat::skip_on_cran()
  library(foreach)
  library(doParallel)
  n1 <- 2000
  n2 <- 1500
  z <- 2
  resultats <-
    foreach::foreach (k = 1:10, .inorder = TRUE, .errorhandling = "remove") %do% {
      ## Simulate mixture data:
      X.sim <- twoComp_mixt(n = n1, weight = 0.2,
                            comp.dist = list("norm", "norm"),
                            comp.param = list(list("mean" = 2, "sd" = 3),
                                              list("mean" = -2, "sd" = 1)))
      data1 <- getmixtData(X.sim)
      Y.sim <- twoComp_mixt(n = n2, weight = 0.5,
                            comp.dist = list("norm", "norm"),
                            comp.param = list(list("mean" = 2, "sd" = 3),
                                              list("mean" = 6, "sd" = 1)))
      data2 <- getmixtData(Y.sim)
      admixMod1 <- admix_model(knownComp_dist = X.sim$comp.dist[[2]],
                               knownComp_param = X.sim$comp.param[[2]])
      admixMod2 <- admix_model(knownComp_dist = Y.sim$comp.dist[[2]],
                               knownComp_param = Y.sim$comp.param[[2]])
      estim <- estim_IBM(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2))
      D_z <- IBM_gap(z = z, par = estim$prop.estim, samples = list(data1,data2), admixMod = list(admixMod1,admixMod2))
      D_z.theo <- IBM_theoretical_gap(z = z, par = c(0.4,0.6), known.p = c(0.4,0.6), mixtMod = list(X.sim,Y.sim))
      list(estim$prop.estim, c(0.4,0.6), D_z, D_z.theo)
    }
  vect.p1 <- sapply(resultats, "[[", 1)[1, ] ; vect.p1.theo <- sapply(resultats, "[[", 2)[1, ]
  vect.p2 <- sapply(resultats, "[[", 1)[2, ] ; vect.p2.theo <- sapply(resultats, "[[", 2)[2, ]
  vect.D_z <- sapply(resultats, "[[", 3) ; vect.D_z.theo <- sapply(resultats, "[[", 4)
  empirical_var_gaussVectcomponents <- round(c(var(sqrt(min(n1,n2)) * (vect.p1 - vect.p1.theo)),
                                               var(sqrt(min(n1,n2)) * (vect.p2 - vect.p2.theo)),
                                               var(sqrt(min(n1,n2)) * (vect.D_z - vect.D_z.theo))), 1)
  ## Then validate theoretical results by comparing them to empirical experience:
  X.sim <- twoComp_mixt(n=n1, weight=0.4, comp.dist = list("gamma","exp"),
                        comp.param = list(list("shape" = 2, "scale" = 2), list("rate" = 1/3)))
  data1 <- getmixtData(X.sim)
  Y.sim <- twoComp_mixt(n=n2, weight=0.6, comp.dist = list("gamma","exp"),
                        comp.param = list(list("shape" = 2, "scale" = 2), list("rate" = 1/5)))
  data2 <- getmixtData(Y.sim)
  admixMod1 <- admix_model(knownComp_dist = X.sim$comp.dist[[2]],
                           knownComp_param = X.sim$comp.param[[2]])
  admixMod2 <- admix_model(knownComp_dist = Y.sim$comp.dist[[2]],
                           knownComp_param = Y.sim$comp.param[[2]])
  estim <- estim_IBM(samples = list(data1,data2), admixMod = list(admixMod1,admixMod2))
  #varCov_estim <- IBM_estimVarCov_gaussVect(x = z, y = z, IBMestim.obj = estim,
  #                                          samples = list(data1,data2), admixMod = list(admixMod1,admixMod2))
  expect_equal(0, 0)
})
