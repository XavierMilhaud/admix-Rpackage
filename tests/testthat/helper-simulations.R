create_test_samples <- function() {
  set.seed(123)
  mixt1 <- twoComp_mixt(n = 120, weight = 0.4, comp.dist = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 100, weight = 0.6, comp.dist = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 1, sd = 1)))
  list(
    samples = list(get_mixture_data(mixt1), get_mixture_data(mixt2)),
    admix = list(admix_model(knownComp_dist = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]]),
                 admix_model(knownComp_dist = mixt2$comp.dist[[2]], knownComp_param = mixt2$comp.param[[2]]))
  )
}
