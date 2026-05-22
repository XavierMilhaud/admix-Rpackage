
# ------------------------------------------------------------------------------
# Helpers : données simulées réutilisées dans tous les tests
# ------------------------------------------------------------------------------

make_PS_data <- function(seed = 42) {
  set.seed(seed)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.4,
                        comp.dist  = list("gamma", "exp"),
                        comp.param = list(list(shape = 2, scale = 0.5),
                                          list(rate = 0.25)))
  list(
    data1     = get_mixture_data(mixt1),
    admixMod1 = admix_model(knownComp_dist  = mixt1$comp.dist[[2]],
                            knownComp_param = mixt1$comp.param[[2]])
  )
}

# ==============================================================================
# 1. estim_PS — Validation des entrées (erreurs attendues)
# ==============================================================================

test_that("estim_PS : erreur si admixMod n'est pas de la bonne classe", {
  d <- make_PS_data()
  expect_error(
    estim_PS(samples  = d$data1,
             admixMod = list("not_a_model"),
             method   = "fixed"),
    regexp = "admixMod"
  )
})


test_that("estim_PS : erreur si samples n'est pas un vecteur", {
  d <- make_PS_data()
  expect_error(
    estim_PS(samples  = as.data.frame(d$data1),
             admixMod = d$admixMod1,
             method   = "fixed"),
    regexp = "numerical vector"
  )
})


test_that("estim_PS : erreur si method est NULL", {
  d <- make_PS_data()
  expect_error(
    estim_PS(samples  = d$data1,
             admixMod = d$admixMod1,
             method   = NULL),
    regexp = "method"
  )
})


# ==============================================================================
# 2. estim_PS — Structure du retour, méthode 'fixed'
# ==============================================================================

test_that("estim_PS (fixed) : retourne un objet de classe 'estim_PS' et 'admix_estim'", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_s3_class(res, "estim_PS")
  expect_s3_class(res, "admix_estim")
})


test_that("estim_PS (fixed) : tous les champs attendus sont présents", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expected_fields <- c("n_populations", "population_sizes", "admixture_models",
                       "estimation_method", "estimated_mixing_weights",
                       "Fs.hat", "dist.out", "c.n", "alp.Lwr", "n",
                       "method", "call")
  expect_true(all(expected_fields %in% names(res)))
})


test_that("estim_PS (fixed) : n_populations vaut 1", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_equal(res$n_populations, 1L)
})


test_that("estim_PS (fixed) : population_sizes correspond à la taille de l'échantillon original", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  # population_sizes reflète n avant transformation
  expect_equal(res$n, length(d$data1))
})


test_that("estim_PS (fixed) : estimated_mixing_weights est un scalaire numérique dans [0, 1]", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_length(res$estimated_mixing_weights, 1)
  expect_true(is.numeric(res$estimated_mixing_weights))
  expect_gte(res$estimated_mixing_weights, 0)
  expect_lte(res$estimated_mixing_weights, 1)
})


test_that("estim_PS (fixed) : c.n est un scalaire numérique positif", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_length(res$c.n, 1)
  expect_true(is.numeric(res$c.n))
  expect_gt(res$c.n, 0)
})


test_that("estim_PS (fixed) : alp.Lwr est NULL (réservé à la méthode lwr.bnd)", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_null(res$alp.Lwr)
})


test_that("estim_PS (fixed) : cv.out est NULL", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_null(res$cv.out)
})


test_that("estim_PS (fixed) : estimation_method mentionne Patra and Sen", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_match(res$estimation_method, "Patra and Sen", fixed = TRUE)
})


test_that("estim_PS (fixed) : $call est conservé", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_false(is.null(res$call))
})


test_that("estim_PS (fixed) : Fs.hat contient les éléments x et y", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_true(!is.null(res$Fs.hat$x))
  expect_true(!is.null(res$Fs.hat$y))
  expect_equal(length(res$Fs.hat$x), length(res$Fs.hat$y))
})


test_that("estim_PS (fixed) : dist.out est de classe 'PS_dist_fun'", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  expect_s3_class(res$dist.out, "PS_dist_fun")
})


# ==============================================================================
# 3. estim_PS — Méthode 'lwr.bnd'
# ==============================================================================

test_that("estim_PS (lwr.bnd) : retourne un objet de classe 'estim_PS'", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "lwr.bnd")
  expect_s3_class(res, "estim_PS")
})


test_that("estim_PS (lwr.bnd) : alp.Lwr est un scalaire numérique dans [0, 1]", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "lwr.bnd")
  expect_length(res$alp.Lwr, 1)
  expect_true(is.numeric(res$alp.Lwr))
  expect_gte(res$alp.Lwr, 0)
  expect_lte(res$alp.Lwr, 1)
})


test_that("estim_PS (lwr.bnd) : estimated_mixing_weights est NULL", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "lwr.bnd")
  expect_null(res$estimated_mixing_weights)
})


test_that("estim_PS (lwr.bnd) : c.n est NULL", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "lwr.bnd")
  expect_null(res$c.n)
})


# ==============================================================================
# 4. estim_PS — Méthode 'cv'
# ==============================================================================

test_that("estim_PS (cv) : retourne un objet de classe 'estim_PS'", {
  d   <- make_PS_data()
  res <- estim_PS(samples    = d$data1,
                  admixMod   = d$admixMod1,
                  method     = "cv",
                  folds      = 3,
                  reps       = 1,
                  cn.length  = 5,
                  gridsize   = 100)
  expect_s3_class(res, "estim_PS")
})


test_that("estim_PS (cv) : cv.out est non-NULL", {
  d   <- make_PS_data()
  res <- estim_PS(samples    = d$data1,
                  admixMod   = d$admixMod1,
                  method     = "cv",
                  folds      = 3,
                  reps       = 1,
                  cn.length  = 5,
                  gridsize   = 100)
  expect_false(is.null(res$cv.out))
})


test_that("estim_PS (cv) : estimated_mixing_weights est dans [0, 1]", {
  d   <- make_PS_data()
  res <- estim_PS(samples    = d$data1,
                  admixMod   = d$admixMod1,
                  method     = "cv",
                  folds      = 3,
                  reps       = 1,
                  cn.length  = 5,
                  gridsize   = 100)
  expect_gte(res$estimated_mixing_weights, 0)
  expect_lte(res$estimated_mixing_weights, 1)
})


test_that("estim_PS (cv) : c.n est le c_n retenu par validation croisée (scalaire positif)", {
  d   <- make_PS_data()
  res <- estim_PS(samples    = d$data1,
                  admixMod   = d$admixMod1,
                  method     = "cv",
                  folds      = 3,
                  reps       = 1,
                  cn.length  = 5,
                  gridsize   = 100)
  expect_length(res$c.n, 1)
  expect_gt(res$c.n, 0)
})


# ==============================================================================
# 5. print.estim_PS — Comportement de la sortie console
# ==============================================================================

test_that("print.estim_PS : retourne l'objet invisiblement", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  out <- withVisible(print.estim_PS(res))
  expect_false(out$visible)
  expect_identical(out$value, res)
})


test_that("print.estim_PS (fixed) : affiche la taille d'échantillon", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output <- capture.output(print.estim_PS(res))
  expect_true(any(grepl("Sample size", output)))
})


test_that("print.estim_PS (fixed) : affiche le poids estimé", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output <- capture.output(print.estim_PS(res))
  expect_true(any(grepl("Mixing weight", output)))
})


test_that("print.estim_PS (fixed) : affiche c_n", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output <- capture.output(print.estim_PS(res))
  expect_true(any(grepl("c_n", output)))
})


test_that("print.estim_PS (fixed) : mentionne 'fixed c_n' comme méthode de sélection", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output <- capture.output(print.estim_PS(res))
  expect_true(any(grepl("fixed c_n", output, ignore.case = TRUE)))
})


test_that("print.estim_PS (lwr.bnd) : affiche 'lower bound' et non le poids estimé", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "lwr.bnd")
  output <- capture.output(suppressMessages(print.estim_PS(res)))
  expect_true(any(grepl("lower bound|Lower bound", output, ignore.case = TRUE)))
  expect_false(any(grepl("Mixing weight", output)))
})


test_that("print.estim_PS (cv) : mentionne 'cross-validation' comme méthode de sélection", {
  d   <- make_PS_data()
  res <- estim_PS(samples    = d$data1,
                  admixMod   = d$admixMod1,
                  method     = "cv",
                  folds      = 3,
                  reps       = 1,
                  cn.length  = 5,
                  gridsize   = 100)
  output <- capture.output(print.estim_PS(res))
  expect_true(any(grepl("cross-validation", output, ignore.case = TRUE)))
})


# ==============================================================================
# 6. summary.estim_PS — Comportement de la sortie console
# ==============================================================================

test_that("summary.estim_PS : retourne l'objet invisiblement", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  out <- withVisible(summary.estim_PS(res))
  expect_false(out$visible)
  expect_identical(out$value, res)
})


test_that("summary.estim_PS : affiche la taille d'échantillon", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output <- capture.output(summary.estim_PS(res))
  expect_true(any(grepl("Sample size", output)))
})


test_that("summary.estim_PS : affiche la composante connue", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output <- capture.output(summary.estim_PS(res))
  expect_true(any(grepl("Known component", output)))
})


test_that("summary.estim_PS : affiche le poids estimé (méthode fixed)", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output <- capture.output(summary.estim_PS(res))
  expect_true(any(grepl("Mixing weight", output)))
})


test_that("summary.estim_PS : affiche le call si show.call = TRUE (défaut)", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output_with    <- capture.output(summary.estim_PS(res, show.call = TRUE))
  output_without <- capture.output(summary.estim_PS(res, show.call = FALSE))
  expect_true(any(grepl("Call", output_with)))
  expect_false(any(grepl("Call", output_without)))
})


test_that("summary.estim_PS (lwr.bnd) : affiche 'lower bound' et non le poids estimé", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "lwr.bnd")
  output <- capture.output(summary.estim_PS(res))
  expect_true(any(grepl("lower bound", output, ignore.case = TRUE)))
  expect_false(any(grepl("Mixing weight", output)))
})


test_that("summary.estim_PS (fixed) : affiche la méthode de sélection 'fixed c_n'", {
  d   <- make_PS_data()
  res <- estim_PS(samples  = d$data1,
                  admixMod = d$admixMod1,
                  method   = "fixed")
  output <- capture.output(summary.estim_PS(res))
  expect_true(any(grepl("fixed c_n", output, ignore.case = TRUE)))
})


test_that("estim_PS handles NULL c.n", {
  mixt <- twoComp_mixt(n = 100, weight = 0.4, comp.dist = list("norm", "norm"),
                       comp.param = list(list(mean = 0, sd = 1), list(mean = 1, sd = 1)))
  x <- estim_PS(samples = get_mixture_data(mixt),
                admixMod = admix_model(knownComp_dist = mixt$comp.dist[[2]], knownComp_param = mixt$comp.param[[2]]),
                method = "fixed", c.n = NULL, gridsize = 50)

  expect_s3_class(x, "estim_PS")
})


test_that("estimCV_PS checks data type", {
  expect_error(estimCV_PS(data = list(1,2,3), admixMod = admix_model("norm", list(mean = 0, sd = 1))),
               "vector")
})


test_that("estimCV_PS checks cn.length", {
  mixt <- twoComp_mixt(n = 50, weight = 0.5, comp.dist = list("norm", "norm"),
                       comp.param = list(list(mean = 0, sd = 1),list(mean = 1, sd = 1)))
  expect_error(
    estimCV_PS(data = get_mixture_data(mixt),
               admixMod = admix_model(knownComp_dist = mixt$comp.dist[[2]], knownComp_param = mixt$comp.param[[2]]),
               cn.s = NULL, cn.length = NULL),
    "can not be null"
  )
})


test_that("estimCV_PS checks cn.s length", {
  mixt <- twoComp_mixt(n = 50, weight = 0.5, comp.dist = list("norm", "norm"),
                       comp.param = list(list(mean = 0, sd = 1),list(mean = 1, sd = 1)))
  expect_error(
    estimCV_PS(data = get_mixture_data(mixt),
               admixMod = admix_model(knownComp_dist = mixt$comp.dist[[2]],knownComp_param = mixt$comp.param[[2]]),
               cn.s = c(0.1, 0.2), cn.length = 3),
    "different"
  )
})


test_that("cv.score checks tr.data class", {
  expect_error(cv.score(tr.data = list(), test.data = rnorm(10), c.n = 0.1),
               "wrong type")
})


test_that("PS_unknownDensity_estim checks class", {
  expect_error(
    PS_unknownDensity_estim(list()),
    "only works"
  )
})


test_that("PS_unknownDensity_estim decreasing density", {
  mixt <- twoComp_mixt(n = 100, weight = 0.5, comp.dist = list("norm", "norm"),
                       comp.param = list(list(mean = 0, sd = 1),list(mean = 1, sd = 1)))
  res <- estim_PS(samples = get_mixture_data(mixt),
                  admixMod = admix_model(knownComp_dist = mixt$comp.dist[[2]],knownComp_param = mixt$comp.param[[2]]),
                  method = "fixed", gridsize = 50)
  out <- PS_unknownDensity_estim(res, dec.density = TRUE)

  expect_type(out, "list")
})


test_that("PS_unknownDensity_estim increasing density", {
  mixt <- twoComp_mixt(n = 100, weight = 0.5, comp.dist = list("norm", "norm"),
                       comp.param = list(list(mean = 0, sd = 1),list(mean = 1, sd = 1)))
  res <- estim_PS(samples = get_mixture_data(mixt),
                  admixMod = admix_model(knownComp_dist = mixt$comp.dist[[2]],knownComp_param = mixt$comp.param[[2]]),
                  method = "fixed", gridsize = 50)
  out <- PS_unknownDensity_estim(res, dec.density = FALSE)

  expect_type(out, "list")
})
