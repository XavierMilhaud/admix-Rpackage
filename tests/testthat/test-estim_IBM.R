
# ------------------------------------------------------------------------------
# Helpers : données simulées réutilisées dans tous les tests
# ------------------------------------------------------------------------------

## Cas 1 — composantes connues DISTINCTES (cas général, 2D)
make_distinct_data <- function(seed = 42) {
  set.seed(seed)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 400, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 5, sd = 2)))
  list(data1 = get_mixture_data(mixt1),
       data2 = get_mixture_data(mixt2),
       admixMod1 = admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]]),
       admixMod2 = admix_model(knownComp_dist  = mixt2$comp.dist[[2]], knownComp_param = mixt2$comp.param[[2]]))
}

## Cas 2 — composantes connues IDENTIQUES (optimisation 1D)
make_equal_data <- function(seed = 42) {
  set.seed(seed)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 400, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  list(data1 = get_mixture_data(mixt1),
       data2 = get_mixture_data(mixt2),
       admixMod1 = admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]]),
       admixMod2 = admix_model(knownComp_dist  = mixt2$comp.dist[[2]], knownComp_param = mixt2$comp.param[[2]]))
}

# ==============================================================================
# 1. estim_IBM — Validation des entrées (erreurs attendues)
# ==============================================================================

test_that("estim_IBM : erreur si le nombre d'échantillons != 2", {
  d <- make_distinct_data()
  expect_error(
    estim_IBM(samples = list(d$data1), admixMod = list(d$admixMod1, d$admixMod2)),
    regexp = "Wrong number of samples")
  expect_error(
    estim_IBM(samples = list(d$data1, d$data2, d$data1), admixMod = list(d$admixMod1, d$admixMod2)),
    regexp = "Wrong number of samples")
})


test_that("estim_IBM : erreur si admixMod n'est pas de la bonne classe", {
  d <- make_distinct_data()
  expect_error(
    estim_IBM(samples  = list(d$data1, d$data2), admixMod = list("not_a_model", d$admixMod2)),
    regexp = "admixMod")
})


# ==============================================================================
# 2. estim_IBM — Structure de l'objet retourné (composantes connues distinctes)
# ==============================================================================

test_that("estim_IBM : retourne un objet de classe 'estim_IBM' et 'admix_estim'", {
  skip_on_cran()
  d   <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_s3_class(res, "estim_IBM")
  expect_s3_class(res, "admix_estim")
})

test_that("estim_IBM : tous les champs attendus sont présents", {
  skip_on_cran()
  d   <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)

  expected_fields <- c("n_populations", "population_sizes", "admixture_models",
                       "estimation_method", "estimated_mixing_weights",
                       "variance_est_p1", "variance_est_p2",
                       "equal.knownComp", "p.X.fixed", "integ.supp", "call")
  expect_true(all(expected_fields %in% names(res)))
})


test_that("estim_IBM : n_populations vaut 2", {
  skip_on_cran()
  d   <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_equal(res$n_populations, 2L)
})


test_that("estim_IBM : population_sizes correspondent aux tailles réelles", {
  skip_on_cran()
  d   <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_equal(res$population_sizes, c(length(d$data1), length(d$data2)))
})


test_that("estim_IBM : estimated_mixing_weights contient 2 valeurs (composantes distinctes)", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_length(res$estimated_mixing_weights, 2)
  expect_true(is.numeric(res$estimated_mixing_weights))
})


test_that("estim_IBM : equal.knownComp est FALSE pour des composantes distinctes", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_false(res$equal.knownComp)
})


test_that("estim_IBM : p.X.fixed est NULL pour des composantes distinctes", {
  skip_on_cran()
  d   <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_null(res$p.X.fixed)
})


test_that("estim_IBM : variances sont NA si compute_var = FALSE (défaut)", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_true(is.na(res$variance_est_p1))
  expect_true(is.na(res$variance_est_p2))
})


test_that("estim_IBM : integ.supp est un vecteur numérique de longueur n.integ (support continu)", {
  skip_on_cran()
  d <- make_distinct_data()
  n_integ <- 150
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = n_integ)
  expect_true(is.numeric(res$integ.supp))
  expect_equal(length(res$integ.supp), n_integ)
})


test_that("estim_IBM : estimation_method mentionne IBM", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_match(res$estimation_method, "IBM", fixed = TRUE)
})


test_that("estim_IBM : l'appel est conservé dans $call", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_true(!is.null(res$call))
})


# ==============================================================================
# 3. estim_IBM — Cas composantes connues IDENTIQUES (optimisation 1D)
# ==============================================================================

test_that("estim_IBM : equal.knownComp est TRUE pour des composantes identiques", {
  skip_on_cran()
  d <- make_equal_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_true(res$equal.knownComp)
})


test_that("estim_IBM : p.X.fixed vaut 0.2 (valeur arbitraire) pour composantes identiques", {
  skip_on_cran()
  d <- make_equal_data()
  res <- estim_IBM(samples = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_equal(res$p.X.fixed, 0.2)
})


test_that("estim_IBM : estimated_mixing_weights contient 1 valeur (composantes identiques)", {
  skip_on_cran()
  d <- make_equal_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  expect_length(res$estimated_mixing_weights, 1)
})


# ==============================================================================
# 4. print.estim_IBM — Comportement de la sortie console
# ==============================================================================

test_that("print.estim_IBM : retourne l'objet invisiblement", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  out <- withVisible(print.estim_IBM(res))
  expect_false(out$visible)
  expect_identical(out$value, res)
})


test_that("print.estim_IBM : affiche les tailles d'échantillons", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(print.estim_IBM(res))
  expect_true(any(grepl("Sample sizes", output)))
})


test_that("print.estim_IBM : affiche les poids estimés (composantes distinctes)", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(print.estim_IBM(res))
  expect_true(any(grepl("Mixing weight", output)))
  expect_gte(sum(grepl("Mixing weight", output)), 2)
})


test_that("print.estim_IBM : mentionne 'distinct known components' si G1 != G2", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(print.estim_IBM(res))
  expect_true(any(grepl("distinct known components", output, ignore.case = TRUE)))
})


test_that("print.estim_IBM : mentionne 'equal known components' si G1 == G2", {
  skip_on_cran()
  d <- make_equal_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(print.estim_IBM(res))
  expect_true(any(grepl("equal known components", output, ignore.case = TRUE)))
})


test_that("print.estim_IBM : n'affiche pas les variances si compute_var = FALSE", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(print.estim_IBM(res))
  expect_false(any(grepl("Variance", output)))
})


# ==============================================================================
# 5. summary.estim_IBM — Comportement de la sortie console
# ==============================================================================

test_that("summary.estim_IBM : retourne l'objet invisiblement", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  out <- withVisible(summary.estim_IBM(res))
  expect_false(out$visible)
  expect_identical(out$value, res)
})


test_that("summary.estim_IBM : affiche le nombre de populations", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(summary.estim_IBM(res))
  expect_true(any(grepl("Number of samples", output)))
})


test_that("summary.estim_IBM : affiche les tailles d'échantillons", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(summary.estim_IBM(res))
  expect_true(any(grepl("Sample sizes", output)))
})


test_that("summary.estim_IBM : affiche le support d'intégration", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(summary.estim_IBM(res))
  expect_true(any(grepl("Support", output, ignore.case = TRUE)))
})


test_that("summary.estim_IBM : affiche le call si show.call = TRUE (défaut)", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output_with <- capture.output(summary.estim_IBM(res, show.call = TRUE))
  output_without <- capture.output(summary.estim_IBM(res, show.call = FALSE))
  expect_true(any(grepl("Call", output_with)))
  expect_false(any(grepl("Call", output_without)))
})


test_that("summary.estim_IBM : affiche les composantes connues de chaque échantillon", {
  skip_on_cran()
  d <- make_distinct_data()
  res <- estim_IBM(samples  = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ  = 200)
  output <- capture.output(summary.estim_IBM(res))
  expect_true(any(grepl("Known component", output)))
})


test_that("estim_IBM calcule la variance quand compute_var = TRUE", {
  skip_on_cran()

  set.seed(1)
  data1 <- c(rnorm(200, 3, 0.5), rnorm(800, 0, 1))
  data2 <- c(rnorm(300, 3, 0.5), rnorm(700, 5, 2))
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=0, sd=1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=5, sd=2))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 100, compute_var = TRUE)

  expect_false(is.na(res$variance_est_p1))
  expect_false(is.na(res$variance_est_p2))
})


test_that("estim_IBM compute_var avec composantes connues égales", {
  skip_on_cran()

  set.seed(1)
  data1 <- c(rnorm(200, 3, 1), rnorm(800, 0, 1))
  data2 <- c(rnorm(300, 5, 1), rnorm(700, 0, 1))
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=0, sd=1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=0, sd=1))
  res_equal_comp <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 100, compute_var = TRUE)

  expect_true(res_equal_comp$equal.knownComp)
  expect_false(is.na(res_equal_comp$variance_est_p2))
})


test_that("print.estim_IBM affiche les variances quand disponibles", {
  skip_on_cran()

  set.seed(1)
  data1 <- c(rnorm(200, 3, 0.5), rnorm(800, 0, 1))
  data2 <- c(rnorm(300, 3, 0.5), rnorm(700, 5, 2))
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=0, sd=1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=5, sd=2))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 100, compute_var = TRUE)
  out <- capture.output(print.estim_IBM(res))

  expect_true(any(grepl("Variance", out)))
})


test_that("summary.estim_IBM couvre la branche equal.knownComp", {
  skip_on_cran()

  set.seed(1)
  data1 <- c(rnorm(200, 3, 1), rnorm(800, 0, 1))
  data2 <- c(rnorm(300, 5, 1), rnorm(700, 0, 1))
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=0, sd=1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=0, sd=1))
  res_equal_comp <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 100, compute_var = TRUE)
  out <- capture.output(summary.estim_IBM(res_equal_comp))
  expect_true(any(grepl("equal known components", out)))
})


test_that("summary.estim_IBM affiche les variances", {
  skip_on_cran()

  set.seed(1)
  data1 <- c(rnorm(200, 3, 0.5), rnorm(800, 0, 1))
  data2 <- c(rnorm(300, 3, 0.5), rnorm(700, 5, 2))
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=0, sd=1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean=5, sd=2))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 100, compute_var = TRUE)
  out <- capture.output(summary.estim_IBM(res))
  expect_true(any(grepl("Variance", out)))
})


test_that("IBM_hessian_contrast fonctionne avec support discret", {
  skip_on_cran()

  set.seed(1)
  # données multinom
  mixt1 <- twoComp_mixt(n = 500, weight = 0.8, comp.dist = list("multinom","multinom"),
                        comp.param = list(list(size=1, prob=c(0.2,0.3,0.5)),
                                          list(size=1, prob=c(0.1,0.6,0.3))))
  data1 <- get_mixture_data(mixt1)
  mixt2 <- twoComp_mixt(n = 500, weight = 0.3, comp.dist = list("multinom","multinom"),
                        comp.param = list(list(size=1, prob=c(0.2,0.3,0.5)),
                                          list(size=1, prob=c(0.7,0.1,0.2))))
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size=1, prob=c(0.1,0.6,0.3)))
  mod2 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size=1, prob=c(0.7,0.1,0.2)))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2))
  expect_s3_class(res, "estim_IBM")
})


test_that("IBM_normalization_term couvre le cas n1 <= n2 avec composantes distinctes", {
  skip_on_cran()

  set.seed(42)
  # n1 < n2 : échantillon 1 plus petit que l'échantillon 2
  mixt1 <- twoComp_mixt(n = 500, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 1000, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 5, sd = 2)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 5, sd = 2))
  # compute_var = TRUE est la seule façon d'atteindre IBM_normalization_term
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 200, compute_var = TRUE)

  expect_s3_class(res, "estim_IBM")
  expect_false(res$equal.knownComp)
  expect_true(length(res$population_sizes) == 2)
  expect_true(res$population_sizes[1] < res$population_sizes[2])  # garantit n1 < n2
  expect_false(is.na(res$variance_est_p1))
  expect_false(is.na(res$variance_est_p2))
  expect_true(res$variance_est_p1 >= 0)
  expect_true(res$variance_est_p2 >= 0)
})


test_that("IBM_normalization_term couvre le cas n1 <= n2 avec composantes connues égales", {
  skip_on_cran()

  set.seed(42)
  # Composantes connues identiques (même loi, mêmes paramètres)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 1000, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 5, sd = 2), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  # Composantes connues identiques : même distribution, mêmes paramètres
  mod1 <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 200, compute_var = TRUE)

  expect_true(res$equal.knownComp)                                  # G1 == G2
  expect_true(res$population_sizes[1] < res$population_sizes[2])   # n1 < n2
  expect_false(is.na(res$variance_est_p2))                         # calcul abouti
  expect_true(res$variance_est_p2 >= 0)
  # Dans ce cas, p1 est fixé arbitrairement à 0.2 et non estimé
  expect_equal(res$p.X.fixed, 0.2)
  expect_equal(res$variance_est_p1, 0)                          # p1 non estimé → variance nulle
})


# ==============================================================================
# Branches d'optimisation dans estim_IBM
# ==============================================================================

test_that("estim_IBM : retry Nelder-Mead", {
  skip_on_cran()

  # On force optim() à retourner une try-error au premier appel via mock,
  # puis un résultat avec par > 5 pour déclencher le switch vers L-BFGS-B, et finalement un résultat valide.
  call_count <- 0L
  local_mocked_bindings(
    optim = function(...) {
      call_count <<- call_count + 1L
      if (call_count <= 3L) stop("forced optim failure")   # déclenche try-error → retry NM
      list(par = c(0.5, 0.6))                              # résultat valide au 4e appel
    },
    .package = "stats"
  )
  d <- make_distinct_data()
  # Ne doit pas planter : les retries absorbent les erreurs
  res <- estim_IBM(samples = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ = 50)
  expect_s3_class(res, "estim_IBM")
})


test_that("estim_IBM : stop() final quand tous les algos échouent", {
  # optim() échoue systématiquement → after retries → stop()
  local_mocked_bindings(
    optim = function(...) stop("forced failure always"),
    .package = "stats"
  )
  d <- make_distinct_data()
  expect_error(
    estim_IBM(samples = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ = 50),
    "whatever the optimization algorithm"
  )
})


test_that("estim_IBM : retry L-BFGS-B et repli Nelder-Mead couverts", {
  # Premier appel NM : renvoie par > 5 → switch vers L-BFGS-B
  # L-BFGS-B échoue plusieurs fois → repli NM → résultat valide
  call_count <- 0L
  local_mocked_bindings(
    optim = function(par, method, ...) {
      call_count <<- call_count + 1L
      if (method == "Nelder-Mead" && call_count == 1L)
        return(list(par = c(10, 10)))       # par > 5 → déclenche L-BFGS-B
      if (method == "L-BFGS-B")
        stop("L-BFGS-B failure")            # L-BFGS-B échoue toujours
      list(par = c(0.5, 0.6))              # repli NM final réussit
    },
    .package = "stats"
  )
  d <- make_distinct_data()
  res <- estim_IBM(samples = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ = 50)
  expect_s3_class(res, "estim_IBM")
})


# ==============================================================================
# estim_IBM branche else NULL (dim inattendue de var_estimators)
# ==============================================================================

test_that("estim_IBM : branche else NULL quand dim(var_estimators) inattendue", {
  skip_on_cran()
  # On mock IBM_estimVarCov_gaussVect pour renvoyer une matrice 4x4
  # (ni 3x3 ni 2x2) → les deux if/else if ratent → else NULL → var_est_p1/p2 non assignés
  # Cela provoque une erreur à L.169 car var_est_p1 n'existe pas.
  local_mocked_bindings(
    IBM_estimVarCov_gaussVect = function(...) matrix(1, nrow = 4, ncol = 4),
    .package = "admix"
  )
  d <- make_distinct_data()
  # var_est_p1 et var_est_p2 ne seront pas créés → erreur à l'accès L.169
  expect_error(
    estim_IBM(samples = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), n.integ = 50, compute_var = TRUE),
    "object 'var_est_p1' not found"
  )
})


test_that("print.estim_IBM : affiche variance_est_p2 quand equal.knownComp = TRUE (L.206)", {
  skip_on_cran()

  set.seed(1)
  data1 <- c(rnorm(500, 3, 1), rnorm(500, 0, 1))
  data2 <- c(rnorm(500, 5, 1), rnorm(500, 0, 1))
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 100, compute_var = TRUE)

  expect_true(res$equal.knownComp)
  expect_false(is.na(res$variance_est_p2))

  out <- capture.output(print.estim_IBM(res))
  expect_true(any(grepl("Variance.*second", out, ignore.case = TRUE)))
})


test_that("IBM_normalization_term et IBM_mat_L : n1 > n2 avec composantes égales", {
  skip_on_cran()

  set.seed(42)
  # n1 > n2 : on inverse les tailles
  mixt1 <- twoComp_mixt(n = 1000, weight = 0.4, comp.dist = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 500,  weight = 0.6, comp.dist = list("norm", "norm"),
                        comp.param = list(list(mean = 5, sd = 2),  list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), n.integ = 100, compute_var = TRUE)

  expect_true(res$equal.knownComp)
  expect_true(res$population_sizes[1] > res$population_sizes[2])  # garantit n1 > n2
  expect_false(is.na(res$variance_est_p2))
  expect_equal(res$p.X.fixed, 0.2)
})


test_that("estim_IBM avec compute_var=TRUE couvre les branches discrètes de IBM_Sigma1/2 et IBM_hessian_contrast", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.8, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.3, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.7, 0.1, 0.2))))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))
  mod2 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, prob = c(0.7, 0.1, 0.2)))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), compute_var = TRUE)

  expect_s3_class(res, "estim_IBM")
  expect_false(res$equal.knownComp)
  expect_false(is.na(res$variance_est_p1))
  expect_false(is.na(res$variance_est_p2))
})


test_that("estim_IBM avec compute_var=TRUE couvre les branches discrètes avec composantes égales (G1=G2)", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.3, 0.4, 0.3)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  mixt2 <- twoComp_mixt(n = 800, weight = 0.4, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.5, 0.3, 0.2)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))
  mod2 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))
  res <- estim_IBM(samples = list(data1, data2), admixMod = list(mod1, mod2), compute_var = TRUE)

  expect_true(res$equal.knownComp)
  expect_false(is.na(res$variance_est_p2))
})


# ==============================================================================
# IBM_gap : branche stepfun (pmultinom) et branche G1equalG2 = TRUE
# ==============================================================================

test_that("IBM_gap : branche stepfun couverte avec données multinom et G1 != G2", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.8, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.3, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.7, 0.1, 0.2))))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))
  mod2 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, prob = c(0.7, 0.1, 0.2)))

  # G1 != G2 : fixed.p.X = NULL → branche is.null(fixed.p.X) = TRUE
  expect_message(
    gap <- IBM_gap(z = 2, par = c(0.5, 0.3), samples = list(data1, data2), admixMod = list(mod1, mod2)),
    "fixed.p.X"
  )
  expect_true(is.numeric(gap))
})


test_that("IBM_gap : branche G1equalG2 = TRUE avec stepfun", {
  set.seed(1)
  # Composantes connues identiques + multinomial → stepfun + G1=G2
  mixt1 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.3, 0.4, 0.3)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.4, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.5, 0.3, 0.2)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))
  mod2 <- admix_model(knownComp_dist = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))

  expect_message(
    gap <- IBM_gap(z = 2, par = 0.4, samples = list(data1, data2), admixMod = list(mod1, mod2)),
    "fixed.p.X"
  )
  expect_true(is.numeric(gap))
})


test_that("IBM_theoretical_gap : branche length(par) != 2", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 400, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  res <- IBM_theoretical_gap(z = 1.0, par = 0.5, known.p = c(0.5, 0.7), mixtMod = list(mixt1, mixt2))

  expect_true(is.numeric(res))
  expect_length(res, 1)
})


test_that("IBM_theoretical_gap : branche stepfun avec données multinom", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.8, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.3, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.7, 0.1, 0.2))))
  res <- IBM_theoretical_gap(z = 2, par = c(0.8, 0.3), known.p = c(0.8, 0.3), mixtMod = list(mixt1, mixt2))

  expect_true(is.numeric(res))
})
