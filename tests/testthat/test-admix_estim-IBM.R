# =============================================================================
# Helpers partagés : données simulées réutilisées dans tous les tests
# =============================================================================

## Cas continu — composantes connues DIFFÉRENTES (cas général)
setup_continuous_diff <- function() {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.5,
                        comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 600, weight = 0.7,
                        comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 5, sd = 2)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  admixMod1 <- admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]])
  admixMod2 <- admix_model(knownComp_dist  = mixt2$comp.dist[[2]], knownComp_param = mixt2$comp.param[[2]])
  list(data1 = data1, data2 = data2, admixMod1 = admixMod1, admixMod2 = admixMod2, true_p1 = 0.5, true_p2 = 0.7)
}

## Cas continu — composantes connues IDENTIQUES (optimisation 1D)
setup_continuous_equal <- function() {
  set.seed(123)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.6,
                        comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 600, weight = 0.4,
                        comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 5, sd = 2), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  admixMod1 <- admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]])
  admixMod2 <- admix_model(knownComp_dist  = mixt2$comp.dist[[2]], knownComp_param = mixt2$comp.param[[2]])
  list(data1 = data1, data2 = data2, admixMod1 = admixMod1, admixMod2 = admixMod2, true_p2 = 0.4)
}

## Cas discret — loi multinomiale
setup_discrete <- function() {
  set.seed(7)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.8,
                        comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  mixt2 <- twoComp_mixt(n = 600, weight = 0.3,
                        comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.7, 0.1, 0.2))))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  admixMod1 <- admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]])
  admixMod2 <- admix_model(knownComp_dist  = mixt2$comp.dist[[2]], knownComp_param = mixt2$comp.param[[2]])
  list(data1 = data1, data2 = data2, admixMod1 = admixMod1, admixMod2 = admixMod2, true_p1 = 0.8, true_p2 = 0.3)
}


# =============================================================================
# 1. Tests de validation des entrées (estim_IBM)
# =============================================================================

test_that("estim_IBM — erreur si nombre de samples != 2", {
  ctx <- setup_continuous_diff()
  expect_error(
    estim_IBM(samples  = list(ctx$data1), admixMod = list(ctx$admixMod1, ctx$admixMod2)),
    regexp = "Wrong number of samples"
  )
})

test_that("estim_IBM — erreur si admixMod n'est pas de classe admix_model", {
  ctx <- setup_continuous_diff()
  expect_error(
    estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list("not_a_model", ctx$admixMod2)),
    regexp = "admixMod"
  )
})

test_that("estim_IBM — erreur si les deux admixMod sont mal spécifiés", {
  ctx <- setup_continuous_diff()
  expect_error(
    estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(list(a = 1), list(b = 2))),
    regexp = "admixMod"
  )
})


# =============================================================================
# 2. Tests de structure de sortie (estim_IBM)
# =============================================================================

test_that("estim_IBM — retourne un objet de classe estim_IBM et admix_estim", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_s3_class(res, "estim_IBM")
  expect_s3_class(res, "admix_estim")
})

test_that("estim_IBM — la liste de résultats contient tous les attributs attendus", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expected_names <- c("n_populations", "population_sizes", "admixture_models",
                      "estimation_method", "estimated_mixing_weights",
                      "variance_est_p1", "variance_est_p2",
                      "equal.knownComp", "p.X.fixed", "integ.supp", "call")
  expect_true(all(expected_names %in% names(res)))
})

test_that("estim_IBM — n_populations vaut 2", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_equal(res$n_populations, 2L)
})

test_that("estim_IBM — population_sizes correspond aux tailles réelles des samples", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_equal(res$population_sizes, c(length(ctx$data1), length(ctx$data2)))
})

test_that("estim_IBM — estimation_method est 'Inversion Best Matching (IBM)'", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_equal(res$estimation_method, "Inversion Best Matching (IBM)")
})

test_that("estim_IBM — equal.knownComp est FALSE pour des composantes différentes", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_false(res$equal.knownComp)
})

test_that("estim_IBM — equal.knownComp est TRUE pour des composantes identiques", {
  ctx <- setup_continuous_equal()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_true(res$equal.knownComp)
})

test_that("estim_IBM — p.X.fixed vaut 0.2 quand composantes égales", {
  ctx <- setup_continuous_equal()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_equal(res$p.X.fixed, 0.2)
})

test_that("estim_IBM — p.X.fixed est NULL quand composantes différentes", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_null(res$p.X.fixed)
})

test_that("estim_IBM — integ.supp est un vecteur numérique trié", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_true(is.numeric(res$integ.supp))
  expect_true(all(diff(res$integ.supp) >= 0))
})

test_that("estim_IBM — variance NA quand compute_var = FALSE (défaut)", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  expect_true(is.na(res$variance_est_p1))
  expect_true(is.na(res$variance_est_p2))
})


# =============================================================================
# 3. Tests de cohérence statistique des estimations (cas continu, composantes differentes)
# =============================================================================

test_that("estim_IBM — poids estimés sont des scalaires numériques dans ]0,1[", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 300)
  w <- res$estimated_mixing_weights
  expect_true(is.numeric(w))
  expect_length(w, 2)
  expect_true(all(w > 0 & w < 1))
})

test_that("estim_IBM — poids estimés proches des vraies valeurs (tolérance 15 %)", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 300)
  w <- res$estimated_mixing_weights
  expect_equal(w[1], ctx$true_p1, tolerance = 0.15)
  expect_equal(w[2], ctx$true_p2, tolerance = 0.15)
})


# =============================================================================
# 4. Tests cas discret (loi multinomiale)
# =============================================================================

test_that("estim_IBM — fonctionne pour le cas discret (multinomial)", {
  ctx <- setup_discrete()
  expect_no_error(estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2)))
})

test_that("estim_IBM — classe correcte en cas discret", {
  ctx <- setup_discrete()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2))
  expect_s3_class(res, "estim_IBM")
})

test_that("estim_IBM — poids estimés dans ]0,1[ en cas discret", {
  ctx <- setup_discrete()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2))
  w <- res$estimated_mixing_weights
  expect_true(all(w > 0 & w < 1))
})


# =============================================================================
# Helper : capture la sortie console de cat() ET print() de façon fiable
# =============================================================================

capture_cat <- function(expr) {
  paste(capture.output(expr, type = "output"), collapse = "\n")
}

# =============================================================================
# 5. Tests de print.estim_IBM
# =============================================================================

test_that("print.estim_IBM — retourne invisiblement l'objet", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  ret <- withVisible(print(res))
  expect_identical(ret$value, res)
  expect_false(ret$visible)
})

test_that("print.estim_IBM — affiche la section 'Method: IBM'", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  out <- capture.output(print(res))
  expect_true(any(grepl("IBM", out)))
})

test_that("print.estim_IBM — la sortie contient le nom de la fonction appelée", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  out <- paste(capture.output(print(res)), collapse = " ")
  expect_true(grepl("estim_IBM", out))
})

test_that("print.estim_IBM — produit au moins une ligne non vide", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  out <- capture.output(print(res))
  expect_true(any(nchar(trimws(out)) > 0))
})


# =============================================================================
# 6. Tests de summary.estim_IBM
# =============================================================================

test_that("summary.estim_IBM — retourne invisiblement l'objet", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  ret <- withVisible(summary(res))
  expect_identical(ret$value, res)
  expect_false(ret$visible)
})

test_that("summary.estim_IBM — produit une sortie non vide", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  out <- capture.output(summary(res))
  expect_true(any(nchar(trimws(out)) > 0))
})

test_that("summary.estim_IBM — la sortie mentionne IBM", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  out <- paste(capture.output(summary(res)), collapse = " ")
  expect_true(grepl("IBM", out))
})

test_that("summary.estim_IBM — la sortie contient le Call", {
  ctx <- setup_continuous_diff()
  res <- estim_IBM(samples  = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2), n.integ  = 200)
  out <- paste(capture.output(summary(res)), collapse = " ")
  expect_true(grepl("Call|estim_IBM", out))
})


# =============================================================================
# 7. Tests de IBM_empirical_contrast
# =============================================================================

test_that("IBM_empirical_contrast — retourne un scalaire numérique non-négatif", {
  ctx <- setup_continuous_diff()
  set.seed(1)
  G <- stats::runif(300, min = min(c(ctx$data1, ctx$data2)), max = max(c(ctx$data1, ctx$data2)))
  val <- IBM_empirical_contrast(par = c(0.5, 0.7), samples = list(ctx$data1, ctx$data2),
                                admixMod = list(ctx$admixMod1, ctx$admixMod2), G = G, fixed.p.X = NULL)
  expect_length(val, 1)
  expect_true(is.numeric(val))
  expect_gte(val, 0)
})

test_that("IBM_empirical_contrast — contraste minimal proche de 0 aux vraies valeurs", {
  ctx <- setup_continuous_diff()
  set.seed(2)
  G <- stats::runif(300, min = min(c(ctx$data1, ctx$data2)), max = max(c(ctx$data1, ctx$data2)))
  val_true  <- IBM_empirical_contrast(par = c(ctx$true_p1, ctx$true_p2), samples  = list(ctx$data1, ctx$data2),
                                      admixMod = list(ctx$admixMod1, ctx$admixMod2), G = G, fixed.p.X = NULL)
  val_wrong <- IBM_empirical_contrast(par = c(0.1, 0.9), samples  = list(ctx$data1, ctx$data2),
                                      admixMod = list(ctx$admixMod1, ctx$admixMod2), G = G, fixed.p.X = NULL)
  expect_lt(val_true, val_wrong)
})

test_that("IBM_empirical_contrast — erreur si fixed.p.X NULL et G1 = G2", {
  ctx <- setup_continuous_equal()
  set.seed(3)
  G <- stats::runif(200, min = min(c(ctx$data1, ctx$data2)), max = max(c(ctx$data1, ctx$data2)))
  expect_error(
    IBM_empirical_contrast(par = 0.4, samples = list(ctx$data1, ctx$data2),
                           admixMod = list(ctx$admixMod1, ctx$admixMod2), G = G, fixed.p.X = NULL)
  )
})


# =============================================================================
# 8. Tests de IBM_gap
# =============================================================================

test_that("IBM_gap — retourne un scalaire numérique", {
  ctx <- setup_continuous_diff()
  val <- IBM_gap(z = 2.0, par = c(0.5, 0.7), samples = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2))
  expect_length(val, 1)
  expect_true(is.numeric(val))
})

test_that("IBM_gap — gap proche de 0 aux vraies valeurs (signe non garanti)", {
  ctx <- setup_continuous_diff()
  val <- IBM_gap(z = median(ctx$data1), par = c(ctx$true_p1, ctx$true_p2), samples = list(ctx$data1, ctx$data2),
                 admixMod = list(ctx$admixMod1, ctx$admixMod2))
  expect_lt(abs(val), 0.5)   # valeur raisonnable pour des échantillons de taille 500-600
})

test_that("IBM_gap — cas G1=G2 : utilise fixed.p.X = 0.2 et retourne un scalaire", {
  ctx <- setup_continuous_equal()

  expect_message(
    val <- IBM_gap(z = 0, par = 0.4, samples = list(ctx$data1, ctx$data2), admixMod = list(ctx$admixMod1, ctx$admixMod2)),
    regexp = "fixed.p.X"
  )
  expect_length(val, 1)
  expect_true(is.numeric(val))
})


# =============================================================================
# 9. Tests de IBM_theoretical_gap
# =============================================================================

test_that("IBM_theoretical_gap — retourne un scalaire numérique", {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 600, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 5, sd = 2)))
  val <- IBM_theoretical_gap(z = 2.8, par = c(0.5, 0.7), known.p = c(0.5, 0.7), mixtMod = list(mixt1, mixt2))
  expect_length(val, 1)
  expect_true(is.numeric(val))
})

test_that("IBM_theoretical_gap — gap théorique nul aux vraies valeurs (identifiability)", {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 600, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 5, sd = 2)))
  # Quand par = known.p, les composantes inconnues sont identiques => gap = 0
  val <- IBM_theoretical_gap(z = 1.0, par = c(0.5, 0.7), known.p = c(0.5, 0.7), mixtMod = list(mixt1, mixt2))
  expect_equal(val, 0, tolerance = 1e-10)
})

test_that("IBM_theoretical_gap — erreur si known.p est NULL", {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 200, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 200, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 5, sd = 2)))
  expect_error(
    IBM_theoretical_gap(z = 1.0, par = c(0.5, 0.7), known.p = NULL, mixtMod = list(mixt1, mixt2)),
    regexp = "known.p"
  )
})

test_that("IBM_theoretical_gap — erreur si mixtMod n'est pas de classe twoComp_mixt", {
  expect_error(
    IBM_theoretical_gap(z = 1.0, par = c(0.5, 0.7), known.p = c(0.5, 0.7), mixtMod = list(list(a = 1), list(b = 2))),
    regexp = "mixtMod"
  )
})
