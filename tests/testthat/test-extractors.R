# ------------------------------------------------------------------------------
# Helpers : objets réutilisés dans plusieurs blocs
# ------------------------------------------------------------------------------

## Données simulées de base (deux mélanges gaussiens)
make_base_mixtures <- function(seed = 42) {
  set.seed(seed)
  mixt1 <- twoComp_mixt(n = 380, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 350, weight = 0.85, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = -1, sd = 1)))
  list(mixt1 = mixt1,
       mixt2 = mixt2,
       data1 = get_mixture_data(mixt1),
       data2 = get_mixture_data(mixt2),
       admixMod1 = admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]]),
       admixMod2 = admix_model(knownComp_dist  = mixt2$comp.dist[[2]], knownComp_param = mixt2$comp.param[[2]]))
}

## Objet admix_estim (une population, méthode BVdk)
make_admix_estim <- function(seed = 42) {
  d <- make_base_mixtures(seed)
  admix_estim(samples = list(d$data1), admixMod = list(d$admixMod1), est_method = "BVdk")
}

## Objet admix_test (deux populations, méthode poly)
make_admix_test_poly <- function(seed = 42) {
  d <- make_base_mixtures(seed)
  admix_test(samples = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), conf_level = 0.95,
             test_method = "poly", ask_poly_param = FALSE, support = "Real")
}

## Objet admix_test (deux populations, méthode icv — fournit tabulated_dist)
make_admix_test_icv <- function(seed = 42) {
  d <- make_base_mixtures(seed)
  admix_test(samples = list(d$data1, d$data2), admixMod = list(d$admixMod1, d$admixMod2), conf_level = 0.95,
             test_method = "icv", n_sim_tab = 20)
}

# ==============================================================================
# 1. get_mixture_data
# ==============================================================================

test_that("get_mixture_data : retourne un vecteur numérique", {
  d <- make_base_mixtures()
  expect_true(is.numeric(d$data1))
  expect_true(is.vector(d$data1))
})

test_that("get_mixture_data : longueur égale à n", {
  set.seed(1)
  mixt <- twoComp_mixt(n = 200, weight = 0.5, comp.dist  = list("norm", "norm"),
                       comp.param = list(list(mean = 0, sd = 1), list(mean = 2, sd = 1)))
  expect_equal(length(get_mixture_data(mixt)), 200L)
})

test_that("get_mixture_data : erreur sur un objet de mauvaise classe", {
  # UseMethod sans méthode applicable lève une erreur générique de dispatch —
  # on vérifie uniquement qu'une erreur est bien déclenchée, sans matcher le message.
  expect_error(get_mixture_data(list(a = 1)))
})


# ==============================================================================
# 2. get_known_component
# ==============================================================================

test_that("get_known_component.admix_estim : retourne une liste", {
  skip_on_cran()
  est <- make_admix_estim()
  res <- get_known_component(est)
  expect_true(is.list(res))
})


test_that("get_known_component.admix_estim : longueur = nombre de populations", {
  skip_on_cran()
  est <- make_admix_estim()
  res <- get_known_component(est)
  expect_equal(length(res), 1L)
})


test_that("get_known_component.admix_estim : erreur sur mauvaise classe", {
  expect_error(get_known_component.admix_estim(list()), regexp = "admix_estim")
})


test_that("get_known_component.IBM_test : retourne les modèles admixture", {
  skip_on_cran()
  tst <- make_admix_test_icv()
  # IBM_test est produit par admix_test avec méthode icv
  if (inherits(tst, "IBM_test")) {
    res <- get_known_component(tst)
    expect_true(!is.null(res))
  } else {
    skip("L'objet retourné n'est pas de classe IBM_test avec ces paramètres")
  }
})


test_that("get_known_component.IBM_test : erreur sur mauvaise classe", {
  expect_error(get_known_component.IBM_test(list()), regexp = "IBM_test")
})


test_that("get_known_component.orthobasis_test : erreur sur mauvaise classe", {
  expect_error(get_known_component.orthobasis_test(list()), regexp = "orthobasis_test")
})

test_that("get_known_component.gaussianity_test : erreur sur mauvaise classe", {
  expect_error(get_known_component.gaussianity_test(list()), regexp = "gaussianity_test")
})

test_that("get_known_component.admix_cluster : erreur sur mauvaise classe", {
  expect_error(get_known_component.admix_cluster(list()), regexp = "admix_cluster")
})


# ==============================================================================
# 3. get_mixing_weights
# ==============================================================================

test_that("get_mixing_weights.admix_estim : retourne un vecteur numérique", {
  skip_on_cran()
  est <- make_admix_estim()
  res <- get_mixing_weights(est)
  expect_true(is.numeric(res))
})

test_that("get_mixing_weights.admix_estim : une valeur par population", {
  skip_on_cran()
  est <- make_admix_estim()
  res <- get_mixing_weights(est)
  expect_length(res, 1L)
})

test_that("get_mixing_weights.admix_estim : poids dans [0, 1]", {
  skip_on_cran()
  est <- make_admix_estim()
  res <- get_mixing_weights(est)
  expect_true(all(res >= 0 & res <= 1))
})

test_that("get_mixing_weights.admix_estim : erreur sur mauvaise classe", {
  expect_error(get_mixing_weights.admix_estim(list()), regexp = "admix_estim")
})

test_that("get_mixing_weights.gaussianity_test : erreur sur mauvaise classe", {
  expect_error(get_mixing_weights.gaussianity_test(list()), regexp = "gaussianity_test")
})

test_that("get_mixing_weights.orthobasis_test : erreur sur mauvaise classe", {
  expect_error(get_mixing_weights.orthobasis_test(list()), regexp = "orthobasis_test")
})


# ==============================================================================
# 4. reject_nullHyp
# ==============================================================================

test_that("reject_nullHyp : retourne un booléen", {
  skip_on_cran()
  tst <- make_admix_test_poly()
  res <- reject_nullHyp(tst)
  expect_true(is.logical(res))
  expect_length(res, 1L)
})

test_that("reject_nullHyp.gaussianity_test : erreur sur mauvaise classe", {
  expect_error(reject_nullHyp.gaussianity_test(list()), regexp = "gaussianity_test")
})

test_that("reject_nullHyp.orthobasis_test : erreur sur mauvaise classe", {
  expect_error(reject_nullHyp.orthobasis_test(list()), regexp = "orthobasis_test")
})

test_that("reject_nullHyp.IBM_test : erreur sur mauvaise classe", {
  expect_error(reject_nullHyp.IBM_test(list()), regexp = "IBM_test")
})


# ==============================================================================
# 5. which_rank
# ==============================================================================

test_that("which_rank : retourne un scalaire numérique entier positif", {
  skip_on_cran()
  tst <- make_admix_test_poly()
  res <- which_rank(tst)
  expect_true(is.numeric(res))
  expect_length(res, 1L)
  expect_gte(res, 1)
})

test_that("which_rank.gaussianity_test : erreur sur mauvaise classe", {
  expect_error(which_rank.gaussianity_test(list()), regexp = "gaussianity_test")
})

test_that("which_rank.orthobasis_test : erreur sur mauvaise classe", {
  expect_error(which_rank.orthobasis_test(list()), regexp = "orthobasis_test")
})

test_that("which_rank.IBM_test : erreur sur mauvaise classe", {
  expect_error(which_rank.IBM_test(list()), regexp = "IBM_test")
})


# ==============================================================================
# 6. get_tabulated_dist
# ==============================================================================

test_that("get_tabulated_dist.IBM_test : retourne un vecteur numérique trié", {
  skip_on_cran()
  tst <- make_admix_test_icv()
  if (inherits(tst, "IBM_test")) {
    res <- get_tabulated_dist(tst)
    expect_true(is.numeric(res))
    expect_true(all(diff(res) >= 0))   # trié en ordre croissant
  } else {
    skip("L'objet n'est pas de classe IBM_test")
  }
})

test_that("get_tabulated_dist.IBM_test : erreur sur mauvaise classe", {
  expect_error(get_tabulated_dist.IBM_test(list()), regexp = "IBM_test")
})

test_that("get_tabulated_dist.admix_cluster : erreur sur mauvaise classe", {
  expect_error(get_tabulated_dist.admix_cluster(list()), regexp = "admix_cluster")
})


# ==============================================================================
# 7. get_discrepancy_rank  &  get_discrepancy_matrix  (test à 2 échantillons)
# ==============================================================================

test_that("get_discrepancy_rank.IBM_test : affiche un message pour le cas 2 échantillons", {
  skip_on_cran()
  tst <- make_admix_test_icv()
  if (inherits(tst, "IBM_test")) {
    # En test à 2 échantillons, discrepancy_rank est NA → message affiché
    out <- capture.output(get_discrepancy_rank(tst))
    expect_true(any(grepl("Two-sample|matrix|rank", out, ignore.case = TRUE)))
  } else {
    skip("L'objet n'est pas de classe IBM_test")
  }
})

test_that("get_discrepancy_rank.IBM_test : erreur sur mauvaise classe", {
  expect_error(get_discrepancy_rank.IBM_test(list()), regexp = "IBM_test")
})

test_that("get_discrepancy_matrix.IBM_test : affiche un message pour le cas 2 échantillons", {
  skip_on_cran()
  tst <- make_admix_test_icv()
  if (inherits(tst, "IBM_test")) {
    out <- capture.output(get_discrepancy_matrix(tst))
    expect_true(any(grepl("Two-sample|matrix", out, ignore.case = TRUE)))
  } else {
    skip("L'objet n'est pas de classe IBM_test")
  }
})

test_that("get_discrepancy_matrix.IBM_test : erreur sur mauvaise classe", {
  expect_error(get_discrepancy_matrix.IBM_test(list()), regexp = "IBM_test")
})

test_that("get_discrepancy_matrix.admix_cluster : erreur sur mauvaise classe", {
  expect_error(get_discrepancy_matrix.admix_cluster(list()), regexp = "admix_cluster")
})


# ==============================================================================
# 8. get_statistic_components
# ==============================================================================

test_that("get_statistic_components.IBM_test : affiche un message pour le cas 2 échantillons", {
  skip_on_cran()
  tst <- make_admix_test_icv()
  if (inherits(tst, "IBM_test")) {
    # En test à 2 échantillons, statistic_name peut être NA
    out <- tryCatch(get_statistic_components(tst), error = function(e) NULL)
    expect_true(!is.null(out))
  } else {
    skip("L'objet n'est pas de classe IBM_test")
  }
})

test_that("get_statistic_components.IBM_test : erreur sur mauvaise classe", {
  expect_error(get_statistic_components.IBM_test(list()), regexp = "IBM_test")
})


# ==============================================================================
# 9. get_cluster_members  &  get_cluster_sizes
#    (tests sur la gestion d'erreur uniquement — clustering coûteux à simuler)
# ==============================================================================

test_that("get_cluster_members.admix_cluster : erreur sur mauvaise classe", {
  expect_error(get_cluster_members.admix_cluster(list()), regexp = "admix_cluster")
})

test_that("get_cluster_sizes.admix_cluster : erreur sur mauvaise classe", {
  expect_error(get_cluster_sizes.admix_cluster(list()), regexp = "admix_cluster")
})


# ------------------------------------------------------------------------------
# Tests d'intégration légers sur get_cluster_members / get_cluster_sizes
# exécutés uniquement si un objet admix_cluster est disponible rapidement
# ------------------------------------------------------------------------------

test_that("get_cluster_members : retourne une liste si objet valide", {
  # Construction d'un faux objet de classe admix_cluster pour tester la structure
  fake_cluster <- structure(
    list(clusters = list(c(1, 2), c(3)), clust_sizes = c(2L, 1L), admixture_models = list(),
         discrepancy_matrix = matrix(NA), tab_distributions  = list()),
    class = "admix_cluster"
  )
  res <- get_cluster_members(fake_cluster)
  expect_true(is.list(res))
  expect_length(res, 2L)
})

test_that("get_cluster_sizes : retourne un vecteur entier si objet valide", {
  fake_cluster <- structure(
    list(clusters = list(c(1, 2), c(3)), clust_sizes = c(2L, 1L), admixture_models = list(),
         discrepancy_matrix = matrix(NA), tab_distributions  = list()),
    class = "admix_cluster"
  )
  res <- get_cluster_sizes(fake_cluster)
  expect_true(is.numeric(res) || is.integer(res))
  expect_equal(sum(res), 3L)
})


test_that("get_known_component generic and methods work", {

  expect_error(get_known_component(1), "no applicable method")

  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))

  obj_estim <- list(estim_objects = list(list(admixture_models = admixMod)))
  class(obj_estim) <- "admix_estim"
  expect_equal(get_known_component(obj_estim)[[1]], admixMod)

  obj_gauss <- list(admixture_models = admixMod)
  class(obj_gauss) <- "gaussianity_test"
  expect_equal(get_known_component(obj_gauss), admixMod)

  obj_ortho <- list(admixture_models = admixMod)
  class(obj_ortho) <- "orthobasis_test"
  expect_equal(get_known_component(obj_ortho), admixMod)

  obj_ibm <- list(admixture_models = admixMod)
  class(obj_ibm) <- "IBM_test"
  expect_equal(get_known_component(obj_ibm), admixMod)

  obj_cluster <- list(admixture_models = admixMod)
  class(obj_cluster) <- "admix_cluster"
  expect_equal(get_known_component(obj_cluster), admixMod)
})


test_that("get_known_component methods reject wrong classes", {
  wrong <- structure(list(), class = "wrong")
  expect_error(get_known_component.admix_estim(wrong), "admix_estim")
  expect_error(get_known_component.gaussianity_test(wrong), "gaussianity_test")
  expect_error(get_known_component.orthobasis_test(wrong), "orthobasis_test")
  expect_error(get_known_component.IBM_test(wrong), "IBM_test")
  expect_error(get_known_component.admix_cluster(wrong), "admix_cluster")
})


test_that("get_mixture_data works", {
  expect_error(get_mixture_data(1), "no applicable method")

  obj <- list(mixt.data = 1:5)
  class(obj) <- "twoComp_mixt"
  expect_equal(get_mixture_data(obj), 1:5)

  wrong <- structure(list(), class = "wrong")
  expect_error(get_mixture_data.twoComp_mixt(wrong), "twoComp_mixt")
})


test_that("get_mixing_weights methods work", {
  obj_estim <- list(
    estim_objects = list(list(estimated_mixing_weights = 0.2), list(estimated_mixing_weights = 0.8))
  )
  class(obj_estim) <- "admix_estim"
  expect_equal(get_mixing_weights(obj_estim), c(0.2, 0.8))

  obj_gauss <- list(estimate = c("Weight 1" = 0.4, "other" = 2))
  class(obj_gauss) <- "gaussianity_test"
  expect_equal(get_mixing_weights(obj_gauss), 0.4)

  obj_ortho <- list(weights = c("Weight A" = 0.7, "x" = 3))
  class(obj_ortho) <- "orthobasis_test"
  expect_equal(get_mixing_weights(obj_ortho), 0.7)
})


test_that("reject_nullHyp methods work", {
  obj1 <- list(reject_decision = TRUE)
  class(obj1) <- "gaussianity_test"
  expect_true(reject_nullHyp(obj1))

  obj2 <- list(reject_decision = FALSE)
  class(obj2) <- "orthobasis_test"
  expect_false(reject_nullHyp(obj2))

  obj3 <- list(reject_decision = TRUE)
  class(obj3) <- "IBM_test"
  expect_true(reject_nullHyp(obj3))
})


test_that("which_rank methods work including NA branch", {
  obj1 <- list(selected_rank = 3)
  class(obj1) <- "gaussianity_test"
  expect_equal(which_rank(obj1), 3)

  obj2 <- list(selected_rank = 5)
  class(obj2) <- "orthobasis_test"
  expect_equal(which_rank(obj2), 5)

  obj3 <- list(selected_rank = 2)
  class(obj3) <- "IBM_test"
  expect_equal(which_rank(obj3), 2)

  obj4 <- list(selected_rank = NA)
  class(obj4) <- "IBM_test"
  expect_output(which_rank(obj4), "Two-sample test")
})


test_that("get_tabulated_dist methods work", {
  obj1 <- list(tabulated_dist = c(3,1,2))
  class(obj1) <- "IBM_test"
  expect_equal(get_tabulated_dist(obj1), c(1,2,3))

  obj2 <- list(tab_distributions = list(1:3))
  class(obj2) <- "admix_cluster"
  expect_equal(get_tabulated_dist(obj2), list(1:3))
})


test_that("IBM extractors cover NA branches", {
  obj_rank <- list(discrepancy_rank = matrix(NA, 2, 2))
  class(obj_rank) <- "IBM_test"
  expect_output(get_discrepancy_rank(obj_rank), "Two-sample test")

  obj_disc <- list(discrepancy_matrix = matrix(NA, 2, 2))
  class(obj_disc) <- "IBM_test"
  expect_output(get_discrepancy_matrix(obj_disc), "Two-sample test")

  obj_stat <- list(statistic_name = NA)
  class(obj_stat) <- "IBM_test"
  expect_output(get_statistic_components(obj_stat), "Two-sample test")
})


test_that("IBM extractors return values", {
  obj1 <- list(discrepancy_rank = matrix(1:4, 2))
  class(obj1) <- "IBM_test"
  expect_equal(get_discrepancy_rank(obj1), matrix(1:4, 2))

  obj2 <- list(discrepancy_matrix = matrix(5:8, 2))
  class(obj2) <- "IBM_test"
  expect_equal(get_discrepancy_matrix(obj2), matrix(5:8, 2))

  obj3 <- list(statistic_name = "T1")
  class(obj3) <- "IBM_test"
  expect_equal(get_statistic_components(obj3), "T1")
})


test_that("cluster extractors work", {
  obj <- list(
    clusters = matrix(1:4, 2),
    clust_sizes = c(2,2),
    discrepancy_matrix = matrix(1:9, 3)
  )
  class(obj) <- "admix_cluster"
  expect_equal(get_cluster_members(obj), matrix(1:4, 2))
  expect_equal(get_cluster_sizes(obj), c(2,2))
  expect_equal(get_discrepancy_matrix(obj), matrix(1:9, 3))
})
