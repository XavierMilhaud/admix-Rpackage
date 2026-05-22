
# =============================================================================
# 1. get_distribution_parameters()
# =============================================================================

test_that("get_distribution_parameters — loi normale : retourne 'mean' et 'sd'", {
  params <- get_distribution_parameters("norm")
  expect_true(is.character(params))
  expect_true("mean" %in% params)
  expect_true("sd"   %in% params)
})

test_that("get_distribution_parameters — loi gamma : retourne 'shape' et 'rate'/'scale'", {
  params <- get_distribution_parameters("gamma")
  expect_true(is.character(params))
  expect_true("shape" %in% params)
  # 'rate' ou 'scale' selon la version de R, les deux sont valides
  expect_true(any(c("rate", "scale") %in% params))
})

test_that("get_distribution_parameters — cas spécial 'multinom' : retourne 'size' et 'prob'", {
  params <- get_distribution_parameters("multinom")
  expect_setequal(params, c("size", "prob"))
})

test_that("get_distribution_parameters — cas spécial 'gompertz' : retourne 'shape' et 'rate'", {
  params <- get_distribution_parameters("gompertz")
  expect_setequal(params, c("shape", "rate"))
})

test_that("get_distribution_parameters — loi de Poisson : retourne 'lambda'", {
  params <- get_distribution_parameters("pois")
  expect_true("lambda" %in% params)
})

test_that("get_distribution_parameters — loi exponentielle : retourne 'rate'", {
  params <- get_distribution_parameters("exp")
  expect_true("rate" %in% params)
})

test_that("get_distribution_parameters — loi inconnue : lève une erreur", {
  expect_error(
    get_distribution_parameters("xyz_loi_inexistante"),
    regexp = "Unrecognized distribution"
  )
})

test_that("get_distribution_parameters — retourne un vecteur de chaînes de caractères", {
  params <- get_distribution_parameters("norm")
  expect_type(params, "character")
  expect_gte(length(params), 1L)
})


# =============================================================================
# 2. validate_distribution()
# =============================================================================

test_that("validate_distribution accepts valid normal distribution", {
  expect_invisible(
    validate_distribution(dist = "norm", params = list(mean = 0, sd = 1))
  )
})

test_that("validate_distribution rejects unknown distributions", {
  expect_error(
    validate_distribution(dist = "foobar", params = list(a = 1)),
    "Unrecognized distribution: foobar"
  )
})

test_that("validate_distribution rejects invalid parameter names", {
  expect_error(
    validate_distribution(dist = "norm", params = list(mu = 0, sigma = 1)),
    "Invalid parameter names for distribution `norm`"
  )
})

test_that("validate_distribution rejects missing parameters", {
  expect_error(
    validate_distribution(dist = "norm", params = list(mean = 0)),
    "Invalid parameter names"
  )
})

test_that("validate_distribution rejects extra parameters", {
  expect_error(
    validate_distribution(dist = "norm", params = list(mean = 0, sd = 1, extra = 2)),
    "Invalid parameter names"
  )
})

test_that("validate_distribution accepts gamma alternative parameterizations", {
  expect_invisible(
    validate_distribution(dist = "gamma", params = list(shape = 1, rate = 2))
  )
  expect_invisible(
    validate_distribution(dist = "gamma", params = list(shape = 1, scale = 2))
  )
})

test_that("validate_distribution rejects invalid gamma parameterization", {
  expect_error(
    validate_distribution(dist = "gamma", params = list(rate = 1, scale = 2)),
    "Invalid parameter names"
  )
})

test_that("validate_distribution accepts nbinom alternative parameterizations", {
  expect_invisible(
    validate_distribution(dist = "nbinom", params = list(size = 1, prob = 0.5))
  )
  expect_invisible(
    validate_distribution(dist = "nbinom", params = list(size = 1, mu = 3))
  )
})

test_that("validate_distribution accepts multinom special case", {
  expect_invisible(
    validate_distribution(dist = "multinom", params = list(size = 1, prob = c(0.2, 0.8)))
  )
})

test_that("validate_distribution accepts gompertz special case", {
  expect_invisible(
    validate_distribution(dist = "gompertz", params = list(shape = 1, rate = 2))
  )
})

test_that("validate_distribution ignores parameter order", {
  expect_invisible(
    validate_distribution(dist = "norm", params = list(sd = 1, mean = 0))
  )
})

test_that("validate_distribution reports expected parameter names", {
  expect_error(
    validate_distribution(dist = "gamma", params = list(shape = 1)),
    regexp = paste0("Expected:", ".*shape, rate", ".*shape, scale")
  )
})

test_that("validate_distribution returns invisibly TRUE", {
  res <- validate_distribution(dist = "norm", params = list(mean = 0, sd = 1))
  expect_true(res)
})

test_that("validate_distribution returns invisibly", {
  expect_invisible(
    validate_distribution(dist = "norm", params = list(mean = 0, sd = 1))
  )
})


# =============================================================================
# 3. distribution_type()
# =============================================================================

test_that("distribution_type — 'norm' est 'Continuous'", {
  expect_equal(distribution_type("norm"), "Continuous")
})

test_that("distribution_type — 'weibull' est 'Continuous'", {
  expect_equal(distribution_type("weibull"), "Continuous")
})

test_that("distribution_type — 'gamma' est 'Continuous'", {
  expect_equal(distribution_type("gamma"), "Continuous")
})

test_that("distribution_type — 'exp' est 'Continuous'", {
  expect_equal(distribution_type("exp"), "Continuous")
})

test_that("distribution_type — 'pois' est 'Discrete'", {
  expect_equal(distribution_type("pois"), "Discrete")
})

test_that("distribution_type — 'binom' est 'Discrete'", {
  expect_equal(distribution_type("binom"), "Discrete")
})

test_that("distribution_type — 'geom' est 'Discrete'", {
  expect_equal(distribution_type("geom"), "Discrete")
})

test_that("distribution_type — 'nbinom' est 'Discrete'", {
  expect_equal(distribution_type("nbinom"), "Discrete")
})

test_that("distribution_type — 'hyper' est 'Discrete'", {
  expect_equal(distribution_type("hyper"), "Discrete")
})

test_that("distribution_type — 'multinom' est 'Multivariate'", {
  expect_equal(distribution_type("multinom"), "Multivariate")
})

test_that("distribution_type — retourne toujours une chaîne de longueur 1", {
  for (dist in c("norm", "pois", "multinom", "gamma", "binom")) {
    result <- distribution_type(dist)
    expect_length(result, 1L)
    expect_type(result, "character")
  }
})


# =============================================================================
# 4. detect_support_type()
# =============================================================================

# Données continues : peu de doublons (taux < 3 %)
continuous_sample <- function(n = 500, seed = 1) {
  set.seed(seed); stats::rnorm(n)
}

# Données discrètes : beaucoup de doublons (taux >> 3 %)
discrete_sample <- function(n = 500, seed = 2) {
  set.seed(seed); sample(1:5, n, replace = TRUE)
}

test_that("detect_support_type — échantillon continu seul : retourne 'Continuous'", {
  expect_equal(detect_support_type(continuous_sample()), "Continuous")
})

test_that("detect_support_type — échantillon discret seul : retourne 'Discrete'", {
  expect_equal(detect_support_type(discrete_sample()), "Discrete")
})

test_that("detect_support_type — deux échantillons continus : retourne 'Continuous'", {
  expect_equal(detect_support_type(continuous_sample(seed=1), continuous_sample(seed=2)), "Continuous")
})

test_that("detect_support_type — deux échantillons discrets : retourne 'Discrete'", {
  expect_equal(detect_support_type(discrete_sample(seed=1), discrete_sample(seed=2)), "Discrete")
})

test_that("detect_support_type — un discret et un continu : retourne 'Continuous'", {
  # La logique source : si l'un des deux est continu, le résultat est 'Continuous'
  expect_equal(detect_support_type(continuous_sample(), discrete_sample()), "Continuous")
})

test_that("detect_support_type — sample2 = NULL : se comporte comme un seul échantillon", {
  expect_equal(
    detect_support_type(continuous_sample(), NULL),
    detect_support_type(continuous_sample())
  )
})

test_that("detect_support_type — retourne une chaîne de longueur 1", {
  result <- detect_support_type(continuous_sample())
  expect_length(result, 1L)
  expect_type(result, "character")
})

test_that("detect_support_type — valeur dans {'Continuous', 'Discrete'}", {
  r1 <- detect_support_type(continuous_sample())
  r2 <- detect_support_type(discrete_sample())
  expect_true(r1 %in% c("Continuous", "Discrete"))
  expect_true(r2 %in% c("Continuous", "Discrete"))
})

test_that("detect_support_type — ignore les NA dans le calcul du taux de doublons", {
  s <- c(continuous_sample(400), rep(NA, 100))
  expect_equal(detect_support_type(s), "Continuous")
})


# =============================================================================
# 5. is_equal_knownComp()
# =============================================================================

test_that("is_equal_knownComp — deux modèles identiques (norm) : retourne TRUE", {
  m1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  m2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  expect_true(is_equal_knownComp(m1, m2))
})

test_that("is_equal_knownComp — même loi, paramètres différents : retourne FALSE", {
  m1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  m2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 2, sd = 1))
  expect_false(is_equal_knownComp(m1, m2))
})

test_that("is_equal_knownComp — lois différentes : retourne FALSE", {
  m1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  m2 <- admix_model(knownComp_dist = "exp",  knownComp_param = list(rate = 1))
  expect_false(is_equal_knownComp(m1, m2))
})

test_that("is_equal_knownComp — deux modèles exp identiques : retourne TRUE", {
  m1 <- admix_model(knownComp_dist = "exp", knownComp_param = list(rate = 2))
  m2 <- admix_model(knownComp_dist = "exp", knownComp_param = list(rate = 2))
  expect_true(is_equal_knownComp(m1, m2))
})

test_that("is_equal_knownComp — même loi, sd différent : retourne FALSE", {
  m1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  m2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 2))
  expect_false(is_equal_knownComp(m1, m2))
})

test_that("is_equal_knownComp — retourne un scalaire logique", {
  m1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  m2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  result <- is_equal_knownComp(m1, m2)
  expect_type(result, "logical")
  expect_length(result, 1L)
})

test_that("is_equal_knownComp — erreur si admixMod1 n'est pas de classe admix_model", {
  m2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  expect_error(
    is_equal_knownComp(list(comp.dist = "norm"), m2),
    regexp = "admixMod1"
  )
})

test_that("is_equal_knownComp — erreur si admixMod2 n'est pas de classe admix_model", {
  m1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  expect_error(
    is_equal_knownComp(m1, list(comp.dist = "norm")),
    regexp = "admixMod2"
  )
})

test_that("is_equal_knownComp — multinom avec mêmes probabilités : retourne TRUE", {
  m1 <- admix_model(knownComp_dist = "multinom",
                    knownComp_param = list(size = 1, prob = c(0.2, 0.5, 0.3)))
  m2 <- admix_model(knownComp_dist = "multinom",
                    knownComp_param = list(size = 1, prob = c(0.2, 0.5, 0.3)))
  expect_true(is_equal_knownComp(m1, m2))
})

test_that("is_equal_knownComp — multinom avec probabilités différentes : retourne FALSE", {
  m1 <- admix_model(knownComp_dist = "multinom",
                    knownComp_param = list(size = 1, prob = c(0.2, 0.5, 0.3)))
  m2 <- admix_model(knownComp_dist = "multinom",
                    knownComp_param = list(size = 1, prob = c(0.2, 0.4, 0.4)))
  expect_false(is_equal_knownComp(m1, m2))
})


# =============================================================================
# 6. estimVarCov_empProcess()   [fonction interne]
# =============================================================================

test_that("estimVarCov_empProcess — retourne un scalaire numérique", {
  skip_on_cran()
  set.seed(42)
  obs <- stats::rnorm(200)
  val <- admix:::estimVarCov_empProcess(x = 0, y = 0, obs.data = obs)
  expect_length(val, 1L)
  expect_true(is.numeric(val))
})

test_that("estimVarCov_empProcess — valeur toujours dans [0, 0.25] (borne de L*(1-L))", {
  skip_on_cran()
  set.seed(42)
  obs <- stats::rnorm(200)
  for (pt in c(-2, -1, 0, 0.5, 1, 2)) {
    val <- admix:::estimVarCov_empProcess(x = pt, y = pt, obs.data = obs)
    expect_gte(val, 0)
    expect_lte(val, 0.25)
  }
})

test_that("estimVarCov_empProcess — toujours non-négatif (propriété de covariance)", {
  skip_on_cran()
  set.seed(1)
  obs <- stats::rnorm(300)
  for (pair in list(c(-1, 1), c(0, 2), c(-2, 0), c(1, 3))) {
    val <- admix:::estimVarCov_empProcess(x = pair[1], y = pair[2], obs.data = obs)
    expect_gte(val, 0)
  }
})

test_that("estimVarCov_empProcess — symétrie : f(x,y) == f(y,x)", {
  skip_on_cran()
  set.seed(5)
  obs <- stats::rnorm(200)
  v_xy <- admix:::estimVarCov_empProcess(x = -1, y = 1, obs.data = obs)
  v_yx <- admix:::estimVarCov_empProcess(x =  1, y = -1, obs.data = obs)
  expect_equal(v_xy, v_yx)
})

test_that("estimVarCov_empProcess — x = y = point central de la ECDF : valeur cohérente", {
  skip_on_cran()
  # Pour la médiane empirique m, L(m) ≈ 0.5 => L(m)*(1-L(m)) ≈ 0.25 (max)
  set.seed(7)
  obs <- stats::rnorm(1000)
  med <- stats::median(obs)
  val <- admix:::estimVarCov_empProcess(x = med, y = med, obs.data = obs)
  expect_gte(val, 0.20)   # proche de 0.25 pour n grand
  expect_lte(val, 0.25)
})

test_that("estimVarCov_empProcess — points très éloignés des données : covariance quasi-nulle", {
  skip_on_cran()
  set.seed(9)
  obs <- stats::rnorm(500)
  # x très grand => L(x) ≈ 1 => L(min(x,y))*(1-L(max(x,y))) ≈ 0
  val <- admix:::estimVarCov_empProcess(x = 100, y = 200, obs.data = obs)
  expect_lt(val, 1e-6)
})

test_that("estimVarCov_empProcess — formule de Donsker vérifiée manuellement", {
  skip_on_cran()
  set.seed(42)
  obs  <- stats::rnorm(500)
  ecdf_obs <- stats::ecdf(obs)
  x <- 0.3; y <- 1.1
  expected <- ecdf_obs(min(x, y)) * (1 - ecdf_obs(max(x, y)))
  obtained <- admix:::estimVarCov_empProcess(x = x, y = y, obs.data = obs)
  expect_equal(obtained, expected, tolerance = 1e-12)
})


# =============================================================================
# Tests unitaires pour poly_orthonormal_basis()
#
# Comportement observé de la fonction :
#   - "Real"              => liste de (deg+1) polynômes (Hermite He)
#   - "Positive"          => liste de (deg+1) polynômes (Laguerre)
#   - "Bounded.continuous"=> liste de (deg+1) polynômes (Legendre)
#   - "Integer"           => une FONCTION de signature function(x, m, deg)
#                           qui retourne une matrice (length(x), deg)
#   - support invalide    => stop() avec message explicite
# =============================================================================

# =============================================================================
# Helpers
# =============================================================================

# Raccourci pour appeler la fonction interne
poly_basis <- function(...) admix:::poly_orthonormal_basis(...)

# Évalue un polynôme (objet de classe 'polynomial') en un point scalaire
eval_poly <- function(p, x) {
  coefs <- as.numeric(p)           # coefficients croissants : a0, a1, a2, ...
  sum(coefs * x^(seq_along(coefs) - 1))
}


test_that("orthoBasis_coef fonctionne avec le support 'Positive'", {
  set.seed(42)
  sample_pos <- stats::rexp(n = 500, rate = 1)
  coef <- orthoBasis_coef(data = sample_pos, supp = "Positive", degree = 3, m = NULL, bounds = NULL)

  expect_length(coef, 3)
  expect_true(all(sapply(coef, is.numeric)))
})


test_that("orthoBasis_coef fonctionne avec le support 'Bounded.continuous' sans bounds", {
  set.seed(42)
  sample_bounded <- stats::runif(n = 500, min = 2, max = 5)
  coef <- orthoBasis_coef(data = sample_bounded, supp = "Bounded.continuous", degree = 3, m = NULL, bounds = NULL)

  expect_length(coef, 3)
  expect_true(all(sapply(coef, is.numeric)))
})


test_that("orthoBasis_coef fonctionne avec le support 'Bounded.continuous' avec bounds explicites", {
  set.seed(42)
  sample_bounded <- stats::runif(n = 500, min = 2, max = 5)
  # other : list( list(min bounds) , list(max bounds) ) — le 2e élément donne les bornes sup
  bound <- c(0,6)
  coef <- orthoBasis_coef(data = sample_bounded, supp = "Bounded.continuous", degree = 3, bounds = bound)

  expect_length(coef, 3)
  expect_true(all(sapply(coef, is.numeric)))
})


test_that("orthoBasis_coef stoppe avec le support 'Integer'", {
  expect_error(
    orthoBasis_coef(data = rpois(100, lambda = 3), supp = "Integer", degree = 3, m = 3, bounds = NULL),
    "Still to implement"
  )
})


test_that("orthoBasis_coef stoppe avec un support invalide", {
  expect_error(
    orthoBasis_coef(data = rnorm(100), supp = "Bounded.discrete", degree = 3, m = NULL, bounds = NULL),
    "support"
  )
})


# =============================================================================
# 1. Support "Real" — polynômes de Hermite (He)
# =============================================================================

test_that("poly_orthonormal_basis Real — retourne une liste", {
  res <- poly_basis("Real", deg = 3, x = NULL, m = NULL)
  expect_type(res, "list")
})

test_that("poly_orthonormal_basis stoppe avec un support invalide", {
  expect_error(
    poly_orthonormal_basis(support = "Bounded.discrete", deg = 3, x = NULL, m = NULL),
    "correct argument"
  )
})

test_that("poly_orthonormal_basis Real — longueur de la liste = deg + 1", {
  for (d in c(1, 3, 5)) {
    res <- poly_basis("Real", deg = d, x = NULL, m = NULL)
    expect_length(res, d + 1L)
  }
})

test_that("poly_orthonormal_basis Real — premier élément est le polynôme constant 1", {
  res <- poly_basis("Real", deg = 3, x = NULL, m = NULL)
  # He_0(x) = 1 : évaluation en plusieurs points doit donner 1
  for (pt in c(-2, 0, 1, 3)) {
    expect_equal(eval_poly(res[[1]], pt), 1, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Real — deuxième élément est He_1(x) = x", {
  res <- poly_basis("Real", deg = 3, x = NULL, m = NULL)
  for (pt in c(-2, 0, 1, 3)) {
    expect_equal(eval_poly(res[[2]], pt), pt, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Real — troisième élément est He_2(x) = x^2 - 1", {
  res <- poly_basis("Real", deg = 3, x = NULL, m = NULL)
  for (pt in c(-2, 0, 1, 3)) {
    expect_equal(eval_poly(res[[3]], pt), pt^2 - 1, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Real — quatrième élément est He_3(x) = x^3 - 3x", {
  res <- poly_basis("Real", deg = 3, x = NULL, m = NULL)
  for (pt in c(-2, 0, 1, 3)) {
    expect_equal(eval_poly(res[[4]], pt), pt^3 - 3 * pt, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Real — chaque élément est un objet polynomial", {
  res <- poly_basis("Real", deg = 4, x = NULL, m = NULL)
  for (i in seq_along(res)) {
    expect_true(inherits(res[[i]], "polynomial"))
  }
})

test_that("poly_orthonormal_basis Real — deg = 1 retourne une liste de longueur 2", {
  res <- poly_basis("Real", deg = 1, x = NULL, m = NULL)
  expect_length(res, 2L)
})


# =============================================================================
# 2. Support "Positive" — polynômes de Laguerre
# =============================================================================

test_that("poly_orthonormal_basis Positive — retourne une liste", {
  res <- poly_basis("Positive", deg = 3, x = NULL, m = NULL)
  expect_type(res, "list")
})

test_that("poly_orthonormal_basis Positive — longueur = deg + 1", {
  for (d in c(1, 3, 5)) {
    res <- poly_basis("Positive", deg = d, x = NULL, m = NULL)
    expect_length(res, d + 1L)
  }
})

test_that("poly_orthonormal_basis Positive — premier élément est L_0(x) = 1", {
  res <- poly_basis("Positive", deg = 3, x = NULL, m = NULL)
  for (pt in c(0, 1, 2, 5)) {
    expect_equal(eval_poly(res[[1]], pt), 1, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Positive — deuxième élément est L_1(x) = 1 - x", {
  res <- poly_basis("Positive", deg = 3, x = NULL, m = NULL)
  for (pt in c(0, 1, 2, 5)) {
    expect_equal(eval_poly(res[[2]], pt), 1 - pt, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Positive — troisième élément est L_2(x) = 1 - 2x + x²/2", {
  res <- poly_basis("Positive", deg = 3, x = NULL, m = NULL)
  for (pt in c(0, 1, 2, 5)) {
    expect_equal(eval_poly(res[[3]], pt), 1 - 2*pt + 0.5*pt^2, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Positive — chaque élément est un objet polynomial", {
  res <- poly_basis("Positive", deg = 4, x = NULL, m = NULL)
  for (i in seq_along(res)) {
    expect_true(inherits(res[[i]], "polynomial"))
  }
})


# =============================================================================
# 3. Support "Bounded.continuous" — polynômes de Legendre
# =============================================================================

test_that("poly_orthonormal_basis Bounded.continuous — retourne une liste", {
  res <- poly_basis("Bounded.continuous", deg = 3, x = NULL, m = NULL)
  expect_type(res, "list")
})

test_that("poly_orthonormal_basis Bounded.continuous — longueur = deg + 1", {
  for (d in c(1, 3, 5)) {
    res <- poly_basis("Bounded.continuous", deg = d, x = NULL, m = NULL)
    expect_length(res, d + 1L)
  }
})

test_that("poly_orthonormal_basis Bounded.continuous — premier élément est P_0(x) = 1", {
  res <- poly_basis("Bounded.continuous", deg = 3, x = NULL, m = NULL)
  for (pt in c(-1, 0, 0.5, 1)) {
    expect_equal(eval_poly(res[[1]], pt), 1, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Bounded.continuous — deuxième élément est P_1(x) = x", {
  res <- poly_basis("Bounded.continuous", deg = 3, x = NULL, m = NULL)
  for (pt in c(-1, 0, 0.5, 1)) {
    expect_equal(eval_poly(res[[2]], pt), pt, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Bounded.continuous — troisième élément est P_2(x) = -0.5 + 1.5x²", {
  res <- poly_basis("Bounded.continuous", deg = 3, x = NULL, m = NULL)
  for (pt in c(-1, 0, 0.5, 1)) {
    expect_equal(eval_poly(res[[3]], pt), -0.5 + 1.5*pt^2, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Bounded.continuous — quatrième élément est P_3(x) = -1.5x + 2.5x³", {
  res <- poly_basis("Bounded.continuous", deg = 3, x = NULL, m = NULL)
  for (pt in c(-1, 0, 0.5, 1)) {
    expect_equal(eval_poly(res[[4]], pt), -1.5*pt + 2.5*pt^3, tolerance = 1e-10)
  }
})

test_that("poly_orthonormal_basis Bounded.continuous — chaque élément est un objet polynomial", {
  res <- poly_basis("Bounded.continuous", deg = 4, x = NULL, m = NULL)
  for (i in seq_along(res)) {
    expect_true(inherits(res[[i]], "polynomial"))
  }
})


# =============================================================================
# 4. Support "Integer" — polynômes de Charlier (retourne une fonction)
# =============================================================================

test_that("poly_orthonormal_basis Integer — retourne une fonction (closure)", {
  res <- poly_basis("Integer", deg = 3, x = 2, m = 3)
  expect_true(is.function(res))
  expect_equal(typeof(res), "closure")
})

test_that("poly_orthonormal_basis Integer — la fonction retourne une matrice", {
  res <- poly_basis("Integer", deg = 3, x = 2, m = 3)
  mat <- res(x = c(1, 2, 3), m = 3, deg = 3)
  expect_true(is.matrix(mat))
})

test_that("poly_orthonormal_basis Integer — dimensions de la matrice : (length(x), deg)", {
  res <- poly_basis("Integer", deg = 3, x = 2, m = 3)
  x_vals <- c(1, 2, 3)
  mat <- res(x = x_vals, m = 3, deg = 3)
  expect_equal(dim(mat), c(length(x_vals), 3L))
})

test_that("poly_orthonormal_basis Integer — valeurs numériques finies", {
  res <- poly_basis("Integer", deg = 3, x = 2, m = 3)
  mat <- res(x = c(1, 2, 3), m = 3, deg = 3)
  expect_true(all(is.finite(mat)))
})

test_that("poly_orthonormal_basis Integer — valeurs reproduisent la sortie observée", {
  res  <- poly_basis("Integer", deg = 3, x = 2, m = 3)
  mat  <- res(x = c(1, 2, 3), m = 3, deg = 3)
  # Valeurs exactes issues de la sortie de débogage fournie
  expected <- matrix(
    c( 1.1547005,  0.7071068,  0.0000000,
       0.5773503, -0.2357023, -0.7071068,
       0.0000000, -0.7071068, -0.4714045),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  expect_equal(mat, expected, tolerance = 1e-6)
})

test_that("poly_orthonormal_basis Integer — fonctionne pour un seul point x", {
  res <- poly_basis("Integer", deg = 3, x = 2, m = 3)
  mat <- res(x = 2, m = 3, deg = 3)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(1L, 3L))
})

test_that("poly_orthonormal_basis Integer — dimensions correctes avec deg = 5", {
  res  <- poly_basis("Integer", deg = 5, x = 2, m = 3)
  mat  <- res(x = 1:4, m = 3, deg = 5)
  expect_equal(dim(mat), c(4L, 5L))
})


# =============================================================================
# 5. Support invalide
# =============================================================================

test_that("poly_orthonormal_basis — support invalide lève une erreur explicite", {
  expect_error(
    poly_basis("Invalid", deg = 3, x = NULL, m = NULL),
    regexp = "Please give a correct argument for the support"
  )
})

test_that("poly_orthonormal_basis — support vide lève une erreur", {
  expect_error(
    poly_basis("", deg = 3, x = NULL, m = NULL)
  )
})


# =============================================================================
# 6. Propriétés transversales (indépendantes du support)
# =============================================================================

test_that("poly_orthonormal_basis — deg = 0 retourne une liste de longueur 1 pour Real", {
  res <- poly_basis("Real", deg = 0, x = NULL, m = NULL)
  expect_length(res, 1L)
  expect_equal(eval_poly(res[[1]], 5), 1, tolerance = 1e-10)
})

test_that("poly_orthonormal_basis — résultats identiques à deux appels successifs (déterminisme)", {
  res1 <- poly_basis("Real",               deg = 4, x = NULL, m = NULL)
  res2 <- poly_basis("Real",               deg = 4, x = NULL, m = NULL)
  res3 <- poly_basis("Positive",           deg = 4, x = NULL, m = NULL)
  res4 <- poly_basis("Positive",           deg = 4, x = NULL, m = NULL)
  res5 <- poly_basis("Bounded.continuous", deg = 4, x = NULL, m = NULL)
  res6 <- poly_basis("Bounded.continuous", deg = 4, x = NULL, m = NULL)

  for (i in seq_along(res1)) {
    expect_equal(as.numeric(res1[[i]]), as.numeric(res2[[i]]))
    expect_equal(as.numeric(res3[[i]]), as.numeric(res4[[i]]))
    expect_equal(as.numeric(res5[[i]]), as.numeric(res6[[i]]))
  }
})

test_that("poly_orthonormal_basis — les listes pour Real/Positive/Bounded.continuous sont distinctes", {
  r <- poly_basis("Real",               deg = 3, x = NULL, m = NULL)
  p <- poly_basis("Positive",           deg = 3, x = NULL, m = NULL)
  b <- poly_basis("Bounded.continuous", deg = 3, x = NULL, m = NULL)
  # Le polynôme de degré 2 (indice 3) diffère entre les trois bases
  expect_false(isTRUE(all.equal(as.numeric(r[[3]]), as.numeric(p[[3]]))))
  expect_false(isTRUE(all.equal(as.numeric(r[[3]]), as.numeric(b[[3]]))))
  expect_false(isTRUE(all.equal(as.numeric(p[[3]]), as.numeric(b[[3]]))))
})


# =============================================================================
#  knownComp_to_uniform()   [fonction interne]
# =============================================================================

test_that("knownComp_to_uniform() rejects invalid admixMod argument", {
  data <- rnorm(100)
  expect_error(knownComp_to_uniform(data = data, admixMod = list()), regexp = "not correctly specified")
  expect_error(knownComp_to_uniform(data = data, admixMod = "wrong"), regexp = "not correctly specified")
})


# ── Continuous support ────────────────────────────────────────────────────────

test_that("knownComp_to_uniform() returns values in [0,1] for continuous data (Normal known component)", {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.4, comp.dist = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  admixMod1 <- admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]])
  result <- knownComp_to_uniform(data = data1, admixMod = admixMod1)

  expect_true(all(result >= 0 & result <= 1))
})

test_that("knownComp_to_uniform() returns a numeric vector of the same length as input (continuous)", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  admixMod1 <- admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]])
  result <- knownComp_to_uniform(data = data1, admixMod = admixMod1)

  expect_type(result, "double")
  expect_length(result, length(data1))
})

test_that("knownComp_to_uniform() output is approximately Uniform(0,1) when data is purely from the known component (continuous)", {
  # If all data come from g (the known component), the CDF transform F_g(X) ~ Uniform(0,1).
  set.seed(123)
  pure_g    <- rnorm(2000, mean = 0, sd = 1)
  admixMod1 <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  result <- knownComp_to_uniform(data = pure_g, admixMod = admixMod1)

  ks_test <- ks.test(result, "punif", 0, 1)
  expect_gt(ks_test$p.value, 0.05)
})

test_that("knownComp_to_uniform() works with an Exponential known component (continuous)", {
  set.seed(7)
  mixt2 <- twoComp_mixt(n = 500, weight = 0.3, comp.dist  = list("norm", "exp"),
                        comp.param = list(list(mean = 5, sd = 1), list(rate = 1)))
  data2 <- get_mixture_data(mixt2)
  admixMod2 <- admix_model(knownComp_dist  = mixt2$comp.dist[[2]], knownComp_param = mixt2$comp.param[[2]])
  result <- knownComp_to_uniform(data = data2, admixMod = admixMod2)

  expect_true(all(result >= 0 & result <= 1))
  expect_length(result, length(data2))
})

test_that("knownComp_to_uniform() is deterministic for continuous data (no randomness)", {
  set.seed(99)
  mixt1 <- twoComp_mixt(n = 200, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  admixMod1 <- admix_model(knownComp_dist  = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]])
  result1 <- knownComp_to_uniform(data = data1, admixMod = admixMod1)
  result2 <- knownComp_to_uniform(data = data1, admixMod = admixMod1)

  expect_equal(result1, result2)
})


# ── Discrete support ──────────────────────────────────────────────────────────

test_that("knownComp_to_uniform() returns values in [0,1] for discrete data (Poisson known component)", {
  set.seed(10)
  mixt3 <- twoComp_mixt(n = 500, weight = 0.4, comp.dist  = list("pois", "pois"),
                        comp.param = list(list(lambda = 5), list(lambda = 2)))
  data3     <- get_mixture_data(mixt3)
  admixMod3 <- admix_model(knownComp_dist  = mixt3$comp.dist[[2]], knownComp_param = mixt3$comp.param[[2]])
  result <- knownComp_to_uniform(data = data3, admixMod = admixMod3)

  expect_true(all(result >= 0 & result <= 1))
})

test_that("knownComp_to_uniform() returns a numeric vector of the same length as input (discrete)", {
  set.seed(11)
  mixt3 <- twoComp_mixt(n = 300, weight = 0.5, comp.dist  = list("pois", "pois"),
                        comp.param = list(list(lambda = 5), list(lambda = 2)))
  data3     <- get_mixture_data(mixt3)
  admixMod3 <- admix_model(knownComp_dist  = mixt3$comp.dist[[2]], knownComp_param = mixt3$comp.param[[2]])
  result <- knownComp_to_uniform(data = data3, admixMod = admixMod3)

  expect_type(result, "double")
  expect_length(result, length(data3))
})

test_that("knownComp_to_uniform() introduces randomness for discrete data (runif jittering)", {
  # For discrete data the function uses runif() internally, so two calls on the
  # same data with different seeds must yield different results.
  set.seed(20)
  mixt3 <- twoComp_mixt(n = 300, weight = 0.5, comp.dist  = list("pois", "pois"),
                        comp.param = list(list(lambda = 5), list(lambda = 2)))
  data3     <- get_mixture_data(mixt3)
  admixMod3 <- admix_model(knownComp_dist  = mixt3$comp.dist[[2]], knownComp_param = mixt3$comp.param[[2]])
  set.seed(1); result1 <- knownComp_to_uniform(data = data3, admixMod = admixMod3)
  set.seed(2); result2 <- knownComp_to_uniform(data = data3, admixMod = admixMod3)

  expect_false(identical(result1, result2))
})

test_that("knownComp_to_uniform() PMF normalization holds for discrete data (Binomial known component)", {
  # Verifies the function runs without error and the output stays in [0,1]
  # for a Binomial known component, which exercises the pmf / sum(pmf) branch.
  set.seed(30)
  mixt4 <- twoComp_mixt(n = 400, weight = 0.35, comp.dist  = list("pois", "binom"),
                        comp.param = list(list(lambda = 5), list(size = 10, prob = 0.4)))
  data4     <- get_mixture_data(mixt4)
  admixMod4 <- admix_model(knownComp_dist  = mixt4$comp.dist[[2]], knownComp_param = mixt4$comp.param[[2]])

  expect_no_error({result <- knownComp_to_uniform(data = data4, admixMod = admixMod4)})

  result <- knownComp_to_uniform(data = data4, admixMod = admixMod4)
  expect_true(all(result >= 0 & result <= 1))
})


test_that("knownComp_to_uniform couvre la branche dmultinom (support discret)", {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.6, comp.dist  = list("multinom", "multinom"),
    comp.param = list(list(size = 1, prob = c(0.2, 0.5, 0.3)), list(size = 1, prob = c(0.1, 0.6, 0.3))))
  data1 <- get_mixture_data(mixt1)
  mod1 <- admix_model(knownComp_dist  = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))
  result <- knownComp_to_uniform(data = data1, admixMod = mod1)

  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
  expect_length(result, length(data1))
})
