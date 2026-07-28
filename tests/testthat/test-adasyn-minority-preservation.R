test_that("sby_adasyn_ratio abaixo de um sempre aumenta a classe rara", {
  dat <- data.frame(
    x1 = c(-4, -3, -2, -1, seq_len(10)),
    x2 = c(1, 2, 1, 2, seq_len(10) + 10),
    y = factor(c(rep("rare", 4), rep("common", 10)),
               levels = c("rare", "common"))
  )
  x <- data.matrix(dat[c("x1", "x2")])

  matrix_out <- sby_adasyn_matrix(
    x, dat$y, sby_adasyn_ratio = 0.4, sby_seed = 17L,
    sby_knn_engine = "FNN"
  )
  tabular_out <- sby_adasyn(
    y ~ ., dat, sby_adasyn_ratio = 0.4, sby_seed = 17L,
    sby_knn_engine = "FNN"
  )
  balance_out <- sby_balance_matrix(
    x, dat$y, sby_strategy = "adasyn", sby_adasyn_ratio = 0.4,
    sby_seed = 17L, sby_knn_engine = "FNN"
  )

  # O ratio 0.4 vira 1.4: floor(4 * 1.4) = 5 raras no resultado.
  expect_equal(sum(matrix_out$sby_y_vector == "rare"), 5L)
  expect_equal(sum(tabular_out$y == "rare"), 5L)
  expect_equal(sum(balance_out$sby_y_vector == "rare"), 5L)
  expect_gte(sum(matrix_out$sby_y_vector == "rare"), sum(dat$y == "rare"))
})

test_that("a quantidade final ADASYN sempre usa um mais o ratio", {
  target <- factor(c(rep("rare", 5), rep("common", 12)),
                   levels = c("rare", "common"))

  generated <- sby_compute_minority_expansion_count(target, 0.4)
  expected_final <- floor(5 * (1 + 0.4))

  expect_equal(5L + generated, expected_final)
  expect_equal(generated, 2L)
})

test_that("ratios nulos produzem inercia e negativos sao rejeitados", {
  target <- factor(c(rep("rare", 5), rep("common", 12)))

  expect_equal(sby_compute_minority_expansion_count(target, 0), 0L)
  expect_equal(sby_compute_majority_retention_count(target, 0), 12L)
  expect_error(sby_compute_minority_expansion_count(target, -1), "negativo")
  expect_error(sby_compute_majority_retention_count(target, -0.5), "negativo")
})

test_that("ADASYN e ADANEAR rejeitam qualquer reducao da classe rara", {
  input <- factor(c(rep("rare", 4), rep("common", 10)))
  preserved <- factor(c(rep("rare", 5), rep("common", 6)))
  reduced <- factor(c(rep("rare", 3), rep("common", 6)))

  expect_invisible(sby_assert_minority_not_reduced(input, preserved))
  expect_error(
    sby_assert_minority_not_reduced(input, reduced, "teste"),
    regexp = "classe rara original nunca pode ser reduzida"
  )
})

test_that("ADANEAR preserva a classe rara original com ratio 0.4", {
  dat <- data.frame(
    x1 = c(-4, -3, -2, -1, seq_len(10)),
    x2 = c(1, 2, 1, 2, seq_len(10) + 10),
    y = factor(c(rep("rare", 4), rep("common", 10)),
               levels = c("rare", "common"))
  )
  x <- data.matrix(dat[c("x1", "x2")])

  matrix_out <- sby_adanear_matrix(
    x, dat$y, sby_adasyn_ratio = 0.4, sby_nearmiss_ratio = 0.5,
    sby_seed = 23L, sby_knn_engine = "FNN"
  )
  tabular_out <- sby_adanear(
    y ~ ., dat, sby_adasyn_ratio = 0.4, sby_nearmiss_ratio = 0.5,
    sby_seed = 23L, sby_knn_engine = "FNN"
  )

  expect_equal(sum(matrix_out$sby_y_vector == "rare"), 5L)
  expect_equal(sum(tabular_out$y == "rare"), 5L)
  expect_gte(sum(matrix_out$sby_y_vector == "rare"), sum(dat$y == "rare"))
})

test_that("NearMiss preserva o papel raro original apos ADASYN inverter contagens", {
  dat <- data.frame(
    x1 = c(-5, -4, -3, seq_len(7)),
    x2 = c(1, 2, 3, seq_len(7) + 10),
    y = factor(c(rep("rare", 3), rep("common", 7)),
               levels = c("rare", "common"))
  )
  x <- data.matrix(dat[c("x1", "x2")])

  out <- sby_adanear_matrix(
    x, dat$y,
    sby_adasyn_ratio = 2,
    sby_nearmiss_ratio = 0.5,
    sby_seed = 29L,
    sby_knn_engine = "FNN"
  )

  # ADASYN: floor(3 * (1 + 2)) = 9 raras. NearMiss remove apenas a maioria.
  expect_equal(sum(out$sby_y_vector == "rare"), 9L)
  expect_equal(sum(out$sby_y_vector == "common"), 4L)
  expect_gte(sum(out$sby_y_vector == "rare"), sum(dat$y == "rare"))
})

test_that("os papeis de classe expoem codigos de nivel que indexam o factor", {
  # Regressao: as rotas HPC selecionam as linhas raras por codigo de nivel.
  # Quando o campo nao existia, as.integer(NULL) produzia integer(0) e o
  # which() devolvia zero linhas, descartando toda a classe rara original.
  y <- factor(c(rep("rare", 4), rep("common", 10)),
              levels = c("rare", "common"))
  roles <- sby_binary_class_counts_fast(y)

  expect_type(roles$sby_minority_level, "integer")
  expect_type(roles$sby_majority_level, "integer")
  expect_length(roles$sby_minority_level, 1L)
  expect_length(roles$sby_majority_level, 1L)
  expect_identical(
    levels(y)[roles$sby_minority_level], roles$sby_minority_label
  )
  expect_identical(
    levels(y)[roles$sby_majority_level], roles$sby_majority_label
  )
  expect_equal(
    sum(as.integer(y) == roles$sby_minority_level), roles$sby_minority_count
  )
  expect_equal(
    sum(as.integer(y) == roles$sby_majority_level), roles$sby_majority_count
  )

  # Tambem vale quando a classe rara nao e o primeiro nivel do factor.
  y_inv <- factor(c(rep("common", 10), rep("rare", 4)),
                  levels = c("common", "rare"))
  roles_inv <- sby_binary_class_counts_fast(y_inv)
  expect_identical(roles_inv$sby_minority_label, "rare")
  expect_equal(sum(as.integer(y_inv) == roles_inv$sby_minority_level), 4L)
})

test_that("as rotas HPC tratam sby_adasyn_ratio como acrescimo sobre a rara", {
  skip_if_not(sby_adanear_hpc_available())
  set.seed(2024)
  dat <- data.frame(
    id = seq_len(50) * 1.0,
    x1 = c(rnorm(10, -3), rnorm(40, 3)),
    y = factor(c(rep("rare", 10), rep("common", 40)),
               levels = c("rare", "common"))
  )
  rare_ids <- dat$id[dat$y == "rare"]

  # ratio 0.4 vira 1.4: floor(10 * 1.4) = 14 raras, ou seja 4 sinteticas.
  adanear_out <- sby_adanear_hpc(
    dat, y ~ ., sby_adasyn_k = 3, sby_nearmiss_k = 3,
    sby_adasyn_ratio = 0.4, sby_nearmiss_ratio = 1, sby_seed = 42L
  )
  adasyn_out <- sby_adasyn_hpc(
    dat, y ~ ., sby_adasyn_k = 3, sby_adasyn_ratio = 0.4, sby_seed = 42L
  )

  expect_equal(sum(adanear_out$y == "rare"), floor(10 * 1.4))
  expect_equal(sum(adasyn_out$y == "rare"), floor(10 * 1.4))

  # Nenhuma linha rara original pode desaparecer da saida.
  expect_true(all(rare_ids %in% adanear_out$id[adanear_out$y == "rare"]))
  expect_true(all(rare_ids %in% adasyn_out$id[adasyn_out$y == "rare"]))

  # NearMiss atua apenas na maioria: floor(14 * 1) comuns retidos.
  expect_equal(sum(adanear_out$y == "common"), 14L)
})

test_that("sby_nearmiss_hpc reduz somente a maioria e preserva toda a rara", {
  skip_if_not(sby_adanear_hpc_available())
  set.seed(2025)
  dat <- data.frame(
    id = seq_len(50) * 1.0,
    x1 = c(rnorm(10, -3), rnorm(40, 3)),
    y = factor(c(rep("rare", 10), rep("common", 40)),
               levels = c("rare", "common"))
  )
  rare_ids <- dat$id[dat$y == "rare"]

  out <- sby_nearmiss_hpc(
    dat, y ~ ., sby_nearmiss_k = 3, sby_nearmiss_ratio = 0.5, sby_seed = 42L
  )

  expect_equal(sum(out$y == "rare"), 10L)
  expect_true(all(rare_ids %in% out$id[out$y == "rare"]))
  expect_equal(sum(out$y == "common"), floor(10 * 0.5))
})

test_that("sby_adasyn_ratio igual a zero mantem a classe rara inalterada no HPC", {
  skip_if_not(sby_adanear_hpc_available())
  set.seed(2026)
  dat <- data.frame(
    id = seq_len(50) * 1.0,
    x1 = c(rnorm(10, -3), rnorm(40, 3)),
    y = factor(c(rep("rare", 10), rep("common", 40)),
               levels = c("rare", "common"))
  )
  rare_ids <- dat$id[dat$y == "rare"]

  out <- sby_adanear_hpc(
    dat, y ~ ., sby_adasyn_k = 3, sby_nearmiss_k = 3,
    sby_adasyn_ratio = 0, sby_nearmiss_ratio = 1, sby_seed = 42L
  )

  expect_equal(sum(out$y == "rare"), 10L)
  expect_true(all(rare_ids %in% out$id[out$y == "rare"]))
  expect_equal(sum(out$y == "common"), 10L)
})

test_that("NearMiss nao troca nem ignora papeis quando ADASYN iguala as classes", {
  dat <- data.frame(
    x1 = c(-4, -3, -2, -1, seq_len(6)),
    x2 = c(1, 2, 1, 2, seq_len(6) + 10),
    y = factor(c(rep("rare", 4), rep("common", 6)),
               levels = c("rare", "common"))
  )
  x <- data.matrix(dat[c("x1", "x2")])

  out <- sby_adanear_matrix(
    x, dat$y,
    sby_adasyn_ratio = 0.5,
    sby_nearmiss_ratio = 0.5,
    sby_seed = 31L,
    sby_knn_engine = "FNN"
  )

  # ADASYN leva a rara a 6; NearMiss ainda deve reter floor(6 * 0.5) = 3 comuns.
  expect_equal(sum(out$sby_y_vector == "rare"), 6L)
  expect_equal(sum(out$sby_y_vector == "common"), 3L)
})
