test_that("ADASYN HPC aloca sinteticas por densidade majoritaria, nao round-robin", {
  skip_if_not(sby_adanear_hpc_available())

  # Tres grupos raros bem separados no espaco, com regimes de densidade
  # majoritaria distintos. Como os vizinhos minoritarios de cada ponto ficam
  # sempre dentro do proprio grupo, toda sintetica cai no envelope do grupo que
  # a originou, o que torna a atribuicao por proximidade exata.
  set.seed(11)
  g <- 8L
  centers <- c(0, 15, 45)
  rare <- rbind(
    matrix(stats::rnorm(g * 2L, centers[1], 1.5), ncol = 2L),  # cercado de maioria
    matrix(stats::rnorm(g * 2L, centers[2], 1.5), ncol = 2L),  # pouca maioria em volta
    matrix(stats::rnorm(g * 2L, centers[3], 1.5), ncol = 2L)   # isolado
  )
  maj <- rbind(
    matrix(stats::rnorm(60L * 2L, centers[1], 1.5), ncol = 2L),
    matrix(stats::rnorm(4L * 2L, centers[2], 1.5), ncol = 2L)
  )

  dat <- data.frame(
    x1 = c(rare[, 1L], maj[, 1L]),
    x2 = c(rare[, 2L], maj[, 2L]),
    y  = factor(c(rep("rare", 3L * g), rep("common", 64L)),
                levels = c("rare", "common"))
  )

  out <- sby_adasyn_hpc(dat, y ~ ., sby_adasyn_k = 5, sby_adasyn_ratio = 1,
                        sby_seed = 42L)

  rare_out <- out[out$y == "rare", , drop = FALSE]
  expect_equal(nrow(rare_out), 2L * 3L * g)

  # Apenas as sinteticas: as originais aparecem primeiro, na ordem de entrada.
  syn <- utils::tail(rare_out, 3L * g)
  syn_center <- rowMeans(as.matrix(syn[, c("x1", "x2")]))
  group_of <- max.col(-abs(outer(syn_center, centers, "-")))
  per_group <- tabulate(group_of, nbins = 3L)

  expect_equal(sum(per_group), 3L * g)
  # Round-robin distribuiria igualmente entre os tres grupos.
  expect_gt(length(unique(per_group)), 1L)
  expect_gt(per_group[[1L]], per_group[[2L]])
  # Pontos raros sem nenhum vizinho majoritario tem r_i = 0 e nao geram nada.
  expect_identical(per_group[[3L]], 0L)
})

test_that("atalhos HPC aceitam formulas com subconjunto de preditores", {
  skip_if_not(sby_adanear_hpc_available())

  set.seed(202)
  dat <- data.frame(
    x1 = c(stats::rnorm(6, -2), stats::rnorm(18, 2)),
    x2 = c(stats::rnorm(6, -2), stats::rnorm(18, 2)),
    x3 = c(stats::rnorm(6, -2), stats::rnorm(18, 2)),
    y  = factor(c(rep("min", 6), rep("maj", 18)), levels = c("min", "maj"))
  )

  ada <- sby_adasyn_hpc(dat, y ~ x1 + x3, sby_adasyn_k = 3,
                        sby_adasyn_ratio = 0.5, sby_seed = 7L)
  expect_identical(names(ada), c("x1", "x3", "y"))
  expect_gte(sum(ada$y == "min"), sum(dat$y == "min"))

  near <- sby_nearmiss_hpc(dat, y ~ x2 + x3, sby_nearmiss_k = 3,
                           sby_nearmiss_ratio = 1, sby_seed = 7L)
  expect_identical(names(near), c("x2", "x3", "y"))
  expect_identical(sum(near$y == "min"), sum(dat$y == "min"))

  both <- sby_adanear_hpc(dat, y ~ x1 + x2, sby_adasyn_k = 3,
                          sby_nearmiss_k = 3, sby_adasyn_ratio = 0.5,
                          sby_nearmiss_ratio = 1, sby_seed = 7L)
  expect_identical(names(both), c("x1", "x2", "y"))
  expect_gte(sum(both$y == "min"), sum(dat$y == "min"))
})

test_that("validacao das razoes rejeita valores nao finitos", {
  dat <- data.frame(
    x1 = seq_len(14),
    x2 = seq_len(14) / 2,
    y  = factor(c(rep("min", 4), rep("maj", 10)), levels = c("min", "maj"))
  )

  expect_error(sby_adasyn_hpc(dat, y ~ ., sby_adasyn_ratio = Inf),
               regexp = "sby_adasyn_ratio")
  expect_error(sby_nearmiss_hpc(dat, y ~ ., sby_nearmiss_ratio = Inf),
               regexp = "sby_nearmiss_ratio")
  expect_error(sby_adanear_hpc(dat, y ~ ., sby_nearmiss_ratio = Inf),
               regexp = "sby_nearmiss_ratio")
  expect_error(sby_adanear_hpc(dat, y ~ ., sby_adasyn_ratio = Inf),
               regexp = "sby_adasyn_ratio")
})

test_that("classes perfeitamente balanceadas abortam nas duas rotas de papeis", {
  balanced <- factor(rep(c("a", "b"), each = 6L), levels = c("a", "b"))

  expect_error(
    sbyadanear:::sby_binary_class_counts_fast(balanced),
    regexp = "desbalanceadas"
  )
  expect_error(
    sbyadanear:::sby_get_binary_class_roles(balanced),
    regexp = "desbalanceadas"
  )
})

test_that("sby_config_max_threads nao vaza para o estado global do OpenMP", {
  skip_if_not(sby_adanear_hpc_available())

  set.seed(808)
  dat <- data.frame(
    x1 = c(stats::rnorm(6, -2), stats::rnorm(18, 2)),
    x2 = c(stats::rnorm(6, -2), stats::rnorm(18, 2)),
    y  = factor(c(rep("min", 6), rep("maj", 18)), levels = c("min", "maj"))
  )

  read_max_threads <- function(){
    report <- sbyadanear:::sby_call_native("sby_hpc_compile_report_cpp")
    as.integer(report$openmp_max_threads)
  }

  before <- read_max_threads()
  invisible(sby_adanear_hpc(dat, y ~ ., sby_adasyn_k = 3, sby_nearmiss_k = 3,
                            sby_adasyn_ratio = 0.5, sby_nearmiss_ratio = 1,
                            sby_config_max_threads = 1L, sby_seed = 5L))
  expect_identical(read_max_threads(), before)

  invisible(sby_adasyn_hpc(dat, y ~ ., sby_adasyn_k = 3,
                           sby_adasyn_ratio = 0.5,
                           sby_config_max_threads = 1L, sby_seed = 5L))
  expect_identical(read_max_threads(), before)
})

test_that("relatorio de compilacao reporta apenas capacidades reais", {
  skip_if_not(sby_adanear_hpc_available())

  report <- sbyadanear:::sby_call_native("sby_hpc_compile_report_cpp")

  expect_true("mkl_linked" %in% names(report))
  expect_true(is.logical(report$mkl_linked))
  # cascade_lake_native e derivado das macros AVX-512 reais do compilador.
  expect_identical(
    isTRUE(report$cascade_lake_native),
    isTRUE(report$avx512f) && isTRUE(report$avx512cd) &&
      isTRUE(report$avx512bw) && isTRUE(report$avx512dq) &&
      isTRUE(report$avx512vl)
  )
})

test_that("razoes que nao rendem uma linha inteira nao inventam registros", {
  skip_if_not(sby_adanear_hpc_available())

  set.seed(31)
  dat <- data.frame(
    x1 = c(stats::rnorm(4, -2), stats::rnorm(40, 2)),
    x2 = c(stats::rnorm(4, -2), stats::rnorm(40, 2)),
    y  = factor(c(rep("min", 4), rep("maj", 40)), levels = c("min", "maj"))
  )

  # floor(4 * 1.1) - 4 = 0 sinteticas: a rota HPC nao arredonda para 1.
  out <- sby_adasyn_hpc(dat, y ~ ., sby_adasyn_k = 3, sby_adasyn_ratio = 0.1,
                        sby_seed = 3L)
  expect_identical(nrow(out), nrow(dat))
  expect_identical(sum(out$y == "min"), sum(dat$y == "min"))

  # floor(4 * 0.1) = 0 majoritarios retidos: aborta em vez de reter uma linha.
  expect_error(
    sby_nearmiss_hpc(dat, y ~ ., sby_nearmiss_k = 3,
                     sby_nearmiss_ratio = 0.1, sby_seed = 3L),
    regexp = "reteve zero linhas"
  )
})
