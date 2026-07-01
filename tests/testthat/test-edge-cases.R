# Edge cases that exercise the contract of the matrix and tabular APIs.

test_that("matrix API rejects n_minority < 2 in both ADASYN and NearMiss", {
  x <- matrix(rnorm(40), 20, 2)
  storage.mode(x) <- "double"
  y_one <- factor(c("min", rep("maj", 19)), levels = c("min", "maj"))

  expect_error(sby_adasyn_matrix(x, y_one, sby_knn_engine = "FNN"),
               regexp = "minorit", ignore.case = TRUE)
  expect_error(sby_nearmiss_matrix(x, y_one, sby_knn_engine = "FNN"),
               regexp = "minorit", ignore.case = TRUE)
  expect_error(sby_nearmiss_index(x, y_one, sby_knn_engine = "FNN"),
               regexp = "minorit", ignore.case = TRUE)
})

test_that("output of ADASYN has no NA/NaN/Inf even with extreme ratio_over", {
  set.seed(1)
  x <- matrix(rnorm(2000), 1000, 2)
  storage.mode(x) <- "double"
  y <- factor(c(rep("min", 10), rep("maj", 990)), levels = c("min", "maj"))
  r <- sby_adasyn_matrix(x, y, sby_ratio_over = 50, sby_seed = 1L,
                         sby_knn_engine = "FNN")
  expect_false(anyNA(r$sby_x_matrix))
  expect_true(all(is.finite(r$sby_x_matrix)))
})

test_that("output of ADANEAR has no NA/NaN/Inf in extreme regimes", {
  set.seed(1)
  x <- matrix(rnorm(2000), 1000, 2)
  storage.mode(x) <- "double"
  y <- factor(c(rep("min", 10), rep("maj", 990)), levels = c("min", "maj"))
  r <- sby_adanear_matrix(x, y, sby_ratio_over = 5, sby_ratio_under = 1,
                          sby_seed = 1L, sby_knn_engine = "FNN")
  expect_false(anyNA(r$sby_x_matrix))
  expect_true(all(is.finite(r$sby_x_matrix)))
})

test_that("seed makes ADASYN, NearMiss and ADANEAR exactly reproducible", {
  set.seed(123)
  x <- matrix(rnorm(200), 100, 2)
  storage.mode(x) <- "double"
  y <- factor(c(rep("min", 20), rep("maj", 80)), levels = c("min", "maj"))

  ada1 <- sby_adasyn_matrix(x, y, sby_ratio_over = 0.5, sby_seed = 42L,
                            sby_knn_engine = "FNN")
  ada2 <- sby_adasyn_matrix(x, y, sby_ratio_over = 0.5, sby_seed = 42L,
                            sby_knn_engine = "FNN")
  expect_identical(ada1$sby_x_matrix, ada2$sby_x_matrix)
  expect_identical(ada1$sby_y_vector, ada2$sby_y_vector)

  near1 <- sby_nearmiss_matrix(x, y, sby_ratio_under = 0.5, sby_seed = 42L,
                               sby_knn_engine = "FNN")
  near2 <- sby_nearmiss_matrix(x, y, sby_ratio_under = 0.5, sby_seed = 42L,
                               sby_knn_engine = "FNN")
  expect_identical(near1$sby_x_matrix, near2$sby_x_matrix)
  expect_identical(near1$sby_y_vector, near2$sby_y_vector)

  ad1 <- sby_adanear_matrix(x, y, sby_ratio_over = 0.5, sby_ratio_under = 0.8,
                            sby_seed = 42L, sby_knn_engine = "FNN")
  ad2 <- sby_adanear_matrix(x, y, sby_ratio_over = 0.5, sby_ratio_under = 0.8,
                            sby_seed = 42L, sby_knn_engine = "FNN")
  expect_identical(ad1$sby_x_matrix, ad2$sby_x_matrix)
  expect_identical(ad1$sby_y_vector, ad2$sby_y_vector)
})

test_that("matrix API rejects NA, Inf and constant columns", {
  x <- matrix(rnorm(40), 20, 2)
  storage.mode(x) <- "double"
  y <- factor(c(rep("min", 5), rep("maj", 15)), levels = c("min", "maj"))

  x_na <- x; x_na[1, 1] <- NA_real_
  expect_error(sby_adasyn_matrix(x_na, y, sby_knn_engine = "FNN"), "NA")

  x_inf <- x; x_inf[1, 1] <- Inf
  expect_error(sby_adasyn_matrix(x_inf, y, sby_knn_engine = "FNN"), "Inf")

  x_const <- x; x_const[, 1] <- 1
  expect_error(sby_adasyn_matrix(x_const, y, sby_knn_engine = "FNN"),
               "desvio padrao zero|indefinido")
})

test_that("matrix API rejects multi-class and tied targets", {
  x <- matrix(rnorm(40), 20, 2)
  storage.mode(x) <- "double"
  y_three <- factor(rep(c("a", "b", "c", "a", "b"), each = 4),
                    levels = c("a", "b", "c"))
  expect_error(sby_adasyn_matrix(x, y_three, sby_knn_engine = "FNN"),
               regexp = "duas|dois", ignore.case = TRUE)

  y_tied <- factor(rep(c("a", "b"), each = 10), levels = c("a", "b"))
  expect_error(sby_adasyn_matrix(x, y_tied, sby_knn_engine = "FNN"),
               regexp = "desbalanceadas|minorit|2 obs", ignore.case = TRUE)
})

test_that("tabular API rejects logical, character, and ordered factor sensibly", {
  set.seed(1)
  dat <- data.frame(a = rnorm(20), b = rnorm(20))

  dat$y_char <- c(rep("x", 5), rep("y", 15))
  expect_s3_class(sby_adasyn(y_char ~ ., dat, sby_ratio_over = 0.5,
                             sby_seed = 1L, sby_knn_engine = "FNN"),
                  "tbl_df")

  dat$y_logical <- c(rep(FALSE, 5), rep(TRUE, 15))
  expect_s3_class(sby_adasyn(y_logical ~ ., dat, sby_ratio_over = 0.5,
                             sby_seed = 1L, sby_knn_engine = "FNN"),
                  "tbl_df")

  dat$y_ord <- factor(c(rep("L", 5), rep("H", 15)),
                      levels = c("L", "H"), ordered = TRUE)
  expect_s3_class(sby_adasyn(y_ord ~ ., dat, sby_ratio_over = 0.5,
                             sby_seed = 1L, sby_knn_engine = "FNN"),
                  "tbl_df")
})

test_that("tabular API gives a clear error when TARGET is reused as predictor", {
  set.seed(1)
  dat <- data.frame(x1 = rnorm(20), TARGET = rnorm(20),
                    y = factor(c(rep("a", 5), rep("b", 15)),
                               levels = c("a", "b")))
  expect_error(sby_adasyn(y ~ ., dat, sby_ratio_over = 0.5, sby_seed = 1L,
                          sby_knn_engine = "FNN"),
               regexp = "TARGET", ignore.case = FALSE)
})

test_that("exclusive zero ratios are rejected and ratio_under above one is accepted", {
  set.seed(1)
  x <- matrix(rnorm(40), 20, 2)
  storage.mode(x) <- "double"
  y <- factor(c(rep("min", 5), rep("maj", 15)), levels = c("min", "maj"))

  expect_error(sby_adasyn_matrix(x, y, sby_ratio_over = 0, sby_knn_engine = "FNN"),
               regexp = "positivo|positive", ignore.case = TRUE)
  expect_error(sby_nearmiss_matrix(x, y, sby_ratio_under = 0,
                                   sby_knn_engine = "FNN"),
               regexp = "maior que zero|positive", ignore.case = TRUE)
  expect_error(sby_nearmiss_matrix(x, y, sby_ratio_under = -1,
                                   sby_knn_engine = "FNN"),
               regexp = "maior que zero|positive", ignore.case = TRUE)
  expect_equal(
    as.integer(table(sby_nearmiss_matrix(x, y, sby_ratio_under = 2,
                                         sby_knn_engine = "FNN")$sby_y_vector)["maj"]),
    10L
  )
})

test_that("adanear zero ratios skip the corresponding stages", {
  set.seed(1)
  x <- matrix(rnorm(40), 20, 2)
  storage.mode(x) <- "double"
  y <- factor(c(rep("min", 5), rep("maj", 15)), levels = c("min", "maj"))

  no_over <- sby_adanear_matrix(x, y, sby_ratio_over = 0, sby_ratio_under = 1,
                                sby_knn_engine = "FNN", sby_return_scaled = TRUE)
  expect_false(no_over$sby_diagnostics$sby_adasyn_executed)
  expect_true(no_over$sby_diagnostics$sby_nearmiss_executed)
  expect_equal(as.integer(table(no_over$sby_y_vector)["maj"]), 5L)

  no_under <- sby_adanear_matrix(x, y, sby_ratio_over = 0.2, sby_ratio_under = 0,
                                 sby_knn_engine = "FNN", sby_return_scaled = TRUE)
  expect_true(no_under$sby_diagnostics$sby_adasyn_executed)
  expect_false(no_under$sby_diagnostics$sby_nearmiss_executed)
  expect_equal(as.integer(table(no_under$sby_y_vector)["maj"]), 15L)

  neither <- sby_adanear_matrix(x, y, sby_ratio_over = 0, sby_ratio_under = 0,
                                sby_knn_engine = "FNN", sby_return_scaled = TRUE)
  expect_false(neither$sby_diagnostics$sby_adasyn_executed)
  expect_false(neither$sby_diagnostics$sby_nearmiss_executed)
  expect_equal(nrow(neither$sby_x_matrix), nrow(x))
  expect_equal(as.character(neither$sby_y_vector), as.character(y))
})

test_that("k larger than minority count is handled gracefully", {
  set.seed(1)
  x <- matrix(rnorm(40), 20, 2)
  storage.mode(x) <- "double"
  y <- factor(c(rep("min", 3), rep("maj", 17)), levels = c("min", "maj"))
  r <- sby_adasyn_matrix(x, y, sby_ratio_over = 0.5, sby_knn_over_k = 50L,
                         sby_seed = 1L, sby_knn_engine = "FNN")
  expect_true(nrow(r$sby_x_matrix) >= 20L)
})
