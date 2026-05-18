test_that("sby_balance_matrix none and weight keep data unchanged", {
  x <- matrix(as.double(seq_len(40)), ncol = 2)
  y <- factor(c(rep("minor", 5), rep("major", 15)), levels = c("minor", "major"))

  none <- sby_balance_matrix(x, y, sby_strategy = "none")
  weight <- sby_balance_matrix(x, y, sby_strategy = "weight")

  expect_identical(none$sby_x_matrix, x)
  expect_identical(none$sby_y_vector, y)
  expect_identical(weight$sby_x_matrix, x)
  expect_identical(weight$sby_y_vector, y)
  expect_equal(weight$sby_class_ratio_output, 3)
})

test_that("matrix validators reject invalid dense inputs", {
  x <- matrix(as.double(seq_len(20)), ncol = 2)
  y <- factor(c(rep("minor", 3), rep("major", 7)), levels = c("minor", "major"))

  x_na <- x
  x_na[1, 1] <- NA_real_
  expect_error(sby_balance_matrix(x_na, y, sby_strategy = "none"), "NA")

  x_inf <- x
  x_inf[1, 1] <- Inf
  expect_error(sby_balance_matrix(x_inf, y, sby_strategy = "none"), "Inf")

  expect_error(sby_balance_matrix(x, factor(rep(c("a", "b", "c", "a", "b"), each = 2)), sby_strategy = "none"), "dois niveis")
  expect_error(sby_balance_matrix(x, y[-1], sby_strategy = "none"), "comprimento")
})

test_that("sby_adasyn_matrix increases only the minority class", {
  set.seed(123)
  x <- cbind(rnorm(30), rnorm(30))
  storage.mode(x) <- "double"
  y <- factor(c(rep("minor", 6), rep("major", 24)), levels = c("minor", "major"))

  out <- sby_adasyn_matrix(x, y, sby_over_ratio = 0.5, sby_seed = 11, sby_knn_engine = "FNN")

  expect_equal(unname(out$sby_input_class_distribution), c(6L, 24L))
  expect_equal(unname(out$sby_output_class_distribution), c(9L, 24L))
  expect_s3_class(out$sby_y_vector, "factor")
  expect_equal(ncol(out$sby_x_matrix), ncol(x))
})

test_that("sby_nearmiss_index agrees with sby_nearmiss_matrix retained rows", {
  set.seed(321)
  x <- cbind(rnorm(40), rnorm(40))
  storage.mode(x) <- "double"
  y <- factor(c(rep("minor", 8), rep("major", 32)), levels = c("minor", "major"))

  idx <- sby_nearmiss_index(x, y, sby_under_ratio = 0.5, sby_seed = 9, sby_knn_engine = "FNN")
  mat <- sby_nearmiss_matrix(x, y, sby_under_ratio = 0.5, sby_seed = 9, sby_knn_engine = "FNN", sby_return_index = TRUE)

  expect_identical(idx$sby_retained_index, mat$sby_retained_index)
  expect_equal(mat$sby_x_matrix, x[idx$sby_retained_index, , drop = FALSE])
  expect_equal(unname(mat$sby_output_class_distribution), c(8L, 16L))
})

test_that("sby_adanear_matrix preserves original class roles across stages", {
  set.seed(456)
  x <- cbind(rnorm(30), rnorm(30))
  storage.mode(x) <- "double"
  y <- factor(c(rep("minor", 10), rep("major", 20)), levels = c("minor", "major"))

  out <- sby_adanear_matrix(
    x, y,
    sby_over_ratio = 1.5,
    sby_under_ratio = 1,
    sby_seed = 7,
    sby_knn_engine = "FNN",
    sby_audit = TRUE
  )

  expect_equal(out$sby_diagnostics$sby_original_minority_label, "minor")
  expect_equal(out$sby_diagnostics$sby_original_majority_label, "major")
  expect_true(out$sby_output_class_distribution[["major"]] <= 20L)
})

test_that("audit levels control heavy matrix intermediates", {
  set.seed(654)
  x <- cbind(rnorm(24), rnorm(24))
  storage.mode(x) <- "double"
  y <- factor(c(rep("minor", 6), rep("major", 18)), levels = c("minor", "major"))

  none <- sby_adasyn_matrix(x, y, sby_seed = 12, sby_knn_engine = "FNN", sby_audit_level = "none")
  full <- sby_adasyn_matrix(x, y, sby_seed = 12, sby_knn_engine = "FNN", sby_audit_level = "full")
  scaled <- sby_adasyn_matrix(x, y, sby_seed = 12, sby_knn_engine = "FNN", sby_return_scaled = TRUE)

  expect_null(none$sby_balanced_scaled)
  expect_null(none$sby_scaling_info)
  expect_true(is.list(full$sby_scaling_info))
  expect_true(is.list(scaled$sby_balanced_scaled))
})

test_that("tabular wrappers keep class distributions compatible with matrix API", {
  set.seed(777)
  x <- cbind(a = rnorm(36), b = rnorm(36))
  storage.mode(x) <- "double"
  y <- factor(c(rep("minor", 9), rep("major", 27)), levels = c("minor", "major"))
  dat <- tibble::as_tibble(as.data.frame(x))
  dat$TARGET <- y

  mat_ada <- sby_adasyn_matrix(x, y, sby_over_ratio = 0.4, sby_seed = 4, sby_knn_engine = "FNN")
  tab_ada <- sby_adasyn(TARGET ~ ., dat, sby_over_ratio = 0.4, sby_seed = 4, sby_knn_engine = "FNN")
  expect_equal(unname(table(tab_ada$TARGET)), unname(mat_ada$sby_output_class_distribution))

  mat_near <- sby_nearmiss_matrix(x, y, sby_under_ratio = 0.5, sby_seed = 5, sby_knn_engine = "FNN")
  tab_near <- sby_nearmiss(TARGET ~ ., dat, sby_under_ratio = 0.5, sby_seed = 5, sby_knn_engine = "FNN")
  expect_equal(unname(table(tab_near$TARGET)), unname(mat_near$sby_output_class_distribution))

  mat_adanear <- sby_adanear_matrix(x, y, sby_over_ratio = 0.4, sby_under_ratio = 0.5, sby_seed = 6, sby_knn_engine = "FNN")
  tab_adanear <- sby_adanear(TARGET ~ ., dat, sby_over_ratio = 0.4, sby_under_ratio = 0.5, sby_seed = 6, sby_knn_engine = "FNN")
  expect_equal(unname(table(tab_adanear$TARGET)), unname(mat_adanear$sby_output_class_distribution))
})

test_that("ADASYN and ADANEAR reject sparse matrices explicitly", {
  testthat::skip_if_not_installed("Matrix")
  x_sparse <- Matrix::Matrix(matrix(as.double(seq_len(20)), ncol = 2), sparse = TRUE)
  y <- factor(c(rep("minor", 3), rep("major", 7)), levels = c("minor", "major"))

  expect_error(sby_adasyn_matrix(x_sparse, y), "esparsas")
  expect_error(sby_adanear_matrix(x_sparse, y), "esparsas")
})
