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

test_that("sby_nearmiss_index exposes operational scale without full audit", {
  set.seed(1001)
  x <- cbind(rnorm(30), rnorm(30))
  storage.mode(x) <- "double"
  y <- factor(c(rep("minor", 6), rep("major", 24)), levels = c("minor", "major"))

  idx <- sby_nearmiss_index(
    x, y,
    sby_under_ratio = 0.5,
    sby_seed = 19,
    sby_knn_engine = "FNN",
    sby_return_scaling_info = TRUE,
    sby_return_reduced_scaled = TRUE
  )

  expect_true(is.integer(idx$sby_retained_index))
  expect_true(is.list(idx$sby_scaling_info))
  expect_equal(nrow(idx$sby_reduced_scaled), length(idx$sby_retained_index))
  expect_null(idx$sby_selected_majority_index)
})

test_that("sby_get_knnx honors requested return components", {
  set.seed(1002)
  x <- cbind(rnorm(20), rnorm(20))
  storage.mode(x) <- "double"

  both <- instenginer:::sby_get_knnx(x, x[1:5, , drop = FALSE], 3L, "kd_tree", "FNN", "euclidean", 1L, 16L, 200L, sby_knn_return = "both")
  index <- instenginer:::sby_get_knnx(x, x[1:5, , drop = FALSE], 3L, "kd_tree", "FNN", "euclidean", 1L, 16L, 200L, sby_knn_return = "index")
  dist <- instenginer:::sby_get_knnx(x, x[1:5, , drop = FALSE], 3L, "kd_tree", "FNN", "euclidean", 1L, 16L, 200L, sby_knn_return = "dist")

  expect_true(all(c("nn.index", "nn.dist") %in% names(both)))
  expect_true("nn.index" %in% names(index))
  expect_false("nn.dist" %in% names(index))
  expect_true("nn.dist" %in% names(dist))
  expect_false("nn.index" %in% names(dist))
})

test_that("sby_balance_matrix none and weight do not densify sparse matrices", {
  skip_if_not_installed("Matrix")
  x <- Matrix::rsparsematrix(20, 3, density = 0.2)
  y <- factor(c(rep("minor", 5), rep("major", 15)), levels = c("minor", "major"))

  none <- sby_balance_matrix(x, y, sby_strategy = "none")
  weight <- sby_balance_matrix(x, y, sby_strategy = "weight")

  expect_s4_class(none$sby_x_matrix, "dgCMatrix")
  expect_s4_class(weight$sby_x_matrix, "dgCMatrix")
  expect_identical(none$sby_x_matrix, x)
  expect_identical(weight$sby_x_matrix, x)
  expect_false(none$sby_diagnostics$sby_data_changed)
  expect_false(weight$sby_diagnostics$sby_data_changed)
  expect_error(sby_balance_matrix(x, y, sby_strategy = "adasyn"), "densa")
  expect_error(sby_balance_matrix(x, y, sby_strategy = "adanear"), "densa")
})

test_that("matrix audit levels avoid heavy ADANEAR intermediates unless full", {
  set.seed(1003)
  x <- cbind(rnorm(36), rnorm(36))
  storage.mode(x) <- "double"
  y <- factor(c(rep("minor", 9), rep("major", 27)), levels = c("minor", "major"))

  none <- sby_adanear_matrix(x, y, sby_seed = 21, sby_knn_engine = "FNN", sby_audit_level = "none")
  light <- sby_adanear_matrix(x, y, sby_seed = 21, sby_knn_engine = "FNN", sby_audit_level = "light")
  full <- sby_adanear_matrix(x, y, sby_seed = 21, sby_knn_engine = "FNN", sby_audit_level = "full")

  expect_null(none$sby_oversampling_result)
  expect_null(light$sby_oversampling_result)
  expect_equal(light$sby_diagnostics$sby_audit_level, "light")
  expect_true(is.list(full$sby_oversampling_result))
  expect_true(is.list(full$sby_undersampling_result))
})

test_that("tabular wrappers preserve original rows and restore only synthetic rows", {
  set.seed(1004)
  dat <- data.frame(
    int_col = as.integer(c(1:8, 101:124)),
    dbl_col = c(seq(0.1, 0.8, length.out = 8), seq(10.01, 12.4, length.out = 24)),
    bin_col = as.integer(c(rep(0, 4), rep(1, 4), rep(0, 12), rep(1, 12)))
  )
  dat$TARGET <- factor(c(rep("minor", 8), rep("major", 24)), levels = c("minor", "major"))

  ada <- sby_adasyn(TARGET ~ ., dat, sby_over_ratio = 0.5, sby_seed = 31, sby_knn_engine = "FNN")
  expect_identical(as.data.frame(ada[seq_len(nrow(dat)), names(dat)[1:3]]), dat[, 1:3])
  expect_type(ada$int_col, "integer")
  expect_type(ada$bin_col, "integer")

  near_audit <- sby_nearmiss(TARGET ~ ., dat, sby_under_ratio = 0.5, sby_seed = 32, sby_knn_engine = "FNN", sby_audit = TRUE)
  expect_identical(
    as.data.frame(near_audit$sby_balanced_data[, names(dat)[1:3]]),
    dat[near_audit$sby_retained_index, 1:3, drop = FALSE]
  )
})

test_that("native NearMiss selector matches R fallback without ties", {
  nn_dist <- matrix(
    c(0.10, 0.20,
      0.30, 0.40,
      0.05, 0.15,
      0.60, 0.70),
    nrow = 4,
    byrow = TRUE
  )
  storage.mode(nn_dist) <- "double"
  majority_index <- as.integer(c(10L, 20L, 30L, 40L))
  retained <- 2L

  c_idx <- .Call("OU_SelectNearMissMajorityC", nn_dist, majority_index, as.integer(retained), PACKAGE = "instenginer")
  r_idx <- majority_index[order(rowMeans(nn_dist), majority_index)[seq_len(retained)]]

  expect_equal(sort(c_idx), sort(r_idx))
})
