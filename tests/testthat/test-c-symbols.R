# Native C symbols must be reachable via the namespace object created by
# useDynLib(.registration = TRUE) and not by string + PACKAGE (which is
# unreachable because R_useDynamicSymbols(dll, FALSE) is set in
# R_init_sbyadanear). These tests also exercise the kernels directly to
# guard against regressions independent of the R glue layer.

test_that("C symbols are reachable via the namespace object", {
  expect_true(is.function(.Call) || TRUE)  # smoke

  z <- matrix(c(0, 1, 2, 3, 4, 5, 6, 7), 4, 2)
  storage.mode(z) <- "double"
  params <- .Call(sbyadanear:::compute_z_score_params_c, z)
  expect_named(params, c("centers", "scales"))
  expect_length(params$centers, 2)
  expect_length(params$scales, 2)
  expect_equal(params$centers, c(mean(z[, 1]), mean(z[, 2])))
  expect_equal(params$scales, c(sd(z[, 1]), sd(z[, 2])))
})

test_that("C symbol string form is intentionally unreachable", {
  expect_error(
    .Call("compute_z_score_params_c", matrix(0, 2, 2), PACKAGE = "sbyadanear"),
    regexp = "not available"
  )
})

test_that("select_nearmiss_majority_c matches R fallback with ties", {
  # Build a matrix where rows 1, 2 and 4 all have mean 0.15 (tie) and row 3
  # has mean 0.10. The smallest-mean wins; ties broken by smallest original
  # index. With retain=3, the native selector and the R fallback must agree on the
  # retained set.
  nn <- matrix(c(0.10, 0.20,
                 0.10, 0.20,
                 0.05, 0.15,
                 0.10, 0.20),
               nrow = 4, byrow = TRUE)
  storage.mode(nn) <- "double"
  maj <- as.integer(c(10L, 20L, 30L, 40L))

  c_idx <- .Call(sbyadanear:::select_nearmiss_majority_c, nn, maj, 3L)
  r_idx <- maj[order(rowMeans(nn), maj)[1:3]]
  expect_equal(sort(c_idx), sort(r_idx))
})

test_that("drop_self_neighbor_c removes self and preserves order", {
  nbr <- matrix(as.integer(c(
    1, 2, 3, 4,
    5, 2, 3, 4,
    3, 2, 1, 4,
    4, 5, 1, 2
  )), 4, 4, byrow = TRUE)
  self <- as.integer(c(1, 2, 3, 4))

  out <- .Call(sbyadanear:::drop_self_neighbor_c, nbr, self, 3L)
  expect_equal(dim(out), c(4L, 3L))
  # Each row must keep the first 3 valid neighbors that are not 'self'.
  expect_equal(out[1, ], c(2L, 3L, 4L))  # 1 -> dropped, keep 2 3 4
  expect_equal(out[2, ], c(5L, 3L, 4L))  # 2 -> dropped, keep 5 3 4
  expect_equal(out[3, ], c(2L, 1L, 4L))  # 3 -> dropped, keep 2 1 4
  expect_equal(out[4, ], c(5L, 1L, 2L))  # 4 -> dropped, keep 5 1 2
})

test_that("drop_self_neighbor_c errors when not enough valid candidates", {
  nbr <- matrix(as.integer(c(1L, NA_integer_, NA_integer_)), 1, 3)
  self <- 1L
  expect_error(
    .Call(sbyadanear:::drop_self_neighbor_c, nbr, self, 2L),
    regexp = "vizinhos suficientes"
  )
})

test_that("brute_force_knn_c matches FNN::get.knnx(algorithm='brute')", {
  skip_if_not_installed("FNN")
  set.seed(2)
  X <- matrix(rnorm(200 * 8), 200, 8)
  storage.mode(X) <- "double"
  Q <- X[1:50, , drop = FALSE]

  r_c <- .Call(sbyadanear:::brute_force_knn_c, X, Q, 5L)
  r_fnn <- FNN::get.knnx(data = X, query = Q, k = 5L, algorithm = "brute")

  expect_equal(r_c$nn.index, r_fnn$nn.index)
  expect_equal(r_c$nn.dist, r_fnn$nn.dist, tolerance = 1e-8)

  r_idx <- .Call(sbyadanear:::brute_force_knn_index_c, X, Q, 5L)
  r_dst <- .Call(sbyadanear:::brute_force_knn_dist_c, X, Q, 5L)

  expect_equal(r_idx$nn.index, r_fnn$nn.index)
  expect_null(r_idx$nn.dist)
  expect_null(r_dst$nn.index)
  expect_equal(r_dst$nn.dist, r_fnn$nn.dist, tolerance = 1e-8)
})

test_that("generate_synthetic_adasyn_col_c matches generate_synthetic_adasyn_c", {
  set.seed(42)
  minority <- matrix(rnorm(10 * 3), 10, 3)
  storage.mode(minority) <- "double"
  nbr <- matrix(as.integer(c(
    2, 3, 4,
    1, 3, 4,
    1, 2, 4,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3
  )), 10, 3, byrow = TRUE)
  per_row <- as.integer(c(2, 1, 1, 1, 1, 0, 0, 0, 0, 0))

  set.seed(123)
  r_row <- .Call(sbyadanear:::generate_synthetic_adasyn_c,
                 minority, nbr, per_row)
  set.seed(123)
  r_col <- .Call(sbyadanear:::generate_synthetic_adasyn_col_c,
                 minority, nbr, per_row)

  expect_equal(r_row, r_col)
})

test_that("brute_force_knn_c native path is used by sby_get_knnx by default", {
  skip_if_not_installed("FNN")
  set.seed(3)
  X <- matrix(rnorm(100 * 4), 100, 4)
  storage.mode(X) <- "double"
  Q <- X[1:10, , drop = FALSE]

  # Both options reach the same numeric result; just exercise both code paths.
  options(sbyadanear.sby_use_native_brute = TRUE)
  on.exit(options(sbyadanear.sby_use_native_brute = NULL), add = TRUE)
  r_native <- sbyadanear:::sby_get_knnx(
    X, Q, 5L, "brute", "FNN", "euclidean", 1L, 16L, 200L,
    sby_knn_return = "both"
  )

  options(sbyadanear.sby_use_native_brute = FALSE)
  r_fnn <- sbyadanear:::sby_get_knnx(
    X, Q, 5L, "brute", "FNN", "euclidean", 1L, 16L, 200L,
    sby_knn_return = "both"
  )

  expect_equal(r_native$nn.index, r_fnn$nn.index)
  expect_equal(r_native$nn.dist, r_fnn$nn.dist, tolerance = 1e-8)
})

test_that("nearmiss_brute_select_c matches exact NearMiss distance route", {
  set.seed(44)
  minority <- matrix(rnorm(12 * 5), 12, 5)
  majority <- matrix(rnorm(30 * 5), 30, 5)
  storage.mode(minority) <- "double"
  storage.mode(majority) <- "double"
  majority_index <- as.integer(seq(101L, 130L))

  knn <- .Call(sbyadanear:::brute_force_knn_dist_c, minority, majority, 4L)
  selected_reference <- .Call(
    sbyadanear:::select_nearmiss_majority_c,
    knn$nn.dist,
    majority_index,
    11L
  )
  selected_fused <- .Call(
    sbyadanear:::nearmiss_brute_select_c,
    minority,
    majority,
    majority_index,
    4L,
    11L
  )

  expect_equal(sort(selected_fused), sort(selected_reference))
})

test_that("brute_force_knn_c rescues near-zero BLAS cancellation", {
  data <- matrix(c(
    1e8 + 1e-3, 1e8,
    1e8 + 1e-4, 1e8
  ), nrow = 2, byrow = TRUE)
  query <- matrix(c(1e8, 1e8), nrow = 1)
  storage.mode(data) <- "double"
  storage.mode(query) <- "double"

  out <- .Call(sbyadanear:::brute_force_knn_c, data, query, 1L)
  exact_dist <- sqrt(sum((data[2, ] - query[1, ])^2))

  expect_equal(out$nn.index, matrix(2L, nrow = 1))
  expect_equal(out$nn.dist[1, 1], exact_dist, tolerance = 1e-12)
})

test_that("rbind_double_matrix_c matches base rbind and preserves column names", {
  first <- matrix(as.double(1:6), nrow = 3)
  second <- matrix(as.double(7:10), nrow = 2)
  colnames(first) <- c("a", "b")
  colnames(second) <- c("a", "b")

  out <- .Call(sbyadanear:::rbind_double_matrix_c, first, second)

  expect_equal(out, rbind(first, second))
  expect_equal(colnames(out), c("a", "b"))
})
