test_that("native engine matches FNN brute for both/index/dist", {
  skip_if_not(sby_adanear_native_available())
  skip_if_not_installed("FNN")
  set.seed(404)

  x <- matrix(rnorm(80L * 6L), nrow = 80L, ncol = 6L)
  q <- x[seq_len(20L), , drop = FALSE]
  storage.mode(x) <- "double"
  storage.mode(q) <- "double"

  ref <- FNN::get.knnx(data = x, query = q, k = 5L, algorithm = "brute")
  both <- sbyadanear:::sby_get_knnx(
    x, q, 5L, "auto", "native", "euclidean", 1L, 16L, 200L,
    sby_knn_return = "both"
  )
  index <- sbyadanear:::sby_get_knnx(
    x, q, 5L, "auto", "native", "euclidean", 1L, 16L, 200L,
    sby_knn_return = "index"
  )
  dist <- sbyadanear:::sby_get_knnx(
    x, q, 5L, "auto", "native", "euclidean", 1L, 16L, 200L,
    sby_knn_return = "dist"
  )

  expect_equal(both$nn.index, ref$nn.index)
  expect_equal(both$nn.dist, ref$nn.dist, tolerance = 1e-8)
  expect_equal(index$nn.index, ref$nn.index)
  expect_null(index$nn.dist)
  expect_null(dist$nn.index)
  expect_equal(dist$nn.dist, ref$nn.dist, tolerance = 1e-8)
})

test_that("native engine can exclude self-neighbors deterministically", {
  skip_if_not(sby_adanear_native_available())

  x <- matrix(c(
    0, 0,
    0, 0,
    1, 0,
    2, 0,
    3, 0
  ), ncol = 2L, byrow = TRUE)
  storage.mode(x) <- "double"

  out <- sbyadanear:::sby_get_knnx(
    x, x, 2L, "auto", "native", "euclidean", 1L, 16L, 200L,
    sby_knn_query_chunk_size = 2L,
    sby_query_is_data = TRUE,
    sby_exclude_self = TRUE,
    sby_knn_return = "both"
  )

  expect_false(any(out$nn.index[, 1L] == seq_len(nrow(x))))
  expect_equal(out$nn.index[1L, 1L], 2L)
  expect_equal(out$nn.index[2L, 1L], 1L)
  expect_equal(out$nn.dist[1L, 1L], 0)
})

test_that("native engine validates non-finite values", {
  skip_if_not(sby_adanear_native_available())

  x <- matrix(rnorm(20L), ncol = 2L)
  storage.mode(x) <- "double"
  x[1L, 1L] <- Inf

  expect_error(
    sbyadanear:::sby_get_knnx(
      x, x, 1L, "auto", "native", "euclidean", 1L, 16L, 200L,
      sby_knn_return = "both"
    ),
    "NA, NaN, Inf"
  )
})

test_that("FNN auto algorithm is resolved before calling FNN", {
  skip_if_not_installed("FNN")
  set.seed(405)
  x <- matrix(rnorm(20L * 5L), nrow = 20L, ncol = 5L)
  q <- x[seq_len(3L), , drop = FALSE]
  storage.mode(x) <- "double"
  storage.mode(q) <- "double"

  old <- getOption("sbyadanear.sby_use_native_brute")
  on.exit(options(sbyadanear.sby_use_native_brute = old), add = TRUE)
  options(sbyadanear.sby_use_native_brute = FALSE)

  out <- sbyadanear:::sby_get_knnx(
    x, q, 3L, "auto", "FNN", "euclidean", 1L, 16L, 200L,
    sby_knn_return = "both"
  )

  expect_equal(dim(out$nn.index), c(3L, 3L))
  expect_equal(dim(out$nn.dist), c(3L, 3L))
})

test_that("FNN native brute honors exclude_self like explicit native", {
  skip_if_not(sby_adanear_native_available())

  x <- matrix(c(
    0, 0,
    1, 0,
    2, 0,
    3, 0,
    4, 0
  ), ncol = 2L, byrow = TRUE)
  storage.mode(x) <- "double"

  old <- getOption("sbyadanear.sby_use_native_brute")
  on.exit(options(sbyadanear.sby_use_native_brute = old), add = TRUE)
  options(sbyadanear.sby_use_native_brute = TRUE)

  fnn_shortcut <- sbyadanear:::sby_get_knnx(
    x, x, 2L, "brute", "FNN", "euclidean", 1L, 16L, 200L,
    sby_query_is_data = TRUE,
    sby_exclude_self = TRUE,
    sby_knn_return = "both"
  )
  explicit_native <- sbyadanear:::sby_get_knnx(
    x, x, 2L, "brute", "native", "euclidean", 1L, 16L, 200L,
    sby_query_is_data = TRUE,
    sby_exclude_self = TRUE,
    sby_knn_return = "both"
  )

  expect_false(any(fnn_shortcut$nn.index[, 1L] == seq_len(nrow(x))))
  expect_identical(fnn_shortcut$nn.index, explicit_native$nn.index)
  expect_equal(fnn_shortcut$nn.dist, explicit_native$nn.dist, tolerance = 1e-12)
})

test_that("native chunking preserves partial return shapes across multiple chunks", {
  skip_if_not(sby_adanear_native_available())
  set.seed(406)
  x <- matrix(rnorm(17L * 4L), nrow = 17L, ncol = 4L)
  q <- x[seq_len(11L), , drop = FALSE]
  storage.mode(x) <- "double"
  storage.mode(q) <- "double"

  idx <- sbyadanear:::sby_get_knnx(
    x, q, 3L, "brute", "native", "euclidean", 1L, 16L, 200L,
    sby_knn_query_chunk_size = 4L,
    sby_knn_return = "index"
  )
  dst <- sbyadanear:::sby_get_knnx(
    x, q, 3L, "brute", "native", "euclidean", 1L, 16L, 200L,
    sby_knn_query_chunk_size = 4L,
    sby_knn_return = "dist"
  )
  both <- sbyadanear:::sby_get_knnx(
    x, q, 3L, "brute", "native", "euclidean", 1L, 16L, 200L,
    sby_knn_query_chunk_size = 4L,
    sby_knn_return = "both"
  )

  expect_named(idx, "nn.index")
  expect_named(dst, "nn.dist")
  expect_equal(dim(idx$nn.index), c(11L, 3L))
  expect_equal(dim(dst$nn.dist), c(11L, 3L))
  expect_identical(idx$nn.index, both$nn.index)
  expect_equal(dst$nn.dist, both$nn.dist, tolerance = 1e-12)
})

test_that("KNN integer validators reject fractional values", {
  expect_error(sbyadanear:::sby_validate_knn_workers(2.5), "inteiro positivo")
  expect_error(sbyadanear:::sby_validate_knn_query_chunk_size(10.5), "inteiro positivo")
  expect_error(sbyadanear:::sby_validate_hnsw_params(16.5, 200L), "inteiro >= 2")
  expect_error(sbyadanear:::sby_validate_hnsw_params(16L, 200.5), "inteiro positivo")
})
