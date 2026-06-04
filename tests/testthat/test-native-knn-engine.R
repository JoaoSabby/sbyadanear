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
