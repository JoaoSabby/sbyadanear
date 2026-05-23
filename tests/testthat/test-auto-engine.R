# Tests for the new sby_resolve_knn_engine heuristic introduced in 0.4.0.

test_that("auto selects RcppHNSW for non-euclidean metrics", {
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "cosine", 100L, 5L)
    ),
    "RcppHNSW"
  )
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "ip", 100L, 5L)
    ),
    "RcppHNSW"
  )
})

test_that("auto selects RcppHNSW for large high-dimensional bases", {
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "euclidean",
                                            100000L, 100L)
    ),
    "RcppHNSW"
  )
})

test_that("auto selects FNN for small or low-dim bases", {
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "euclidean", 1000L, 5L)
    ),
    "FNN"
  )
  # Many rows but small dimensionality: still FNN.
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "euclidean",
                                            1000000L, 5L)
    ),
    "FNN"
  )
})

test_that("auto preserves an explicit choice", {
  expect_equal(
    sbyadanear:::sby_resolve_knn_engine("FNN", 1L, "euclidean", 1000L, 5L),
    "FNN"
  )
  expect_equal(
    sbyadanear:::sby_resolve_knn_engine("RcppHNSW", 1L, "euclidean",
                                          1000L, 5L),
    "RcppHNSW"
  )
})

test_that("the cells threshold is configurable via option", {
  # By default n=1000, p=60 gives 60000 cells, below the 5e6 default,
  # so auto must pick FNN.
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "euclidean",
                                            1000L, 60L)
    ),
    "FNN"
  )
  # Lower the threshold to 50000 and the same call must pick RcppHNSW
  # because n*p (60000) is now above the threshold AND p >= 50.
  old <- options(sbyadanear.sby_auto_engine_hnsw_min_cells = 50000)
  on.exit(options(old), add = TRUE)
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "euclidean",
                                            1000L, 60L)
    ),
    "RcppHNSW"
  )
})
