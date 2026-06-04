# Tests for the conservative sby_resolve_knn_engine heuristic.

test_that("auto refuses non-euclidean approximate routes unless explicitly allowed", {
  expect_error(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "cosine", 100L, 5L)
    ),
    "automatico exato suporta apenas"
  )
  skip_if_not_installed("RcppHNSW")
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "ip", 100L, 5L, TRUE)
    ),
    "RcppHNSW"
  )
})

test_that("auto selects native when exact native routines are available", {
  skip_if_not(sby_adanear_native_available())
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "euclidean", 100000L, 100L)
    ),
    "native"
  )
  expect_equal(
    suppressMessages(
      sbyadanear:::sby_resolve_knn_engine("auto", 1L, "euclidean", 1000L, 5L)
    ),
    "native"
  )
})

test_that("auto preserves an explicit choice", {
  expect_equal(
    sbyadanear:::sby_resolve_knn_engine("native", 1L, "euclidean", 1000L, 5L),
    "native"
  )
  expect_equal(
    sbyadanear:::sby_resolve_knn_engine("FNN", 1L, "euclidean", 1000L, 5L),
    "FNN"
  )
  expect_equal(
    sbyadanear:::sby_resolve_knn_engine("RcppHNSW", 1L, "euclidean", 1000L, 5L),
    "RcppHNSW"
  )
})
