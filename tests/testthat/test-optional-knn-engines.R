test_that("KernelKnn engine returns the common KNN contract", {
  skip_if_not_installed("KernelKnn")

  set.seed(101)
  x <- matrix(stats::rnorm(60), ncol = 3)
  query <- x[1:4, , drop = FALSE]

  out <- sbyadanear:::sby_get_knnx(
    x, query, 3L, "auto", "KernelKnn", "euclidean", 1L, 16L, 200L,
    sby_knn_return = "both"
  )

  expect_equal(dim(out$nn.index), c(4L, 3L))
  expect_equal(dim(out$nn.dist), c(4L, 3L))
  expect_type(out$nn.index, "integer")
  expect_type(out$nn.dist, "double")
})

test_that("bigKNN engine returns the common KNN contract", {
  skip_if_not_installed("bigKNN")
  skip_if_not_installed("bigmemory")

  set.seed(102)
  x <- matrix(stats::rnorm(60), ncol = 3)
  query <- x[1:4, , drop = FALSE]

  out <- sbyadanear:::sby_get_knnx(
    x, query, 3L, "auto", "bigKNN", "euclidean", 1L, 16L, 200L,
    sby_knn_return = "both"
  )

  expect_equal(dim(out$nn.index), c(4L, 3L))
  expect_equal(dim(out$nn.dist), c(4L, 3L))
  expect_type(out$nn.index, "integer")
  expect_type(out$nn.dist, "double")
})

test_that("external engines reject unsupported metrics explicitly", {
  expect_error(
    sbyadanear:::sby_get_knnx(
      matrix(1:12, ncol = 2) + 0,
      matrix(1:4, ncol = 2) + 0,
      1L, "auto", "KernelKnn", "cosine", 1L, 16L, 200L
    ),
    "KernelKnn"
  )
  expect_error(
    sbyadanear:::sby_get_knnx(
      matrix(1:12, ncol = 2) + 0,
      matrix(1:4, ncol = 2) + 0,
      1L, "auto", "bigKNN", "cosine", 1L, 16L, 200L
    ),
    "bigKNN"
  )
})
