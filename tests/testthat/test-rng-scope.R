test_that("sby_adasyn_matrix preserves caller RNG state", {
  set.seed(2024)
  x <- matrix(rnorm(40), ncol = 2)
  y <- factor(rep(c("minor", "major"), c(6, 14)))

  set.seed(9001)
  old_kind <- RNGkind()
  old_seed <- .Random.seed

  out1 <- sby_adasyn_matrix(
    x,
    y,
    sby_ratio_over = 0.5,
    sby_seed = 42L,
    sby_knn_engine = "FNN"
  )
  out2 <- sby_adasyn_matrix(
    x,
    y,
    sby_ratio_over = 0.5,
    sby_seed = 42L,
    sby_knn_engine = "FNN"
  )

  expect_identical(RNGkind(), old_kind)
  expect_identical(.Random.seed, old_seed)
  expect_equal(out1$sby_x_matrix, out2$sby_x_matrix)
  expect_identical(out1$sby_y_vector, out2$sby_y_vector)
})

test_that("sby_nearmiss_index deterministic route does not consume RNG", {
  set.seed(2025)
  x <- matrix(rnorm(40), ncol = 2)
  y <- factor(rep(c("minor", "major"), c(6, 14)))

  set.seed(9002)
  old_kind <- RNGkind()
  old_seed <- .Random.seed

  idx1 <- sby_nearmiss_index(
    x,
    y,
    sby_ratio_under = 0.5,
    sby_seed = 10L,
    sby_knn_engine = "FNN"
  )
  idx2 <- sby_nearmiss_index(
    x,
    y,
    sby_ratio_under = 0.5,
    sby_seed = 999L,
    sby_knn_engine = "FNN"
  )

  expect_identical(RNGkind(), old_kind)
  expect_identical(.Random.seed, old_seed)
  expect_identical(idx1$sby_retained_index, idx2$sby_retained_index)
})
