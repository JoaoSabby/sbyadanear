test_that("ADASYN returns synthetic row positions in the sby attribute", {
  set.seed(100)
  x <- matrix(rnorm(40), ncol = 2)
  y <- factor(c(rep("minor", 5L), rep("major", 15L)))

  result <- sby_adasyn_matrix(
    x, y,
    sby_adasyn_ratio = 0.6,
    sby_knn_over_k = 2L,
    sby_seed = 10L,
    sby_knn_engine = "FNN"
  )

  expect_identical(attr(result, "sby")$synthetic_rows, 21:23)
  expect_identical(
    result$sby_y_vector[attr(result, "sby")$synthetic_rows],
    factor(rep("minor", 3L), levels = levels(result$sby_y_vector))
  )
})

test_that("ADASYN skipped returns the integer zero sentinel", {
  x <- matrix(seq_len(20), ncol = 2)
  y <- factor(c(rep("minor", 3L), rep("major", 7L)))

  matrix_result <- sby_adasyn_matrix(x, y, sby_adasyn_ratio = 0)
  tabular_result <- sby_adasyn(
    target ~ ., data.frame(target = y, x1 = x[, 1], x2 = x[, 2]),
    sby_adasyn_ratio = 0
  )

  expect_identical(attr(matrix_result, "sby"), list(synthetic_rows = 0L))
  expect_identical(attr(tabular_result, "sby"), list(synthetic_rows = 0L))
})

test_that("ADANEAR reports positions after NearMiss reordering", {
  set.seed(200)
  x <- matrix(rnorm(48), ncol = 2)
  y <- factor(c(rep("minor", 6L), rep("major", 18L)))

  result <- sby_adanear_matrix(
    x, y,
    sby_adasyn_ratio = 0.5,
    sby_nearmiss_ratio = 1,
    sby_knn_over_k = 2L,
    sby_knn_under_k = 2L,
    sby_seed = 20L,
    sby_knn_engine = "FNN",
    sby_return_scaled = TRUE
  )
  synthetic_rows <- attr(result, "sby")$synthetic_rows

  expect_type(synthetic_rows, "integer")
  expect_identical(synthetic_rows, which(result$sby_retained_index > nrow(x)))
  expect_true(length(synthetic_rows) > 0L)
  expect_true(all(result$sby_y_vector[synthetic_rows] == "minor"))
})
