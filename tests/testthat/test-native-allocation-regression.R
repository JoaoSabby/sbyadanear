test_that("native brute path handles repeated zero-distance rescues", {
  skip_if_not(sby_adanear_native_available())

  x <- matrix(
    rep(c(0, 0, 1, 1, 2, 2), each = 2L),
    ncol = 2L,
    byrow = TRUE
  )
  y <- factor(c("minor", "minor", "major", "major", "major", "major"))

  old <- getOption("sbyadanear.sby_use_native_brute")
  on.exit(options(sbyadanear.sby_use_native_brute = old), add = TRUE)
  options(sbyadanear.sby_use_native_brute = TRUE)

  result <- sby_nearmiss_matrix(
    sby_x_matrix = x,
    sby_y_vector = y,
    sby_nearmiss_ratio = 1,
    sby_knn_under_k = 2L,
    sby_seed = 11L,
    sby_knn_algorithm = "brute",
    sby_knn_engine = "FNN",
    sby_knn_workers = 1L
  )

  expect_equal(nrow(result$sby_x_matrix), 4L)
  expect_equal(result$sby_output_class_distribution["minor"], 2L)
  expect_equal(result$sby_output_class_distribution["major"], 2L)
})

test_that("RcppParallel brute path returns stable shapes after vector reuse", {
  skip_if_not(sby_adanear_native_available())

  x <- matrix(as.double(seq_len(30L)), ncol = 3L)
  result <- sby_adasyn_matrix(
    sby_x_matrix = x,
    sby_y_vector = factor(c(rep("minor", 3L), rep("major", 7L))),
    sby_adasyn_ratio = 0.5,
    sby_knn_over_k = 2L,
    sby_seed = 13L,
    sby_knn_algorithm = "brute",
    sby_knn_engine = "FNN",
    sby_knn_parallel_backend = "RcppParallel",
    sby_knn_workers = 2L
  )

  expect_equal(ncol(result$sby_x_matrix), 3L)
  expect_equal(length(result$sby_y_vector), nrow(result$sby_x_matrix))
})
