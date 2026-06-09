test_that("HPC shortcuts honor configurable sampling ratios", {
  set.seed(123)
  dat <- data.frame(
    x1 = c(rnorm(4, -2), rnorm(10, 2)),
    x2 = c(rnorm(4, -2), rnorm(10, 2)),
    y = factor(c(rep("min", 4), rep("maj", 10)), levels = c("min", "maj"))
  )

  ada_low <- sby_adasyn_hpc(dat, y ~ ., sby_k_neighbor_adanear = 3,
                            sby_over_ratio = 0.5)
  ada_high <- sby_adasyn_hpc(dat, y ~ ., sby_k_neighbor_adanear = 3,
                             sby_over_ratio = 1)
  expect_equal(nrow(ada_low), 16L)
  expect_equal(nrow(ada_high), 18L)

  near_equal <- sby_nearmiss_hpc(dat, y ~ ., sby_k_neighbor_nearmiss = 3,
                                 sby_under_ratio = 1)
  near_loose <- sby_nearmiss_hpc(dat, y ~ ., sby_k_neighbor_nearmiss = 3,
                                 sby_under_ratio = 0.5)
  expect_equal(nrow(near_equal), 8L)
  expect_equal(nrow(near_loose), 12L)

  adanear_equal <- sby_adanear_hpc(dat, y ~ ., sby_k_neighbor_adanear = 3,
                                   sby_k_neighbor_nearmiss = 3,
                                   sby_over_ratio = 0.5,
                                   sby_under_ratio = 1)
  adanear_loose <- sby_adanear_hpc(dat, y ~ ., sby_k_neighbor_adanear = 3,
                                   sby_k_neighbor_nearmiss = 3,
                                   sby_over_ratio = 1,
                                   sby_under_ratio = 0.5)
  expect_equal(nrow(adanear_equal), 12L)
  expect_equal(nrow(adanear_loose), 18L)
})

test_that("HPC shortcuts validate sampling ratios before native execution", {
  dat <- data.frame(
    x1 = seq_len(14),
    x2 = seq_len(14) / 2,
    y = factor(c(rep("min", 4), rep("maj", 10)), levels = c("min", "maj"))
  )

  expect_error(sby_adasyn_hpc(dat, y ~ ., sby_over_ratio = 0),
               regexp = "sby_over_ratio")
  expect_error(sby_nearmiss_hpc(dat, y ~ ., sby_under_ratio = 2),
               regexp = "sby_under_ratio")
  expect_error(sby_adanear_hpc(dat, y ~ ., sby_over_ratio = 0),
               regexp = "sby_over_ratio")
  expect_error(sby_adanear_hpc(dat, y ~ ., sby_under_ratio = 2),
               regexp = "sby_under_ratio")
})

test_that("sby_adanear_hpc restores original scale and integer predictor types", {
  skip_if_not(sby_adanear_hpc_available())
  set.seed(321)
  dat <- data.frame(
    int_col = as.integer(c(1:5, 20:34)),
    class_col = factor(c(rep("rare", 5), rep("common", 15)),
                       levels = c("rare", "common")),
    dbl_col = c(rnorm(5, -10, 0.25), rnorm(15, 10, 0.25))
  )

  out <- sby_adanear_hpc(
    dat,
    class_col ~ .,
    sby_k_neighbor_adanear = 3,
    sby_k_neighbor_nearmiss = 3,
    sby_over_ratio = 1,
    sby_under_ratio = 0.5
  )

  expect_s3_class(out, "tbl_df")
  expect_true(is.integer(out$int_col))
  expect_type(out$dbl_col, "double")
  expect_true(all(out$int_col >= min(dat$int_col)))
  expect_true(all(out$int_col <= max(dat$int_col)))
  expect_identical(names(out), names(dat))
  expect_true("class_col" %in% names(out))
  expect_false("TARGET" %in% names(out))
  expect_setequal(levels(out$class_col), levels(dat$class_col))

  expect_true(any(out$dbl_col > 1))
  expect_false(all(abs(out$dbl_col) < 4))
})
