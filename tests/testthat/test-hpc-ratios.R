test_that("HPC shortcuts honor configurable sampling ratios", {
  set.seed(123)
  dat <- data.frame(
    x1 = c(rnorm(4, -2), rnorm(10, 2)),
    x2 = c(rnorm(4, -2), rnorm(10, 2)),
    y = factor(c(rep("min", 4), rep("maj", 10)), levels = c("min", "maj"))
  )

  ada_low <- sby_adasyn_hpc(dat, y ~ ., sby_k_adasyn = 3,
                            sby_ratio_over = 0.5)
  ada_high <- sby_adasyn_hpc(dat, y ~ ., sby_k_adasyn = 3,
                             sby_ratio_over = 1)
  expect_equal(nrow(ada_low), 16L)
  expect_equal(nrow(ada_high), 18L)

  near_equal <- sby_nearmiss_hpc(dat, y ~ ., sby_k_nearmiss = 3,
                                 sby_ratio_under = 1)
  near_loose <- sby_nearmiss_hpc(dat, y ~ ., sby_k_nearmiss = 3,
                                 sby_ratio_under = 0.5)
  expect_equal(nrow(near_equal), 8L)
  expect_equal(nrow(near_loose), 6L)

  adanear_equal <- sby_adanear_hpc(dat, y ~ ., sby_k_adasyn = 3,
                                   sby_k_nearmiss = 3,
                                   sby_ratio_over = 0.5,
                                   sby_ratio_under = 1)
  adanear_loose <- sby_adanear_hpc(dat, y ~ ., sby_k_adasyn = 3,
                                   sby_k_nearmiss = 3,
                                   sby_ratio_over = 1,
                                   sby_ratio_under = 0.5)
  expect_equal(nrow(adanear_equal), 12L)
  expect_equal(nrow(adanear_loose), 12L)
  expect_equal(sum(adanear_equal$y == "min"), 6L)
  expect_equal(sum(adanear_loose$y == "min"), 8L)
})

test_that("HPC shortcuts validate sampling ratios before native execution", {
  dat <- data.frame(
    x1 = seq_len(14),
    x2 = seq_len(14) / 2,
    y = factor(c(rep("min", 4), rep("maj", 10)), levels = c("min", "maj"))
  )

  expect_error(sby_adasyn_hpc(dat, y ~ ., sby_ratio_over = 0),
               regexp = "sby_ratio_over")
  expect_error(sby_nearmiss_hpc(dat, y ~ ., sby_ratio_under = 0),
               regexp = "sby_ratio_under")
  expect_error(sby_nearmiss_hpc(dat, y ~ ., sby_ratio_under = -1),
               regexp = "sby_ratio_under")
  expect_error(sby_adanear_hpc(dat, y ~ ., sby_ratio_over = -1),
               regexp = "sby_ratio_over")
  expect_error(sby_adanear_hpc(dat, y ~ ., sby_ratio_under = -1),
               regexp = "sby_ratio_under")
  expect_equal(nrow(sby_nearmiss_hpc(dat, y ~ ., sby_ratio_under = 2)), 12L)
  expect_equal(nrow(sby_adanear_hpc(dat, y ~ ., sby_ratio_over = 0,
                                    sby_ratio_under = 0)), nrow(dat))
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
    sby_k_adasyn = 3,
    sby_k_nearmiss = 3,
    sby_ratio_over = 1,
    sby_ratio_under = 0.5
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

test_that("sby_adanear_hpc honors fixed seeds", {
  set.seed(456)
  dat <- data.frame(
    x1 = c(rnorm(5, -2), rnorm(15, 2)),
    x2 = c(rnorm(5, -2), rnorm(15, 2)),
    y = factor(c(rep("min", 5), rep("maj", 15)), levels = c("min", "maj"))
  )

  out1 <- sby_adanear_hpc(dat, y ~ ., sby_k_adasyn = 3,
                          sby_k_nearmiss = 3,
                          sby_seed = 99L,
                          sby_ratio_over = 1,
                          sby_ratio_under = 0.5)
  out2 <- sby_adanear_hpc(dat, y ~ ., sby_k_adasyn = 3,
                          sby_k_nearmiss = 3,
                          sby_seed = 99L,
                          sby_ratio_over = 1,
                          sby_ratio_under = 0.5)

  expect_identical(out1, out2)
  expect_error(sby_adanear_hpc(dat, y ~ ., sby_seed = 1.5),
               regexp = "sby_seed")
})
