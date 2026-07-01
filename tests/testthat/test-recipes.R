# End-to-end tests for the recipes::step wrappers around the three sampling
# routines. These tests intentionally exercise prep() -> bake() so that any
# regression in the S3 dispatch (prep / bake / tidy / required_pkgs / print)
# is caught early.

test_that("sby_step_adasyn balances a recipe end-to-end", {
  skip_if_not_installed("recipes")
  set.seed(11)
  dat <- data.frame(
    a = rnorm(40),
    b = rnorm(40),
    TARGET = factor(c(rep("min", 10), rep("maj", 30)),
                    levels = c("min", "maj"))
  )

  rec <- recipes::recipe(TARGET ~ ., data = dat)
  rec <- sby_step_adasyn(
    recipe                  = rec,
    recipes::all_outcomes(),
    sby_adasyn_ratio          = 0.5,
    sby_seed                = 1L,
    sby_knn_engine          = "FNN"
  )
  prepped <- recipes::prep(rec, training = dat)
  # ADASYN changes row count, so the canonical recipes pattern is to read
  # the trained training set via juice() (the step has skip = TRUE by
  # default, so bake(new_data = ...) intentionally returns rows unchanged).
  juiced <- recipes::juice(prepped)

  expect_s3_class(juiced, "tbl_df")
  expect_true("TARGET" %in% names(juiced))
  expect_gt(nrow(juiced), nrow(dat))
  expect_setequal(levels(juiced$TARGET), levels(dat$TARGET))
})

test_that("sby_step_nearmiss reduces the majority class via recipes", {
  skip_if_not_installed("recipes")
  set.seed(12)
  dat <- data.frame(
    a = rnorm(60),
    b = rnorm(60),
    TARGET = factor(c(rep("min", 12), rep("maj", 48)),
                    levels = c("min", "maj"))
  )

  rec <- recipes::recipe(TARGET ~ ., data = dat)
  rec <- sby_step_nearmiss(
    recipe          = rec,
    recipes::all_outcomes(),
    sby_nearmiss_ratio = 0.5,
    sby_seed        = 1L,
    sby_knn_engine  = "FNN"
  )
  prepped <- recipes::prep(rec, training = dat)
  juiced <- recipes::juice(prepped)

  expect_s3_class(juiced, "tbl_df")
  expect_lt(nrow(juiced), nrow(dat))
  expect_setequal(levels(juiced$TARGET), levels(dat$TARGET))
})

test_that("sby_step_adanear runs the full pipeline through recipes", {
  skip_if_not_installed("recipes")
  set.seed(13)
  dat <- data.frame(
    a = rnorm(50),
    b = rnorm(50),
    TARGET = factor(c(rep("min", 10), rep("maj", 40)),
                    levels = c("min", "maj"))
  )

  rec <- recipes::recipe(TARGET ~ ., data = dat)
  rec <- sby_step_adanear(
    recipe                  = rec,
    recipes::all_outcomes(),
    sby_adasyn_ratio          = 0.5,
    sby_nearmiss_ratio         = 0.8,
    sby_seed                = 1L,
    sby_knn_engine          = "FNN"
  )
  prepped <- recipes::prep(rec, training = dat)
  juiced <- recipes::juice(prepped)

  expect_s3_class(juiced, "tbl_df")
  # The original minority label must still be a level in the output.
  expect_true("min" %in% levels(juiced$TARGET))
  expect_true("maj" %in% levels(juiced$TARGET))
})

test_that("tidy() and required_pkgs() work on the three step types", {
  skip_if_not_installed("recipes")
  set.seed(14)
  dat <- data.frame(
    a = rnorm(20),
    b = rnorm(20),
    TARGET = factor(c(rep("min", 5), rep("maj", 15)),
                    levels = c("min", "maj"))
  )
  rec <- recipes::recipe(TARGET ~ ., data = dat) |>
    sby_step_adasyn(recipes::all_outcomes(), sby_seed = 1L)
  td <- generics::tidy(rec$steps[[1]])
  expect_s3_class(td, "data.frame")
  expect_true("sby_sampling_method" %in% names(td))
  pkgs <- generics::required_pkgs(rec$steps[[1]])
  expect_type(pkgs, "character")
  expect_true("recipes" %in% pkgs)
})
