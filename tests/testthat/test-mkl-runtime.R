test_that("mkl env vars are readable and config object is valid", {
  old <- options(sbyadanear.perf_mode = "auto")
  on.exit(options(old), add = TRUE)

  withr::local_envvar(c(
    OMP_NUM_THREADS = "2",
    MKL_NUM_THREADS = "2",
    OPENBLAS_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1"
  ))

  cfg <- sbyadanear:::sby_resolve_oneapi_mkl()
  expect_true(is.list(cfg))
  expect_true(all(c("enabled", "threads") %in% names(cfg)))
  expect_true(is.logical(cfg$enabled))
  expect_true(is.numeric(cfg$threads) || is.integer(cfg$threads))
  expect_identical(Sys.getenv("MKL_NUM_THREADS"), "2")
})
