test_that("oneAPI MKL config is automatic and supports perf_mode", {
  old <- options(sbyadanear.perf_mode = "auto")
  on.exit(options(old), add = TRUE)

  cfg <- sbyadanear:::sby_resolve_oneapi_mkl()
  expect_true(is.list(cfg))
  expect_true(all(c("enabled", "threads") %in% names(cfg)))
})

test_that("manual perf_mode disables automatic MKL thread policy", {
  old <- options(sbyadanear.perf_mode = "manual")
  on.exit(options(old), add = TRUE)

  cfg <- sbyadanear:::sby_resolve_oneapi_mkl()
  expect_false(cfg$enabled)
})
