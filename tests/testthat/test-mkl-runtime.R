test_that("mkl env vars are readable and config object is valid", {
  old <- options(sbyadanear.perf_mode = "auto")
  on.exit(options(old), add = TRUE)

  withr::local_envvar(c(
    OMP_NUM_THREADS = "2",
    MKL_NUM_THREADS = "2"
  ))

  cfg <- sbyadanear:::sby_resolve_oneapi_mkl()
  expect_true(is.list(cfg))
  expect_true(all(c("enabled", "threads") %in% names(cfg)))
  expect_true(is.logical(cfg$enabled))
  expect_true(is.numeric(cfg$threads) || is.integer(cfg$threads))
  expect_identical(Sys.getenv("MKL_NUM_THREADS"), "2")
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "2")
})

test_that("hpc env control restores threading and stripe variables", {
  withr::local_envvar(c(
    MKL_NUM_THREADS = "3",
    OMP_NUM_THREADS = "4",
    MKL_NUM_STRIPES = "7",
    KMP_AFFINITY = "server-value"
  ))

  expect_identical(
    sbyadanear:::sby_hpc_env_keys(),
    c("MKL_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_STRIPES")
  )

  previous <- sbyadanear:::sby_hpc_capture_env()
  sbyadanear:::sby_hpc_apply_env(sby_total_threads = 2L)
  expect_identical(Sys.getenv("MKL_NUM_THREADS"), "2")
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "2")
  expect_identical(Sys.getenv("MKL_NUM_STRIPES"), "2")
  expect_identical(Sys.getenv("KMP_AFFINITY"), "server-value")

  sbyadanear:::sby_hpc_restore_env(previous)
  expect_identical(Sys.getenv("MKL_NUM_THREADS"), "3")
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "4")
  expect_identical(Sys.getenv("MKL_NUM_STRIPES"), "7")
  expect_identical(Sys.getenv("KMP_AFFINITY"), "server-value")
})


test_that("mkl stripe resolver adapts to class shape and matrix size", {
  expect_identical(
    sbyadanear:::sby_hpc_resolve_mkl_num_stripes(16L),
    4L
  )
  expect_identical(
    sbyadanear:::sby_hpc_resolve_mkl_num_stripes(16L, 200L, 100L, 10L),
    8L
  )
  expect_identical(
    sbyadanear:::sby_hpc_resolve_mkl_num_stripes(16L, 1000000L, 100L, 10L),
    16L
  )
  expect_identical(
    sbyadanear:::sby_hpc_resolve_mkl_num_stripes(16L, 100L, 250L, 10L),
    1L
  )
})
