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

test_that("hpc env control captures and restores threading variables", {
  withr::local_envvar(c(
    MKL_NUM_THREADS = "3",
    OMP_NUM_THREADS = "4",
    KMP_AFFINITY = "server-value"
  ))

  expect_identical(
    sbyadanear:::sby_hpc_env_keys(),
    c("MKL_NUM_THREADS", "OMP_NUM_THREADS")
  )

  previous <- sbyadanear:::sby_hpc_capture_env()
  Sys.setenv(MKL_NUM_THREADS = "2", OMP_NUM_THREADS = "2")

  sbyadanear:::sby_hpc_restore_env(previous)
  expect_identical(Sys.getenv("MKL_NUM_THREADS"), "3")
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "4")
  expect_identical(Sys.getenv("KMP_AFFINITY"), "server-value")
})

test_that("hpc thread resolver honours cgroup cpu quotas", {
  quota <- sbyadanear:::sby_hpc_cgroup_cpu_quota()
  expect_true(is.na(quota) || (is.integer(quota) && quota >= 1L))

  threads <- sbyadanear:::sby_hpc_resolve_threads(-1L)
  expect_true(is.integer(threads) && threads >= 1L)
  if (!is.na(quota)) {
    expect_lte(threads, quota)
  }

  # Um teto explicito nunca e ultrapassado, e entradas invalidas caem no
  # comportamento automatico em vez de propagar NA para o motor nativo.
  expect_lte(sbyadanear:::sby_hpc_resolve_threads(1L), 1L)
  expect_true(sbyadanear:::sby_hpc_resolve_threads(NA_integer_) >= 1L)
})
