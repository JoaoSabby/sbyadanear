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

test_that("hpc thread resolver honours the process affinity", {
  affinity <- sbyadanear:::sby_hpc_affinity()
  expect_type(affinity, "integer")
  if (length(affinity)) {
    expect_lte(sbyadanear:::sby_hpc_resolve_threads(8L), length(affinity))
  }
})

test_that("cpu report exposes affinity, runtime limits and thread environment", {
  report <- sby_hpc_cpu_report()
  expect_true(all(c(
    "affinity_cpus", "affinity_cpu_count", "openmp_max_threads",
    "mkl_max_threads", "requested_threads", "effective_threads",
    "affinity_warning", "hpc_environment"
  ) %in% names(report)))
  expect_true(all(c(
    "MKL_NUM_THREADS", "MKL_DOMAIN_NUM_THREADS", "OMP_NUM_THREADS"
  ) %in% names(report$hpc_environment)))
})

test_that("an affinity restricted to one CPU resolves eight threads to one", {
  expect_identical(
    sbyadanear:::sby_hpc_effective_thread_limit(
      sby_requested = 8L,
      sby_physical = 48L,
      sby_quota = NA_integer_,
      sby_affinity_count = 1L
    ),
    1L
  )
})

test_that("native thread guard applies limits and restores runtime state", {
  skip_if_not(sby_adanear_hpc_available())
  before <- sbyadanear:::sby_call_native("sby_hpc_compile_report_cpp")
  available <- before$affinity_cpu_count
  if (is.null(available) || is.na(available)) available <- parallel::detectCores()

  for (threads in c(1L, 2L, 8L)) {
    probe <- sbyadanear:::sby_call_native(
      "sby_hpc_thread_probe_cpp", threads, FALSE
    )
    if (isTRUE(before$openmp)) {
      expect_identical(as.integer(probe$openmp_max_threads),
                       min(threads, as.integer(available)))
    }
    if (isTRUE(before$mkl_linked)) {
      expect_identical(as.integer(probe$mkl_max_threads),
                       min(threads, as.integer(available)))
    } else {
      expect_true(is.na(probe$mkl_max_threads))
    }
  }

  after <- sbyadanear:::sby_call_native("sby_hpc_compile_report_cpp")
  expect_identical(after$openmp_max_threads, before$openmp_max_threads)
  expect_identical(after$mkl_max_threads, before$mkl_max_threads)

  expect_error(
    sbyadanear:::sby_call_native("sby_hpc_thread_probe_cpp", 1L, TRUE),
    "kernel abortado"
  )
  after_abort <- sbyadanear:::sby_call_native("sby_hpc_compile_report_cpp")
  expect_identical(after_abort$openmp_max_threads, before$openmp_max_threads)
  expect_identical(after_abort$mkl_max_threads, before$mkl_max_threads)
})

test_that("oneMKL local limit takes precedence without changing its environment", {
  skip_if_not(sby_adanear_hpc_available())
  report <- sbyadanear:::sby_call_native("sby_hpc_compile_report_cpp")
  skip_if_not(isTRUE(report$mkl_linked))
  withr::local_envvar(c(
    MKL_NUM_THREADS = "46",
    MKL_DOMAIN_NUM_THREADS = "BLAS=46"
  ))
  probe <- sbyadanear:::sby_call_native("sby_hpc_thread_probe_cpp", 2L, FALSE)
  expect_identical(as.integer(probe$mkl_max_threads), 2L)
  expect_identical(Sys.getenv("MKL_NUM_THREADS"), "46")
  expect_identical(Sys.getenv("MKL_DOMAIN_NUM_THREADS"), "BLAS=46")
})
