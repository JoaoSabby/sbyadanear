test_that("all public sampling routes share geometry and domains", {
  d <- fixture(); x <- as.matrix(d[-ncol(d)])
  for(model in 1:3) {
    h <- suppressWarnings(sby_adanear_hpc(d,target~.,sby_adasyn_k=3L,sby_nearmiss_k=2L,
      sby_adasyn_ratio=.75,nearmiss_model=model,sby_seed=41L,sby_audit=TRUE))
    c <- suppressWarnings(sby_adanear(target~.,d,sby_knn_over_k=3L,sby_knn_under_k=2L,
      sby_adasyn_ratio=.75,nearmiss_model=model,sby_seed=41L,sby_audit=TRUE))
    m <- suppressWarnings(sby_adanear_matrix(x,d$target,sby_knn_over_k=3L,sby_knn_under_k=2L,
      sby_adasyn_ratio=.75,nearmiss_model=model,sby_seed=41L,sby_audit=TRUE))
    expect_identical(values(h),values(c$sby_balanced_data))
    expect_equal(as.matrix(h[-ncol(h)]),m$sby_x_matrix,ignore_attr=TRUE)
    expect_identical(h$target,m$sby_y_vector)
  }
  out <- sby_balance_matrix(x,d$target,sby_strategy="none",sby_seed=3L)
  expect_equal(out$sby_x_matrix,x,ignore_attr=TRUE)
  expect_null(attr(out,"audit"))
})

test_that("explicit seed restores caller RNG, including error and absent seed", {
  d <- fixture(); set.seed(9); before <- .Random.seed; kind <- RNGkind()
  a <- sampler(d,sby_adasyn_ratio=1,sby_nearmiss_ratio=0)
  expect_identical(.Random.seed,before); expect_identical(RNGkind(),kind)
  RNGkind("L'Ecuyer-CMRG"); set.seed(8); before2 <- .Random.seed; kind2 <- RNGkind()
  b <- sampler(d,sby_adasyn_ratio=1,sby_nearmiss_ratio=0)
  expect_identical(values(a),values(b)); expect_identical(.Random.seed,before2)
  expect_identical(RNGkind(),kind2)
  separated <- data.frame(x=c(0,.1,10,11,12),target=factor(c("r","r","m","m","m")))
  runtime_before_error <- sby_hpc_cpu_report()$runtime
  expect_error(sby_adasyn_hpc(separated,target~.,sby_adasyn_ratio=1,sby_adasyn_k=1L,sby_seed=1L),"zero")
  expect_identical(.Random.seed,before2)
  runtime_after_error <- sby_hpc_cpu_report()$runtime
  for(nm in c("openmp_max_threads","openmp_dynamic","mkl_max_threads"))
    expect_identical(runtime_before_error[[nm]],runtime_after_error[[nm]])
  rm(".Random.seed",envir=.GlobalEnv)
  sampler(d,sby_adasyn_ratio=1,sby_nearmiss_ratio=0)
  expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
  do.call(RNGkind,as.list(kind)); assign(".Random.seed",before,.GlobalEnv)
})

test_that("thread limits do not change ADASYN output and controls are restored", {
  d <- fixture(); before <- sby_hpc_cpu_report()$runtime
  out <- lapply(c(1L,2L,4L),function(t) sby_adanear_hpc(d,target~.,
    sby_adasyn_ratio=1,nearmiss_model=1L,sby_seed=42L,
    sby_config_max_threads=t,sby_audit=TRUE))
  expect_identical(values(out[[1]]),values(out[[2]]))
  expect_identical(values(out[[1]]),values(out[[3]]))
  after <- sby_hpc_cpu_report()$runtime
  for(nm in c("openmp_max_threads","openmp_dynamic","mkl_max_threads")) expect_identical(before[[nm]],after[[nm]])
  for(i in seq_along(out)) {
    a <- attr(out[[i]],"audit")
    expect_lte(a$threads$effective,c(1L,2L,4L)[i])
    expect_true(all(a$telemetry$omp_observed_threads<=a$threads$effective))
  }
  if(identical(Sys.getenv("SBY_REQUIRE_INTEL"),"true")) {
    expect_match(after$compiler,"Intel")
    expect_match(after$openmp_library,"libiomp5")
    expect_gte(attr(out[[2]],"audit")$threads$effective,2L)
    expect_true(all(attr(out[[2]],"audit")$telemetry$omp_observed_threads[2:5]==2L))
  }
})

test_that("recipes supports all variants and preserves audited data return", {
  d <- fixture()
  for(fn in list(sby_step_nearmiss,sby_step_nearmiss_hpc,sby_step_adanear,sby_step_adanear_hpc)) {
    for(model in 1:3) {
      rec <- recipes::recipe(target~.,data=d)
      rec <- fn(rec,recipes::all_outcomes(),nearmiss_model=model,sby_seed=4L,
        sby_audit=TRUE,skip=FALSE)
      trained <- suppressWarnings(recipes::prep(rec,training=d))
      out <- suppressWarnings(recipes::bake(trained,new_data=d))
      expect_s3_class(out,"data.frame")
      expect_true(all(out$binary %in% 0:1))
      expect_identical(trained$steps[[1]]$parameters$nearmiss_model,model)
      expect_equal(trained$steps[[1]]$audit_log$last$audit$nearmiss$model,model)
      expect_match(trained$steps[[1]]$audit_log$last$sbyaudit$function_name,"step")
    }
  }
})
