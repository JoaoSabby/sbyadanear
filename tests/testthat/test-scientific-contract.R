test_that("NearMiss variants agree with an independent Euclidean oracle", {
  d <- fixture(); x <- as.matrix(d[setdiff(names(d),"target")])
  for(model in 1:3) {
    out <- suppressWarnings(sby_nearmiss_hpc(d,target~.,sby_nearmiss_k=2L,
      sby_nearmiss_ratio=.75,nearmiss_model=model,sby_nearmiss_m=2L,
      sby_seed=3L,sby_audit=TRUE,sby_config_max_threads=2L))
    a <- attr(out,"audit")
    z <- sweep(sweep(x,2,a$initial_scale$sby_center,`-`),2,a$initial_scale$sby_scale,`/`)
    mi <- which(d$target=="rare"); ma <- which(d$target=="major")
    expected <- oracle_nearmiss(z,z[mi,,drop=FALSE],z[ma,,drop=FALSE],2L,model,2L,3L)
    expect_identical(a$original_majority_indices,as.integer(ma[expected]))
    expect_identical(values(out[seq_along(a$original_indices),]),values(d[a$original_indices,]))
  }
})

test_that("NearMiss-1 averages d rather than d squared, including k=all", {
  x <- rbind(c(0,0),c(4,0),c(0,0),c(2,1.5),c(20,20))
  y <- factor(c("rare","rare","major","major","major"))
  a <- sby_nearmiss_index(x,y,sby_nearmiss_ratio=.5,sby_knn_under_k=2L,
    nearmiss_model=1L,sby_scaling_info=list(sby_center=c(0,0),sby_scale=c(1,1)),
    sby_audit=TRUE,sby_seed=1L)
  expect_identical(a$sby_selected_majority_index,3L)
  expect_equal(attr(a,"audit")$nearmiss$scores[1:2],c(2,2.5),tolerance=1e-13)
})

test_that("the full pipeline uses the single initial scale and continuous synthetic rare rows", {
  d <- fixture(); x <- as.matrix(d[-ncol(d)])
  centers <- colMeans(x)
  scales <- sqrt(colMeans(sweep(x,2,centers,`-`)^2)); scales[scales==0] <- 1
  z <- sweep(sweep(x,2,centers,`-`),2,scales,`/`)
  mi <- which(d$target=="rare"); ma <- which(d$target=="major")
  for(model in 1:3) {
    out <- suppressWarnings(sampler(d,sby_adasyn_ratio=1.25,sby_nearmiss_ratio=.4,
      nearmiss_model=model,sby_nearmiss_m=2L,sby_audit=TRUE))
    a <- attr(out,"audit")
    expect_equal(a$initial_scale$sby_center,centers,tolerance=1e-14)
    expect_equal(a$initial_scale$sby_scale,scales,tolerance=1e-14)
    rare <- rbind(z[mi,,drop=FALSE],a$synthetic_standardized)
    expected <- oracle_nearmiss(z,rare,z[ma,,drop=FALSE],2L,model,2L,3L)
    expect_identical(a$original_majority_indices,as.integer(ma[expected]))
    expect_true(all(mi %in% a$original_indices))
    expect_identical(values(out[seq_along(a$original_indices),]),values(d[a$original_indices,]))
    expect_equal(nrow(out),length(a$original_indices)+5L)
  }
})

test_that("ADASYN traces reproduce continuous interpolation and quotas exactly", {
  d <- fixture()
  out <- sampler(d,sby_adasyn_ratio=1.25,sby_nearmiss_ratio=0,sby_audit=TRUE)
  a <- attr(out,"audit"); tr <- a$synthetics
  x <- as.matrix(d[setdiff(names(d),"target")])
  expected <- x[tr$parent,,drop=FALSE]+tr$lambda*(x[tr$neighbor,,drop=FALSE]-x[tr$parent,,drop=FALSE])
  expect_equal(a$synthetic_continuous,expected,tolerance=1e-13,ignore_attr=TRUE)
  expect_equal(sum(a$adasyn$quotas),5L)
  expect_true(all(tr$parent != tr$neighbor))
  expect_true(all(d$target[tr$parent]=="rare" & d$target[tr$neighbor]=="rare"))
  expect_true(all(tr$lambda>=0 & tr$lambda<1))
  expect_true(all(out$binary %in% 0:1))
  expect_type(out$integer,"integer")
  expect_true(all(out$integer==round(out$integer)))
  expect_true(all(out$whole_double==round(out$whole_double)))
  expect_identical(out$constant,rep(1L,nrow(out)))
  for(nm in colnames(x)) expect_true(all(out[[nm]]>=min(d[[nm]]) & out[[nm]]<=max(d[[nm]])))
  expect_identical(values(out[seq_len(nrow(d)),]),values(d))
  expect_identical(levels(out$target),levels(d$target))
  expect_true(is.ordered(out$target))
})

test_that("difficulty and largest-remainder quotas agree with direct original KNN", {
  d <- fixture(); out <- sampler(d,sby_adasyn_ratio=1.25,sby_nearmiss_ratio=0,sby_audit=TRUE)
  a <- attr(out,"audit"); x <- as.matrix(d[setdiff(names(d),"target")])
  z <- sweep(sweep(x,2,a$initial_scale$sby_center,`-`),2,a$initial_scale$sby_scale,`/`)
  mi <- which(d$target=="rare")
  hits <- vapply(mi,function(i) {
    dd <- sqrt(rowSums(sweep(z,2,z[i,],`-`)^2)); dd[i] <- Inf
    sum(d$target[order(dd,seq_along(dd))[1:3]]=="major")
  },integer(1))
  expect_identical(a$adasyn$difficulty_hits,hits)
  exact <- 5*hits/sum(hits); quota <- floor(exact)
  remaining <- 5-sum(quota)
  if(remaining>0) quota[order(-(exact-quota),mi)[seq_len(remaining)]] <- quota[order(-(exact-quota),mi)[seq_len(remaining)]]+1
  expect_equal(a$adasyn$quotas,quota)
})

test_that("stable inverse keeps tiny synthetics inside original rare convex hull", {
  d <- data.frame(x=c(0,1e-8,100,101,102),target=factor(c("r","r","m","m","m")))
  out <- sby_adasyn_hpc(d,target~.,sby_adasyn_ratio=2,sby_adasyn_k=3L,
    sby_seed=2L,sby_audit=TRUE)
  a <- attr(out,"audit")
  expect_true(all(a$synthetic_continuous>=0 & a$synthetic_continuous<=1e-8))
  expect_equal(a$initial_scale$sby_center,mean(d$x),ignore_attr=TRUE)
  expect_equal(a$initial_scale$sby_scale,sqrt(mean((d$x-mean(d$x))^2)),ignore_attr=TRUE)
})

test_that("beta, d_th and zero difficulty policies are explicit", {
  d <- fixture()
  out <- sby_adasyn_hpc(d,target~.,sby_adasyn_beta=.5,sby_seed=3L,sby_audit=TRUE)
  expect_length(attr(out,"sbyaudit")$adasyn,3L)
  expect_error(sby_adasyn_hpc(d,target~.,sby_adasyn_beta=.5,sby_adasyn_ratio=.2),"OU")
  off <- sby_adasyn_hpc(d,target~.,sby_adasyn_ratio=2,sby_adasyn_d_th=.4,sby_seed=1L)
  expect_length(attr(off,"sbyaudit")$adasyn,0L)
  separated <- data.frame(x=c(0,.1,.2,100,101,102,103),target=factor(c(rep("r",3),rep("m",4))))
  expect_error(sby_adasyn_hpc(separated,target~.,sby_adasyn_ratio=1,sby_adasyn_k=1L,sby_seed=2L),"zero")
  uniform <- sby_adasyn_hpc(separated,target~.,sby_adasyn_ratio=1,sby_adasyn_k=1L,
    sby_seed=2L,sby_adasyn_zero_difficulty="uniform",sby_audit=TRUE)
  expect_true(attr(uniform,"audit")$adasyn$uniform_fallback)
  expect_identical(attr(uniform,"audit")$adasyn$quotas,rep(1L,3L))
})

test_that("NearMiss-3 shortage and retention overflow have deliberate outcomes", {
  d <- data.frame(x=c(0,.1,.01,10,11,12,13,14),target=factor(c("r","r",rep("m",6))))
  expect_warning(out <- sby_nearmiss_hpc(d,target~.,sby_nearmiss_ratio=2,
    nearmiss_model=3L,sby_nearmiss_m=1L,sby_seed=1L,sby_audit=TRUE),"insuficientes")
  expect_true(attr(out,"audit")$nearmiss$candidate_shortage)
  expect_length(attr(out,"sbyaudit")$original_majority_indices,1L)
  huge <- sby_nearmiss_hpc(d,target~.,sby_nearmiss_ratio=1e300,sby_seed=1L)
  expect_identical(values(huge),values(d))
  overflow <- sby_nearmiss_hpc(d,target~.,sby_nearmiss_ratio=.Machine$double.xmax,sby_seed=1L)
  expect_identical(values(overflow),values(d))
  expect_error(sby_nearmiss_hpc(d,target~.,sby_nearmiss_ratio=.01,sby_seed=1L),"zero")
})

test_that("attributes describe positions, origins, parameters and honest telemetry", {
  d <- fixture(); out <- sampler(d,sby_adasyn_ratio=.75,nearmiss_model=1L,sby_audit=TRUE)
  s <- attr(out,"sbyaudit"); a <- attr(out,"audit")
  expect_identical(s$package,"sbyadanear")
  expect_identical(s$function_name,"sby_adanear_hpc")
  expect_equal(s$sby_adanear_hpc$sby_seed,17L)
  expect_equal(s$sby_adanear_hpc$nearmiss_model,1L)
  expect_identical(s$adasyn,as.integer(a$synthetics$output_row))
  expect_true(all(out$target[s$adasyn]=="rare"))
  expect_true(all(out$target[s$nearmiss]=="major"))
  expect_identical(values(out[s$nearmiss,]),values(d[s$original_majority_indices,]))
  expect_s3_class(a$telemetry,"tbl_df")
  expect_equal(nrow(a$telemetry),7L)
  expect_true(all(a$telemetry$elapsed_seconds>=0))
  expect_true(all(is.na(a$telemetry$mkl_observed_threads)))
  plain <- sampler(d,sby_adasyn_ratio=.75,nearmiss_model=1L)
  expect_null(attr(plain,"audit"))
  expect_identical(values(plain),values(out))
})

test_that("no-op paths, target types, formula selection and integer-only matrices preserve originals", {
  d <- fixture(); d$target <- as.character(d$target)
  out <- sby_adanear_hpc(d,target~x+integer,sby_adasyn_ratio=0,
    sby_nearmiss_ratio=0,sby_seed=8L,sby_audit=TRUE)
  expect_identical(names(out),c("x","integer","target"))
  expect_identical(values(out),values(d[c("x","integer","target")]))
  expect_false(attr(out,"audit")$adasyn$executed)
  expect_false(attr(out,"audit")$nearmiss$executed)
  expect_true(all(is.na(attr(out,"audit")$adasyn$difficulty_hits)))
  expect_error(sby_nearmiss_hpc(d,target~log(x),sby_seed=2L),"transformações")
  ints <- data.frame(x=as.integer(c(0,1,0,1,0)),target=c(1L,1L,0L,0L,0L))
  no_op <- sby_adanear_hpc(ints,target~.,sby_adasyn_ratio=0,sby_nearmiss_ratio=0,sby_seed=1L)
  expect_identical(values(no_op),values(ints))
  only_one <- data.frame(x=c(.1,1,2,3),target=c("r","m","m","m"))
  expect_error(sby_adasyn_hpc(only_one,target~.,sby_adasyn_ratio=1,sby_seed=1L),"dois")
})

test_that("invalid parameters fail with useful messages before computation", {
  d <- fixture()
  for(v in list(0L,4L,1.5,NA_integer_,c(1L,2L),"3",Inf,NULL))
    expect_error(sby_nearmiss_hpc(d,target~.,nearmiss_model=v),"nearmiss_model")
  for(v in list(0,-2,1.5,NA_real_))
    expect_error(sby_nearmiss_hpc(d,target~.,sby_config_max_threads=v),"sby_config_max_threads")
  expect_error(sby_adasyn_hpc(d,target~.,sby_adasyn_ratio=1e300,sby_seed=1L),"limite")
  expect_error(sby_adanear_matrix(as.matrix(d[-ncol(d)]),d$target,sby_adasyn_ratio=1,sby_max_output_rows=14),"expansão")
  expect_error(sby_adanear_matrix(as.matrix(d[-ncol(d)]),d$target,sby_max_dense_gb=1e-9),"buffers")
  expect_error(sby_nearmiss(target~.,d,sby_knn_engine="FNN"),"oneAPI")
})
