test_that("numeric class labels are matched without lossy text conversion", {
  d <- fixture()
  d$target <- c(rep(1,4),rep(1+2*.Machine$double.eps,10))
  expect_identical(as.character(d$target[1]),as.character(d$target[5]))
  o <- sampler(d,sby_adasyn_ratio=.5,sby_nearmiss_ratio=1,nearmiss_model=1L,sby_audit=TRUE)
  a <- attr(o,"audit")
  expect_identical(a$original_minority_indices,1:4)
  expect_equal(a$counts$minority,4)
  expect_equal(a$counts$majority,10)
  expect_true(all(o$target[attr(o,"sbyaudit")$adasyn]==1))
  m <- sby_nearmiss_matrix(as.matrix(d[-ncol(d)]),d$target,sby_seed=1L,nearmiss_model=1L)
  expect_identical(unname(m$sby_input_class_distribution),c(4L,10L))
})

test_that("recipe formulas preserve non-syntactic outcome column names", {
  d <- fixture(); names(d)[ncol(d)] <- "class label"
  rec <- recipes::recipe(`class label`~.,data=d)
  rec <- sby_step_adanear_hpc(rec,recipes::all_outcomes(),sby_seed=1L,nearmiss_model=1L,skip=FALSE)
  trained <- recipes::prep(rec,training=d)
  o <- recipes::bake(trained,new_data=d)
  expect_true("class label" %in% names(o))
  expect_s3_class(o,"data.frame")
})

test_that("audit-free and audited calls preserve identical sampling values", {
  d <- fixture()
  for(model in 1:3) {
    a <- suppressWarnings(sampler(d,sby_adasyn_ratio=1,nearmiss_model=model,sby_audit=FALSE))
    b <- suppressWarnings(sampler(d,sby_adasyn_ratio=1,nearmiss_model=model,sby_audit=TRUE))
    expect_identical(values(a),values(b))
    expect_null(attr(a,"audit"))
  }
})

test_that("legacy workers and budget types are validated even when inactive", {
  d <- fixture()
  expect_error(sby_adasyn(target~.,d,sby_seed=1L,sby_config_max_threads=1L,sby_knn_workers=0L),"sby_knn_workers")
  expect_error(sby_adanear_matrix(as.matrix(d[-ncol(d)]),d$target,sby_seed=1L,sby_max_dense_gb=1i),"sby_max_dense_gb")
  expect_error(sby_adanear_matrix(as.matrix(d[-ncol(d)]),d$target,sby_seed=1L,sby_max_output_rows=1i),"sby_max_output_rows")
  expect_error(sby_adasyn_matrix(as.matrix(d[-ncol(d)]),d$target,sby_seed=1L,sby_scaling_info=1),"Scaling information")
})

test_that("finite Euclidean distances survive squared overflow and underflow", {
  x <- matrix(c(0,1,10,9,.1,.75),ncol=1,dimnames=list(NULL,"x"))
  y <- c("rare","rare",rep("major",4))
  for(sd in c(1e-200,1e200)) {
    o <- sby_nearmiss_matrix(x,y,sby_scaling_info=list(sby_center=0,sby_scale=sd),
      sby_nearmiss_ratio=.5,sby_knn_under_k=1L,nearmiss_model=1L,sby_seed=1L,sby_audit=TRUE)
    expect_identical(attr(o,"sbyaudit")$original_majority_indices,5L)
    expect_true(all(is.finite(attr(o,"audit")$nearmiss$scores)))
  }
})

test_that("exact neighbor heap handles duplicates and K reaching class size", {
  d <- fixture(); d$x[2] <- d$x[1]; d$z[2] <- d$z[1]
  for(model in 1:3) for(k in c(1L,4L,100L)) {
    o <- sby_nearmiss_hpc(d,target~.,sby_nearmiss_k=k,sby_nearmiss_m=k,
      sby_nearmiss_ratio=.5,nearmiss_model=model,sby_seed=1L,sby_audit=TRUE)
    a <- attr(o,"audit"); z <- scale(as.matrix(d[-ncol(d)]),
      center=a$initial_scale$sby_center,scale=a$initial_scale$sby_scale)
    chosen <- oracle_nearmiss(z,z[1:4,,drop=FALSE],z[5:14,,drop=FALSE],k,model,k,2L)
    expect_identical(a$original_majority_indices,as.integer(chosen+4L))
  }
})

test_that("thread plans use current cgroup and ancestor CPU quotas", {
  root <- tempfile("cgroup-"); dir.create(root)
  membership <- tempfile("membership-")
  on.exit(unlink(c(root,membership),recursive=TRUE),add=TRUE)
  dir.create(file.path(root,"service")); dir.create(file.path(root,"service","job"))
  writeLines("max 100000",file.path(root,"cpu.max"))
  writeLines("150000 100000",file.path(root,"service","cpu.max"))
  writeLines("200000 100000",file.path(root,"service","job","cpu.max"))
  writeLines("0::/service/job",membership)
  quota <- if(exists("sby_quota",mode="function")) get("sby_quota",mode="function") else getFromNamespace("sby_quota","sbyadanear")
  expect_equal(quota(root,membership),1.5)
  writeLines("50000 100000",file.path(root,"service","job","cpu.max"))
  expect_equal(quota(root,membership),.5)
  unlink(file.path(root,"service","job","cpu.max"))
  expect_equal(quota(root,membership),1.5)
  writeLines("max 100000",file.path(root,"service","cpu.max"))
  expect_true(is.na(quota(root,membership)))
  dir.create(file.path(root,"cpu")); dir.create(file.path(root,"cpu","job"))
  writeLines("250000",file.path(root,"cpu","job","cpu.cfs_quota_us"))
  writeLines("100000",file.path(root,"cpu","job","cpu.cfs_period_us"))
  writeLines("2:cpu,cpuacct:/job",membership)
  expect_equal(quota(root,membership),2.5)
})

test_that("unrepresentable distances and invalid scaling fail instructively", {
  x <- matrix(c(-1e308,1e308,0,-9e307,9e307),ncol=1)
  y <- c("rare","rare",rep("major",3))
  expect_error(sby_nearmiss_matrix(x,y,sby_seed=1L,nearmiss_model=1L,
    sby_scaling_info=list(sby_center=0,sby_scale=1)),"Euclidean distances")
  expect_error(sby_nearmiss_matrix(x,y,sby_seed=1L,
    sby_scaling_info=list(sby_center="zero",sby_scale=1)),"Escala fornecida")
})
