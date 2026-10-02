library(sbyadanear)
print(sby_hpc_cpu_report())
stopifnot(sby_hpc_cpu_report()$effective >= 2L)
for (threads in c(1L, 2L)) {
  cat("THREAD_PROBE",threads,"\n")
  probe <- sbyadanear:::sby_native("sby_thread_probe_cpp",threads,1024L)
  print(probe)
  stopifnot(probe$omp_observed_threads == threads,
            probe$mkl_configured_threads == threads,probe$checksum == 2048)
}
set.seed(55)
d <- as.data.frame(matrix(rnorm(64000),ncol=64))
d$target <- factor(c(rep("rare",100),rep("major",900)))
cat("PUBLIC_PIPELINE",1L,"\n")
one <- sby_adanear_hpc(d,target~.,sby_seed=42L,sby_config_max_threads=1L,
                      nearmiss_model=1L,sby_audit=TRUE)
cat("PUBLIC_PIPELINE",2L,"\n")
two <- sby_adanear_hpc(d,target~.,sby_seed=42L,sby_config_max_threads=2L,
                      nearmiss_model=1L,sby_audit=TRUE)
# Compare values; audit timing and thread diagnostics necessarily differ.
stopifnot(identical(lapply(one,identity),lapply(two,identity)),
          all(attr(two,"audit")$telemetry$omp_observed_threads[2:5] == 2L))
print(attr(two,"audit")$telemetry)
