#' @export
prep.step_sby_step_adasyn_hpc <- function(x, training, info = NULL, ...){
  sby_prep_step_sampling(x = x, training = training, info = info, sby_step_name = "sby_step_adasyn_hpc()")
}
#' @export
prep.step_sby_step_nearmiss_hpc <- function(x, training, info = NULL, ...){
  sby_prep_step_sampling(x = x, training = training, info = info, sby_step_name = "sby_step_nearmiss_hpc()")
}
#' @export
prep.step_sby_step_adanear_hpc <- function(x, training, info = NULL, ...){
  sby_prep_step_sampling(x = x, training = training, info = info, sby_step_name = "sby_step_adanear_hpc()")
}
#' @export
bake.step_sby_step_adasyn_hpc <- function(object, new_data, ...){
  sby_bake_step_sampling(object = object, new_data = new_data, sby_step_name = "sby_step_adasyn_hpc()")
}
#' @export
bake.step_sby_step_nearmiss_hpc <- function(object, new_data, ...){
  sby_bake_step_sampling(object = object, new_data = new_data, sby_step_name = "sby_step_nearmiss_hpc()")
}
#' @export
bake.step_sby_step_adanear_hpc <- function(object, new_data, ...){
  sby_bake_step_sampling(object = object, new_data = new_data, sby_step_name = "sby_step_adanear_hpc()")
}
#' @export
tidy.step_sby_step_adasyn_hpc <- function(x, ...){sby_tidy_step_sampling(x)}
#' @export
tidy.step_sby_step_nearmiss_hpc <- function(x, ...){sby_tidy_step_sampling(x)}
#' @export
tidy.step_sby_step_adanear_hpc <- function(x, ...){sby_tidy_step_sampling(x)}
#' @export
print.step_sby_step_adasyn_hpc <- function(x, width = max(20, options()$width - 30), ...){sby_print_step_sampling(x, width, "HPC ADASYN sampling")}
#' @export
print.step_sby_step_nearmiss_hpc <- function(x, width = max(20, options()$width - 30), ...){sby_print_step_sampling(x, width, "HPC NearMiss sampling")}
#' @export
print.step_sby_step_adanear_hpc <- function(x, width = max(20, options()$width - 30), ...){sby_print_step_sampling(x, width, "HPC ADANEAR sampling")}
#' @export
required_pkgs.step_sby_step_adasyn_hpc <- function(x, ...){c("sbyadanear", "recipes")}
#' @export
required_pkgs.step_sby_step_nearmiss_hpc <- function(x, ...){c("sbyadanear", "recipes")}
#' @export
required_pkgs.step_sby_step_adanear_hpc <- function(x, ...){c("sbyadanear", "recipes")}
