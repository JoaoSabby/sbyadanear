# Native symbol bindings for registered .Call interfaces.
# The package disables dynamic symbol lookup, so .Call must receive these
# NativeSymbolInfo objects instead of character names.

sby_load_native_symbol <- function(sby_symbol_name){
  getNativeSymbolInfo(
    name = sby_symbol_name,
    PACKAGE = "sbyadanear"
  )
}

sby_native_symbol_available <- function(sby_symbol_name){
  sby_namespace <- asNamespace("sbyadanear")
  if(!exists(sby_symbol_name, envir = sby_namespace, inherits = FALSE)){
    return(FALSE)
  }

  sby_symbol <- get(
    x = sby_symbol_name,
    envir = sby_namespace,
    inherits = FALSE
  )

  return(inherits(sby_symbol, "NativeSymbolInfo"))
}

check_user_interrupt_c <- sby_load_native_symbol("check_user_interrupt_c")
compute_z_score_params_c <- sby_load_native_symbol("compute_z_score_params_c")
apply_z_score_c <- sby_load_native_symbol("apply_z_score_c")
generate_synthetic_adasyn_c <- sby_load_native_symbol("generate_synthetic_adasyn_c")
generate_synthetic_adasyn_col_c <- sby_load_native_symbol("generate_synthetic_adasyn_col_c")
select_nearmiss_majority_c <- sby_load_native_symbol("select_nearmiss_majority_c")
drop_self_neighbor_c <- sby_load_native_symbol("drop_self_neighbor_c")
brute_force_knn_c <- sby_load_native_symbol("brute_force_knn_c")
brute_force_knn_index_c <- sby_load_native_symbol("brute_force_knn_index_c")
brute_force_knn_dist_c <- sby_load_native_symbol("brute_force_knn_dist_c")
brute_force_knn_rcpp_parallel_c <- sby_load_native_symbol("brute_force_knn_rcpp_parallel_c")
brute_force_knn_native_c <- sby_load_native_symbol("brute_force_knn_native_c")
brute_force_knn_native_parallel_c <- sby_load_native_symbol("brute_force_knn_native_parallel_c")
rcpp_parallel_uses_tbb_c <- sby_load_native_symbol("rcpp_parallel_uses_tbb_c")
nearmiss_brute_select_c <- sby_load_native_symbol("nearmiss_brute_select_c")
rbind_double_matrix_c <- sby_load_native_symbol("rbind_double_matrix_c")
compute_zscore_population_fortran_c <- sby_load_native_symbol("compute_zscore_population_fortran_c")
apply_zscore_fortran_c <- sby_load_native_symbol("apply_zscore_fortran_c")
revert_zscore_fortran_c <- sby_load_native_symbol("revert_zscore_fortran_c")
rbind_matrix_fortran_c <- sby_load_native_symbol("rbind_matrix_fortran_c")
