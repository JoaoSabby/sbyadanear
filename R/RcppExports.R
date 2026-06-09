# Native symbol bindings for registered .Call interfaces.
# Symbols are resolved lazily because package installation byte-compiles R code
# before every registered native routine is available through the loaded DLL.

sby_native_symbol_cache <- new.env(parent = emptyenv())

sby_load_native_symbol <- function(sby_symbol_name){
  if(exists(sby_symbol_name, envir = sby_native_symbol_cache, inherits = FALSE)){
    return(get(
      x = sby_symbol_name,
      envir = sby_native_symbol_cache,
      inherits = FALSE
    ))
  }

  sby_symbol <- getNativeSymbolInfo(
    name = sby_symbol_name,
    PACKAGE = "sbyadanear"
  )

  assign(
    x = sby_symbol_name,
    value = sby_symbol,
    envir = sby_native_symbol_cache
  )

  return(sby_symbol)
}

sby_native_symbol_available <- function(sby_symbol_name){
  tryCatch(
    {
      sby_load_native_symbol(sby_symbol_name)
      TRUE
    },
    error = function(sby_error){
      FALSE
    }
  )
}

delayedAssign(check_user_interrupt_c, sby_load_native_symbol("check_user_interrupt_c"))
delayedAssign(compute_z_score_params_c, sby_load_native_symbol("compute_z_score_params_c"))
delayedAssign(apply_z_score_c, sby_load_native_symbol("apply_z_score_c"))
delayedAssign(generate_synthetic_adasyn_c, sby_load_native_symbol("generate_synthetic_adasyn_c"))
delayedAssign(generate_synthetic_adasyn_col_c, sby_load_native_symbol("generate_synthetic_adasyn_col_c"))
delayedAssign(select_nearmiss_majority_c, sby_load_native_symbol("select_nearmiss_majority_c"))
delayedAssign(drop_self_neighbor_c, sby_load_native_symbol("drop_self_neighbor_c"))
delayedAssign(brute_force_knn_c, sby_load_native_symbol("brute_force_knn_c"))
delayedAssign(brute_force_knn_index_c, sby_load_native_symbol("brute_force_knn_index_c"))
delayedAssign(brute_force_knn_dist_c, sby_load_native_symbol("brute_force_knn_dist_c"))
delayedAssign(brute_force_knn_rcpp_parallel_c, sby_load_native_symbol("brute_force_knn_rcpp_parallel_c"))
delayedAssign(brute_force_knn_native_c, sby_load_native_symbol("brute_force_knn_native_c"))
delayedAssign(brute_force_knn_native_parallel_c, sby_load_native_symbol("brute_force_knn_native_parallel_c"))
delayedAssign(rcpp_parallel_uses_tbb_c, sby_load_native_symbol("rcpp_parallel_uses_tbb_c"))
delayedAssign(nearmiss_brute_select_c, sby_load_native_symbol("nearmiss_brute_select_c"))
delayedAssign(rbind_double_matrix_c, sby_load_native_symbol("rbind_double_matrix_c"))
delayedAssign(compute_zscore_population_fortran_c, sby_load_native_symbol("compute_zscore_population_fortran_c"))
delayedAssign(apply_zscore_fortran_c, sby_load_native_symbol("apply_zscore_fortran_c"))
delayedAssign(revert_zscore_fortran_c, sby_load_native_symbol("revert_zscore_fortran_c"))
delayedAssign(rbind_matrix_fortran_c, sby_load_native_symbol("rbind_matrix_fortran_c"))
