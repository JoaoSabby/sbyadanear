#' Aplicar ADANEAR diretamente sobre matrix double e factor binario
#'
#' @return Lista leve com `sby_x_matrix`, `sby_y_vector`, razoes, distribuicoes e diagnosticos.
#' @export
sby_adanear_matrix <- function(
  sby_x_matrix,
  sby_y_vector,
  sby_over_ratio = 0.2,
  sby_under_ratio = 0.5,
  sby_knn_over_k = 5L,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_audit_level = c("none", "light", "full"),
  sby_return_scaled = FALSE,
  sby_return_original_scale = TRUE,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "FNN", "RcppHNSW"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L,
  sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L,
  sby_memory_guard = TRUE,
  sby_max_output_rows = Inf,
  sby_max_dense_gb = Inf
){
  sby_audit_level <- sby_resolve_audit_level(sby_audit, sby_audit_level)
  sby_audit <- identical(sby_audit_level, "full")
  sby_return_scaled <- sby_validate_logical_scalar(sby_return_scaled, "sby_return_scaled")
  sby_return_original_scale <- sby_validate_logical_scalar(sby_return_original_scale, "sby_return_original_scale")
  sby_x_matrix <- sby_validate_dense_double_matrix(sby_x_matrix = sby_x_matrix)
  if(length(sby_y_vector) != nrow(sby_x_matrix)){
    sby_adanear_abort("'sby_y_vector' deve ter comprimento igual a nrow(sby_x_matrix)")
  }
  sby_class_info_input <- sby_binary_class_counts_fast(sby_y_vector)
  sby_original_roles <- sby_get_binary_class_roles(sby_target_factor = sby_y_vector)

  sby_over_result <- sby_adasyn_matrix(
    sby_x_matrix = sby_x_matrix,
    sby_y_vector = sby_y_vector,
    sby_over_ratio = sby_over_ratio,
    sby_knn_over_k = sby_knn_over_k,
    sby_seed = sby_seed,
    sby_audit = TRUE,
    sby_return_scaled = TRUE,
    sby_return_original_scale = FALSE,
    sby_knn_algorithm = sby_knn_algorithm,
    sby_knn_engine = sby_knn_engine,
    sby_knn_distance_metric = sby_knn_distance_metric,
    sby_knn_workers = sby_knn_workers,
    sby_knn_hnsw_m = sby_knn_hnsw_m,
    sby_knn_hnsw_ef = sby_knn_hnsw_ef,
    sby_memory_guard = sby_memory_guard,
    sby_max_output_rows = sby_max_output_rows,
    sby_max_dense_gb = sby_max_dense_gb
  )

  sby_under_result <- sby_nearmiss_matrix(
    sby_x_matrix = sby_over_result$sby_balanced_scaled$x,
    sby_y_vector = sby_over_result$sby_balanced_scaled$y,
    sby_under_ratio = sby_under_ratio,
    sby_knn_under_k = sby_knn_under_k,
    sby_seed = sby_seed,
    sby_audit = TRUE,
    sby_return_index = TRUE,
    sby_return_scaled = TRUE,
    sby_return_original_scale = sby_return_original_scale,
    sby_scaling_info = sby_over_result$sby_scaling_info,
    sby_input_already_scaled = TRUE,
    sby_fixed_minority_label = sby_original_roles$sby_minority_label,
    sby_fixed_majority_label = sby_original_roles$sby_majority_label,
    sby_knn_algorithm = sby_knn_algorithm,
    sby_knn_engine = sby_knn_engine,
    sby_knn_distance_metric = sby_knn_distance_metric,
    sby_knn_workers = sby_knn_workers,
    sby_knn_hnsw_m = sby_knn_hnsw_m,
    sby_knn_hnsw_ef = sby_knn_hnsw_ef,
    sby_memory_guard = sby_memory_guard
  )

  sby_class_info_output <- sby_binary_class_counts_fast(sby_under_result$sby_y_vector)
  sby_diagnostics <- list(
    sby_method = "adanear",
    sby_input_rows = nrow(sby_x_matrix),
    sby_after_oversampling_rows = nrow(sby_over_result$sby_balanced_scaled$x),
    sby_output_rows = nrow(sby_under_result$sby_x_matrix),
    sby_output_scale = sby_under_result$sby_diagnostics$sby_output_scale,
    sby_original_minority_label = sby_original_roles$sby_minority_label,
    sby_original_majority_label = sby_original_roles$sby_majority_label,
    sby_oversampling_diagnostics = sby_over_result$sby_diagnostics,
    sby_undersampling_diagnostics = sby_under_result$sby_diagnostics
  )

  sby_result <- list(
    sby_x_matrix = sby_under_result$sby_x_matrix,
    sby_y_vector = sby_under_result$sby_y_vector,
    sby_class_ratio_input = sby_class_info_input$sby_class_ratio,
    sby_class_ratio_output = sby_class_info_output$sby_class_ratio,
    sby_input_class_distribution = sby_class_info_input$sby_class_counts,
    sby_output_class_distribution = sby_class_info_output$sby_class_counts,
    sby_diagnostics = sby_diagnostics
  )

  if(isTRUE(sby_audit)){
    sby_result$sby_oversampling_result <- sby_over_result
    sby_result$sby_undersampling_result <- sby_under_result
    sby_result$sby_scaling_info <- sby_over_result$sby_scaling_info
    sby_result$sby_retained_index <- sby_under_result$sby_retained_index
  }
  if(isTRUE(sby_return_scaled)){
    sby_result$sby_balanced_scaled <- sby_under_result$sby_balanced_scaled
    sby_result$sby_scaling_info <- sby_over_result$sby_scaling_info
  }

  return(sby_result)
}
