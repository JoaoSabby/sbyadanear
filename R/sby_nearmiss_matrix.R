#' Aplicar NearMiss-1 diretamente sobre matrix double e factor binario
#'
#' @return Lista leve com `sby_x_matrix`, `sby_y_vector`, razoes, distribuicoes e diagnosticos.
#' @export
sby_nearmiss_matrix <- function(
  sby_x_matrix,
  sby_y_vector,
  sby_under_ratio = 0.5,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_audit_level = c("none", "light", "full"),
  sby_return_index = TRUE,
  sby_return_scaled = FALSE,
  sby_return_original_scale = TRUE,
  sby_scaling_info = NULL,
  sby_input_already_scaled = FALSE,
  sby_fixed_minority_label = NULL,
  sby_fixed_majority_label = NULL,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "FNN", "RcppHNSW"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L,
  sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L,
  sby_memory_guard = TRUE
){
  sby_audit_level <- sby_resolve_audit_level(sby_audit, sby_audit_level)
  sby_audit <- identical(sby_audit_level, "full")
  sby_return_index <- sby_validate_logical_scalar(sby_return_index, "sby_return_index")
  sby_return_scaled <- sby_validate_logical_scalar(sby_return_scaled, "sby_return_scaled")
  sby_return_original_scale <- sby_validate_logical_scalar(sby_return_original_scale, "sby_return_original_scale")
  sby_memory_guard <- sby_validate_logical_scalar(sby_memory_guard, "sby_memory_guard")

  sby_x_matrix <- sby_validate_dense_double_matrix(sby_x_matrix = sby_x_matrix)
  if(length(sby_y_vector) != nrow(sby_x_matrix)){
    sby_adanear_abort("'sby_y_vector' deve ter comprimento igual a nrow(sby_x_matrix)")
  }
  sby_class_info_input <- sby_binary_class_counts_fast(sby_y_vector)

  sby_index_result <- sby_nearmiss_index(
    sby_x_matrix = sby_x_matrix,
    sby_y_vector = sby_y_vector,
    sby_under_ratio = sby_under_ratio,
    sby_knn_under_k = sby_knn_under_k,
    sby_seed = sby_seed,
    sby_scaling_info = sby_scaling_info,
    sby_input_already_scaled = sby_input_already_scaled,
    sby_fixed_minority_label = sby_fixed_minority_label,
    sby_fixed_majority_label = sby_fixed_majority_label,
    sby_knn_algorithm = sby_knn_algorithm,
    sby_knn_engine = sby_knn_engine,
    sby_knn_distance_metric = sby_knn_distance_metric,
    sby_knn_workers = sby_knn_workers,
    sby_knn_hnsw_m = sby_knn_hnsw_m,
    sby_knn_hnsw_ef = sby_knn_hnsw_ef,
    sby_audit = TRUE
  )
  sby_retained_index <- sby_index_result$sby_retained_index
  sby_scaling_info <- sby_index_result$sby_scaling_info

  if(isTRUE(sby_input_already_scaled)){
    sby_reduced_scaled <- sby_x_matrix[sby_retained_index, , drop = FALSE]
  }else{
    sby_x_scaled <- sby_apply_z_score_scaling_matrix(sby_x_matrix, sby_scaling_info)
    sby_reduced_scaled <- sby_x_scaled[sby_retained_index, , drop = FALSE]
  }

  if(isTRUE(sby_return_original_scale)){
    sby_x_out <- sby_revert_z_score_scaling_matrix(sby_reduced_scaled, sby_scaling_info)
    sby_output_scale <- "original"
  }else{
    sby_x_out <- sby_reduced_scaled
    sby_output_scale <- "z_score"
  }
  sby_y_out <- sby_y_vector[sby_retained_index]
  sby_class_info_output <- sby_binary_class_counts_fast(sby_y_out)
  sby_diagnostics <- sby_index_result$sby_diagnostics
  sby_diagnostics$sby_method <- "nearmiss"
  sby_diagnostics$sby_output_scale <- sby_output_scale

  sby_result <- list(
    sby_x_matrix = sby_x_out,
    sby_y_vector = sby_y_out,
    sby_class_ratio_input = sby_class_info_input$sby_class_ratio,
    sby_class_ratio_output = sby_class_info_output$sby_class_ratio,
    sby_input_class_distribution = sby_class_info_input$sby_class_counts,
    sby_output_class_distribution = sby_class_info_output$sby_class_counts,
    sby_diagnostics = sby_diagnostics
  )
  if(isTRUE(sby_return_index) || isTRUE(sby_audit)){
    sby_result$sby_retained_index <- sby_retained_index
  }
  if(isTRUE(sby_audit) || isTRUE(sby_return_scaled)){
    sby_result$sby_scaling_info <- sby_scaling_info
  }
  if(isTRUE(sby_return_scaled)){
    sby_result$sby_balanced_scaled <- list(x = sby_reduced_scaled, y = sby_y_out)
  }
  return(sby_result)
}
