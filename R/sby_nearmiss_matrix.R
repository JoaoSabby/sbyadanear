#' Aplicar NearMiss-1 diretamente sobre matrix double e factor binario
#'
#' @param sby_knn_query_chunk_size Número inteiro positivo que define quantas linhas de consulta KNN são processadas por bloco. O padrão é `1000L`; ajuste para equilibrar overhead de chamadas e pico de memória.
#'
#' @return Lista leve com `sby_x_matrix`, `sby_y_vector`, razoes, distribuicoes e diagnosticos.
#' @export
sby_nearmiss_matrix <- function(
  sby_x_matrix,
  sby_y_vector,
  sby_under_ratio = 1,
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
  sby_knn_engine = c("auto", "native", "FNN", "RcppHNSW", "KernelKnn", "bigKNN"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L,
  sby_knn_parallel_backend = c("parallel", "RcppParallel"),
  sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L,
  sby_knn_query_chunk_size = 1000L,
  sby_memory_guard = TRUE
){
  sby_audit_level <- sby_resolve_audit_level(sby_audit, sby_audit_level)
  sby_audit_full <- identical(sby_audit_level, "full")
  sby_audit_light <- sby_audit_level %in% c("light", "full")
  sby_return_index <- sby_validate_logical_scalar(sby_return_index, "sby_return_index")
  sby_return_scaled <- sby_validate_logical_scalar(sby_return_scaled, "sby_return_scaled")
  sby_return_original_scale <- sby_validate_logical_scalar(sby_return_original_scale, "sby_return_original_scale")
  sby_memory_guard <- sby_validate_logical_scalar(sby_memory_guard, "sby_memory_guard")

  sby_x_matrix <- sby_validate_dense_double_matrix(sby_x_matrix = sby_x_matrix)
  if(length(sby_y_vector) != collapse::fnrow(sby_x_matrix)){
    sby_adanear_abort("'sby_y_vector' deve ter comprimento igual ao numero de linhas de 'sby_x_matrix'")
  }
  sby_class_info_input <- sby_binary_class_counts_fast(sby_y_vector)
  if(isTRUE(sby_memory_guard)){
    sby_check_dense_memory_budget(collapse::fnrow(sby_x_matrix), collapse::fncol(sby_x_matrix), 2L, Inf, "sby_nearmiss_matrix")
  }

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
    sby_knn_parallel_backend = sby_knn_parallel_backend,
    sby_knn_hnsw_m = sby_knn_hnsw_m,
    sby_knn_hnsw_ef = sby_knn_hnsw_ef,
    sby_knn_query_chunk_size = sby_knn_query_chunk_size,
    sby_audit = sby_audit_full,
    sby_audit_level = sby_audit_level,
    sby_return_scaling_info = TRUE,
    sby_return_reduced_scaled = TRUE
  )
  sby_retained_index <- sby_index_result$sby_retained_index
  sby_scaling_info <- sby_index_result$sby_scaling_info

  sby_reduced_scaled <- sby_index_result$sby_reduced_scaled

  if(isTRUE(sby_return_original_scale)){
    sby_x_out <- sby_revert_z_score_scaling_matrix(sby_reduced_scaled, sby_scaling_info, sby_engine = sby_knn_engine)
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
  if(isTRUE(sby_audit_light)){
    sby_result$sby_diagnostics$sby_audit_level <- sby_audit_level
  }
  if(isTRUE(sby_return_index) || isTRUE(sby_audit_full)){
    sby_result$sby_retained_index <- sby_retained_index
  }
  if(isTRUE(sby_audit_full) || isTRUE(sby_return_scaled)){
    sby_result$sby_scaling_info <- sby_scaling_info
  }
  if(isTRUE(sby_return_scaled)){
    sby_result$sby_balanced_scaled <- list(x = sby_reduced_scaled, y = sby_y_out)
  }
  return(sby_result)
}
