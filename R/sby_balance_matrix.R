#' Balancear matrix double e factor binario por estrategia industrial
#'
#' @param sby_knn_query_chunk_size Número inteiro positivo que define quantas linhas de consulta KNN são processadas por bloco. O padrão é `1000L`; ajuste para equilibrar overhead de chamadas e pico de memória.
#'
#' @return Lista leve com matrix, alvo, razoes de classe e diagnosticos.
#'
#' @export
sby_balance_matrix <- function(
  sby_x_matrix,
  sby_y_vector,
  sby_strategy = c("none", "weight", "adasyn", "nearmiss", "adanear", "adanearWeight"),
  sby_adasyn_ratio = 0.2,
  sby_nearmiss_ratio = 1,
  sby_knn_over_k = 5L,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_audit_level = c("none", "light", "full"),
  sby_return_index = TRUE,
  sby_return_scaled = FALSE,
  sby_return_original_scale = TRUE,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "native", "FNN", "RcppHNSW", "KernelKnn", "bigKNN"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L,
  sby_knn_parallel_backend = c("parallel", "RcppParallel"),
  sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L,
  sby_knn_query_chunk_size = 1000L,
  sby_memory_guard = TRUE,
  sby_max_output_rows = Inf,
  sby_max_dense_gb = Inf
){
  sby_audit_level <- sby_resolve_audit_level(sby_audit, sby_audit_level)
  sby_audit_full <- identical(sby_audit_level, "full")
  sby_strategy <- match.arg(sby_strategy)

  if(sby_strategy %in% c("none", "weight")){
    sby_validate_matrix_like_row_count(sby_x_matrix, sby_y_vector)
    sby_class_info_input <- sby_binary_class_counts_fast(sby_y_vector)
    return(list(
      sby_x_matrix = sby_x_matrix,
      sby_y_vector = sby_y_vector,
      sby_class_ratio_input = sby_class_info_input$sby_class_ratio,
      sby_class_ratio_output = sby_class_info_input$sby_class_ratio,
      sby_input_class_distribution = sby_class_info_input$sby_class_counts,
      sby_output_class_distribution = sby_class_info_input$sby_class_counts,
      sby_diagnostics = list(
        sby_method = sby_strategy,
        sby_input_rows = collapse::fnrow(sby_x_matrix),
        sby_output_rows = collapse::fnrow(sby_x_matrix),
        sby_data_changed = FALSE
      )
    ))
  }

  sby_x_matrix <- sby_validate_dense_double_matrix(sby_x_matrix = sby_x_matrix)
  if(length(sby_y_vector) != collapse::fnrow(sby_x_matrix)){
    sby_adanear_abort("'sby_y_vector' deve ter comprimento igual ao numero de linhas de 'sby_x_matrix'")
  }
  sby_class_info_input <- sby_binary_class_counts_fast(sby_y_vector)

  if(identical(sby_strategy, "adasyn")){
    return(sby_adasyn_matrix(
      sby_x_matrix = sby_x_matrix,
      sby_y_vector = sby_y_vector,
      sby_adasyn_ratio = sby_adasyn_ratio,
      sby_knn_over_k = sby_knn_over_k,
      sby_seed = sby_seed,
      sby_audit = sby_audit_full,
      sby_audit_level = sby_audit_level,
      sby_return_scaled = sby_return_scaled,
      sby_return_original_scale = sby_return_original_scale,
      sby_knn_algorithm = sby_knn_algorithm,
      sby_knn_engine = sby_knn_engine,
      sby_knn_distance_metric = sby_knn_distance_metric,
      sby_knn_workers = sby_knn_workers,
      sby_knn_parallel_backend = sby_knn_parallel_backend,
      sby_knn_hnsw_m = sby_knn_hnsw_m,
      sby_knn_hnsw_ef = sby_knn_hnsw_ef,
      sby_knn_query_chunk_size = sby_knn_query_chunk_size,
      sby_memory_guard = sby_memory_guard,
      sby_max_output_rows = sby_max_output_rows,
      sby_max_dense_gb = sby_max_dense_gb
    ))
  }

  if(identical(sby_strategy, "nearmiss")){
    return(sby_nearmiss_matrix(
      sby_x_matrix = sby_x_matrix,
      sby_y_vector = sby_y_vector,
      sby_nearmiss_ratio = sby_nearmiss_ratio,
      sby_knn_under_k = sby_knn_under_k,
      sby_seed = sby_seed,
      sby_audit = sby_audit_full,
      sby_audit_level = sby_audit_level,
      sby_return_index = sby_return_index,
      sby_return_scaled = sby_return_scaled,
      sby_return_original_scale = sby_return_original_scale,
      sby_knn_algorithm = sby_knn_algorithm,
      sby_knn_engine = sby_knn_engine,
      sby_knn_distance_metric = sby_knn_distance_metric,
      sby_knn_workers = sby_knn_workers,
      sby_knn_parallel_backend = sby_knn_parallel_backend,
      sby_knn_hnsw_m = sby_knn_hnsw_m,
      sby_knn_hnsw_ef = sby_knn_hnsw_ef,
      sby_knn_query_chunk_size = sby_knn_query_chunk_size,
      sby_memory_guard = sby_memory_guard
    ))
  }

  sby_result <- sby_adanear_matrix(
    sby_x_matrix = sby_x_matrix,
    sby_y_vector = sby_y_vector,
    sby_adasyn_ratio = sby_adasyn_ratio,
    sby_nearmiss_ratio = sby_nearmiss_ratio,
    sby_knn_over_k = sby_knn_over_k,
    sby_knn_under_k = sby_knn_under_k,
    sby_seed = sby_seed,
    sby_audit = sby_audit_full,
    sby_audit_level = sby_audit_level,
    sby_return_scaled = sby_return_scaled,
    sby_return_original_scale = sby_return_original_scale,
    sby_knn_algorithm = sby_knn_algorithm,
    sby_knn_engine = sby_knn_engine,
    sby_knn_distance_metric = sby_knn_distance_metric,
    sby_knn_workers = sby_knn_workers,
    sby_knn_parallel_backend = sby_knn_parallel_backend,
    sby_knn_hnsw_m = sby_knn_hnsw_m,
    sby_knn_hnsw_ef = sby_knn_hnsw_ef,
    sby_knn_query_chunk_size = sby_knn_query_chunk_size,
    sby_memory_guard = sby_memory_guard,
    sby_max_output_rows = sby_max_output_rows,
    sby_max_dense_gb = sby_max_dense_gb
  )
  if(identical(sby_strategy, "adanearWeight")){
    sby_result$sby_diagnostics$sby_method <- "adanearWeight"
    sby_result$sby_diagnostics$sby_weight_guidance <- "Use sby_class_ratio_output para scale_pos_weight pos-balanceamento."
  }
  return(sby_result)
}
