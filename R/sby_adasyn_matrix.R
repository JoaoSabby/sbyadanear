#' Aplicar ADASYN diretamente sobre matrix double e factor binario
#'
#' @return Lista leve com `sby_x_matrix`, `sby_y_vector`, razoes, distribuicoes e diagnosticos.
#' @export
sby_adasyn_matrix <- function(
  sby_x_matrix,
  sby_y_vector,
  sby_over_ratio = 0.2,
  sby_knn_over_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_audit_level = c("none", "light", "full"),
  sby_return_scaled = FALSE,
  sby_return_original_scale = TRUE,
  sby_scaling_info = NULL,
  sby_input_already_scaled = FALSE,
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
  sby_adanear_check_user_interrupt()

  sby_audit_level <- sby_resolve_audit_level(sby_audit, sby_audit_level)
  sby_audit <- identical(sby_audit_level, "full")
  sby_return_scaled <- sby_validate_logical_scalar(sby_return_scaled, "sby_return_scaled")
  sby_return_original_scale <- sby_validate_logical_scalar(sby_return_original_scale, "sby_return_original_scale")
  sby_input_already_scaled <- sby_validate_logical_scalar(sby_input_already_scaled, "sby_input_already_scaled")
  sby_memory_guard <- sby_validate_logical_scalar(sby_memory_guard, "sby_memory_guard")

  sby_x_matrix <- sby_validate_dense_double_matrix(sby_x_matrix = sby_x_matrix)
  if(length(sby_y_vector) != nrow(sby_x_matrix)){
    sby_adanear_abort("'sby_y_vector' deve ter comprimento igual a nrow(sby_x_matrix)")
  }
  sby_class_info_input <- sby_binary_class_counts_fast(sby_y_vector)

  sby_validate_seed(sby_seed = sby_seed)
  sby_knn_over_k <- sby_validate_positive_integer_scalar(sby_knn_over_k, "sby_knn_over_k")
  sby_knn_algorithm <- match.arg(sby_knn_algorithm)
  sby_knn_engine <- match.arg(sby_knn_engine)
  sby_knn_distance_metric <- match.arg(sby_knn_distance_metric)
  sby_knn_workers <- sby_validate_knn_workers(sby_knn_workers)
  sby_hnsw_params <- sby_validate_hnsw_params(sby_knn_hnsw_m, sby_knn_hnsw_ef)
  sby_knn_hnsw_m <- sby_hnsw_params$sby_knn_hnsw_m
  sby_knn_hnsw_ef <- sby_hnsw_params$sby_knn_hnsw_ef
  sby_knn_engine <- sby_resolve_knn_engine(sby_knn_engine, sby_knn_workers)
  sby_knn_algorithm <- sby_resolve_knn_algorithm(sby_knn_algorithm, NCOL(sby_x_matrix), sby_knn_engine)

  sby_synthetic_count <- sby_compute_minority_expansion_count(sby_y_vector, sby_over_ratio)
  sby_output_rows <- nrow(sby_x_matrix) + sby_synthetic_count
  if(is.finite(sby_max_output_rows) && sby_output_rows > sby_max_output_rows){
    sby_adanear_abort("'sby_max_output_rows' seria excedido pelo ADASYN")
  }
  if(isTRUE(sby_memory_guard)){
    sby_check_dense_memory_budget(sby_output_rows, ncol(sby_x_matrix), 2L, sby_max_dense_gb, "sby_adasyn_matrix")
  }

  if(isTRUE(sby_input_already_scaled)){
    if(is.null(sby_scaling_info)){
      sby_adanear_abort("'sby_scaling_info' e obrigatorio quando 'sby_input_already_scaled = TRUE'")
    }
    sby_validate_scaling_info(sby_scaling_info, NCOL(sby_x_matrix))
    sby_x_scaled <- sby_x_matrix
  }else{
    if(is.null(sby_scaling_info)){
      sby_scaling_info <- sby_compute_z_score_params(sby_x_matrix)
    }else{
      sby_validate_scaling_info(sby_scaling_info, NCOL(sby_x_matrix))
    }
    sby_x_scaled <- sby_apply_z_score_scaling_matrix(sby_x_matrix, sby_scaling_info)
  }

  set.seed(sby_seed)
  sby_adasyn_result <- sby_generate_adasyn_samples(
    sby_x_scaled = sby_x_scaled,
    sby_target_factor = sby_y_vector,
    sby_synthetic_count = sby_synthetic_count,
    sby_knn_over_k = sby_knn_over_k,
    sby_knn_algorithm = sby_knn_algorithm,
    sby_knn_engine = sby_knn_engine,
    sby_knn_distance_metric = sby_knn_distance_metric,
    sby_knn_workers = sby_knn_workers,
    sby_knn_hnsw_m = sby_knn_hnsw_m,
    sby_knn_hnsw_ef = sby_knn_hnsw_ef
  )
  colnames(sby_adasyn_result$x) <- colnames(sby_x_matrix)

  if(isTRUE(sby_return_original_scale)){
    sby_x_out <- sby_revert_z_score_scaling_matrix(sby_adasyn_result$x, sby_scaling_info)
    sby_output_scale <- "original"
  }else{
    sby_x_out <- sby_adasyn_result$x
    sby_output_scale <- "z_score"
  }
  sby_y_out <- as.factor(sby_adasyn_result$y)
  sby_class_info_output <- sby_binary_class_counts_fast(sby_y_out)

  sby_diagnostics <- list(
    sby_method = "adasyn",
    sby_input_rows = nrow(sby_x_matrix),
    sby_output_rows = nrow(sby_x_out),
    sby_generated_rows = sby_synthetic_count,
    sby_output_scale = sby_output_scale,
    sby_knn_engine = sby_knn_engine,
    sby_knn_algorithm = sby_knn_algorithm,
    sby_knn_distance_metric = sby_knn_distance_metric,
    sby_knn_workers = sby_knn_workers,
    sby_minority_label = sby_class_info_input$sby_minority_label,
    sby_majority_label = sby_class_info_input$sby_majority_label
  )

  sby_result <- list(
    sby_x_matrix = sby_x_out,
    sby_y_vector = sby_y_out,
    sby_class_ratio_input = sby_class_info_input$sby_class_ratio,
    sby_class_ratio_output = sby_class_info_output$sby_class_ratio,
    sby_input_class_distribution = sby_class_info_input$sby_class_counts,
    sby_output_class_distribution = sby_class_info_output$sby_class_counts,
    sby_diagnostics = sby_diagnostics
  )

  if(isTRUE(sby_audit) || isTRUE(sby_return_scaled)){
    sby_result$sby_scaling_info <- sby_scaling_info
  }
  if(isTRUE(sby_return_scaled)){
    sby_result$sby_balanced_scaled <- list(x = sby_adasyn_result$x, y = sby_y_out)
  }

  return(sby_result)
}
