#' Selecionar indices NearMiss-1 diretamente sobre matrix double e factor binario
#'
#' @return Lista com indices retidos, distribuicoes e diagnosticos.
#' @export
sby_nearmiss_index <- function(
  sby_x_matrix,
  sby_y_vector,
  sby_under_ratio = 0.5,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
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
  sby_audit = FALSE,
  sby_audit_level = c("none", "light", "full")
){
  sby_adanear_check_user_interrupt()

  sby_audit_level <- sby_resolve_audit_level(sby_audit, sby_audit_level)
  sby_audit <- identical(sby_audit_level, "full")
  sby_input_already_scaled <- sby_validate_logical_scalar(sby_input_already_scaled, "sby_input_already_scaled")
  sby_x_matrix <- sby_validate_dense_double_matrix(sby_x_matrix = sby_x_matrix)
  if(length(sby_y_vector) != nrow(sby_x_matrix)){
    sby_adanear_abort("'sby_y_vector' deve ter comprimento igual a nrow(sby_x_matrix)")
  }
  sby_class_info_input <- sby_binary_class_counts_fast(sby_y_vector)

  sby_validate_seed(sby_seed = sby_seed)
  sby_knn_under_k <- sby_validate_positive_integer_scalar(sby_knn_under_k, "sby_knn_under_k")
  sby_knn_algorithm <- match.arg(sby_knn_algorithm)
  sby_knn_engine <- match.arg(sby_knn_engine)
  sby_knn_distance_metric <- match.arg(sby_knn_distance_metric)
  sby_knn_workers <- sby_validate_knn_workers(sby_knn_workers)
  sby_hnsw_params <- sby_validate_hnsw_params(sby_knn_hnsw_m, sby_knn_hnsw_ef)
  sby_knn_hnsw_m <- sby_hnsw_params$sby_knn_hnsw_m
  sby_knn_hnsw_ef <- sby_hnsw_params$sby_knn_hnsw_ef
  sby_knn_engine <- sby_resolve_knn_engine(sby_knn_engine, sby_knn_workers)
  sby_knn_algorithm <- sby_resolve_knn_algorithm(sby_knn_algorithm, NCOL(sby_x_matrix), sby_knn_engine)

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

  if(identical(sby_class_info_input$sby_minority_count, sby_class_info_input$sby_majority_count)){
    sby_retained_index <- seq_len(nrow(sby_x_matrix))
    sby_selected_majority_index <- integer(0L)
  }else{
    sby_class_roles <- sby_get_binary_class_roles(
      sby_target_factor = sby_y_vector,
      sby_minority_label = sby_fixed_minority_label,
      sby_majority_label = sby_fixed_majority_label
    )
    sby_minority_index <- which(sby_y_vector == sby_class_roles$sby_minority_label)
    sby_majority_index <- which(sby_y_vector == sby_class_roles$sby_majority_label)
    sby_minority_matrix <- sby_x_scaled[sby_minority_index, , drop = FALSE]
    sby_majority_matrix <- sby_x_scaled[sby_majority_index, , drop = FALSE]

    sby_retained_majority_count <- sby_compute_majority_retention_count(
      sby_target_factor = sby_y_vector,
      sby_under_ratio = sby_under_ratio,
      sby_minority_label = sby_class_roles$sby_minority_label,
      sby_majority_label = sby_class_roles$sby_majority_label
    )
    sby_effective_k <- min(as.integer(sby_knn_under_k), nrow(sby_minority_matrix))
    if(sby_effective_k < 1L){
      sby_adanear_abort("Sem linhas minoritarias suficientes para NearMiss")
    }

    set.seed(sby_seed)
    sby_knn_result <- sby_get_knnx(
      sby_data = sby_minority_matrix,
      sby_query = sby_majority_matrix,
      sby_k = sby_effective_k,
      sby_knn_algorithm = sby_knn_algorithm,
      sby_knn_engine = sby_knn_engine,
      sby_knn_distance_metric = sby_knn_distance_metric,
      sby_knn_workers = sby_knn_workers,
      sby_knn_hnsw_m = sby_knn_hnsw_m,
      sby_knn_hnsw_ef = sby_knn_hnsw_ef,
      sby_knn_return = "dist"
    )

    if(sby_adanear_native_available()){
      sby_selected_majority_index <- .Call(
        OU_SelectNearMissMajorityC,
        sby_knn_result$nn.dist,
        as.integer(sby_majority_index),
        as.integer(sby_retained_majority_count)
      )
    }else{
      sby_mean_distances <- rowMeans(sby_knn_result$nn.dist)
      sby_selected_order <- order(sby_mean_distances, decreasing = FALSE)
      sby_selected_majority_index <- sby_majority_index[sby_selected_order[seq_len(sby_retained_majority_count)]]
    }
    sby_retained_index <- sort(c(sby_minority_index, sby_selected_majority_index))
  }

  sby_y_out <- sby_y_vector[sby_retained_index]
  sby_class_info_output <- sby_binary_class_counts_fast(sby_y_out)
  sby_diagnostics <- list(
    sby_method = "nearmiss_index",
    sby_input_rows = nrow(sby_x_matrix),
    sby_output_rows = length(sby_retained_index),
    sby_removed_rows = nrow(sby_x_matrix) - length(sby_retained_index),
    sby_knn_engine = sby_knn_engine,
    sby_knn_algorithm = sby_knn_algorithm,
    sby_knn_distance_metric = sby_knn_distance_metric,
    sby_knn_workers = sby_knn_workers
  )

  sby_result <- list(
    sby_retained_index = as.integer(sby_retained_index),
    sby_input_class_distribution = sby_class_info_input$sby_class_counts,
    sby_output_class_distribution = sby_class_info_output$sby_class_counts,
    sby_diagnostics = sby_diagnostics
  )
  if(isTRUE(sby_audit)){
    sby_result$sby_selected_majority_index <- as.integer(sby_selected_majority_index)
    sby_result$sby_scaling_info <- sby_scaling_info
  }
  return(sby_result)
}
