#' Atalho HPC do oversampling ADASYN
#'
#' @description
#' `sby_adasyn_hpc()` e o atalho de alto desempenho do oversampling ADASYN. Ele
#' gera as linhas sinteticas estritamente no espaco padronizado por z-score,
#' aproveitando a Vector Statistics Library para as estatisticas iniciais, o
#' cblas_sgemm para a matriz de distancias e a interpolacao lambda com
#' vdrnguniform. A reversao final do z-score usa as unidades FMA do hardware
#' durante a montagem zero-copy do tibble no C++. A funcao representa um atalho:
#' a rotina original `sby_adasyn()` continua disponivel.
#'
#' @details
#' Esta funcao substitui internamente a rota denominada "native" como caminho
#' rapido do ADASYN. Quando o motor HPC consolidado nao esta disponivel, recorre
#' de forma transparente a `sby_adasyn()` com `sby_knn_engine = "native"` para
#' garantir que tudo continue funcionando. Apenas `MKL_NUM_THREADS`, `OMP_NUM_THREADS` e
#' `MKL_NUM_STRIPES` sao controladas temporariamente e restauradas por um bloco
#' `on.exit()` inflexivel.
#'
#' @param .data Data frame, tibble ou matriz com a coluna de desfecho e os
#'   preditores numericos referenciados em `formula`.
#' @param formula Formula no formato `alvo ~ preditores`.
#' @param sby_k_neighbor_adanear Numero inteiro positivo de vizinhos usado pela
#'   etapa ADASYN. O padrao e `3`.
#' @param sby_over_ratio Fator relativo de expansao da classe minoritaria. O
#'   padrao e `0.2`.
#' @param sby_config_max_threads Numero inteiro de threads do motor HPC. O padrao
#'   `-1` detecta os nucleos fisicos disponiveis.
#'
#' @return Tibble balanceado com classe `c("tbl_df", "tbl", "data.frame")`.
#' @export
sby_adasyn_hpc <- function(
  .data,
  formula,
  sby_k_neighbor_adanear = 3,
  sby_over_ratio = 0.2,
  sby_config_max_threads = -1
){
  sby_adanear_check_user_interrupt()

  sby_previous_env <- sby_hpc_capture_env()
  on.exit(sby_hpc_restore_env(sby_previous_env), add = TRUE)

  sby_total_threads <- sby_hpc_resolve_threads(sby_config_max_threads)

  sby_formula_data <- sby_extract_formula_data(sby_formula = formula, sby_data = .data)
  sby_predictor_data <- sby_formula_data$sby_predictor_data
  sby_target_vector <- sby_formula_data$sby_target_vector

  sby_validate_sampling_inputs(sby_predictor_data, sby_target_vector, sby_seed = 1L)
  sby_x_matrix <- sby_adanear_as_numeric_matrix(sby_predictor_data)
  sby_column_names <- sby_adanear_get_column_names(sby_predictor_data)
  colnames(sby_x_matrix) <- sby_column_names
  sby_target_factor <- as.factor(sby_target_vector)
  sby_class_counts <- sby_binary_class_counts_fast(sby_target_factor)
  sby_hpc_apply_env(
    sby_total_threads = sby_total_threads,
    sby_majority_count = sby_class_counts$sby_majority_count,
    sby_minority_count = sby_class_counts$sby_minority_count
  )
  sby_target_name <- sby_formula_data$sby_target_name

  sby_k_neighbor_adanear <- sby_validate_positive_integer_scalar(
    sby_k_neighbor_adanear, "sby_k_neighbor_adanear"
  )
  sby_compute_minority_expansion_count(sby_target_factor, sby_over_ratio)

  if(sby_adanear_hpc_available()){
    sby_balanced_data <- sby_call_native(
      "sby_adasyn_hpc_cpp",
      sby_x_matrix,
      sby_target_factor,
      as.integer(sby_k_neighbor_adanear),
      as.numeric(sby_over_ratio),
      as.integer(sby_total_threads),
      sby_column_names,
      sby_target_name,
      levels(sby_target_factor)
    )
    return(sby_balanced_data)
  }

  sby_balanced_data <- sby_adasyn(
    sby_formula = formula,
    sby_data = .data,
    sby_over_ratio = sby_over_ratio,
    sby_knn_over_k = sby_k_neighbor_adanear,
    sby_audit = FALSE,
    sby_restore_types = FALSE,
    sby_knn_engine = "native",
    sby_knn_distance_metric = "euclidean"
  )
  return(sby_balanced_data)
}
####
## Fim
#
