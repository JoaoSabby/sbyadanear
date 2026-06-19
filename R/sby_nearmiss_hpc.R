#' Atalho HPC do undersampling NearMiss-1
#'
#' @description
#' `sby_nearmiss_hpc()` e o atalho de alto desempenho do undersampling NearMiss-1.
#' Ranqueia e retem as linhas majoritarias estritamente no espaco padronizado via
#' MKL VSL e cblas_sgemm. O C++ retorna apenas os indices das linhas retidas;
#' a reconstrucao do tibble ocorre na camada R diretamente a partir dos dados
#' originais, sem aritmética alguma.
#'
#' @details
#' Controla temporariamente apenas `MKL_NUM_THREADS`, `OMP_NUM_THREADS` e
#' `MKL_NUM_STRIPES`, restaurando os valores originais por `on.exit()` inflexivel.
#' O ambiente e configurado antes de qualquer operacao matricial.
#'
#' @param .data Data frame ou tibble com a coluna de desfecho e preditores
#'   numericos referenciados em `formula`.
#' @param formula Formula no formato `alvo ~ preditores`.
#' @param sby_k_neighbor_nearmiss Numero inteiro positivo de vizinhos do
#'   NearMiss-1. Padrao: `7`.
#' @param sby_under_ratio Razao minima minoria/maioria apos o NearMiss-1.
#'   Padrao: `0.5`.
#' @param sby_config_max_threads Numero inteiro de threads do motor HPC. `-1`
#'   detecta os nucleos fisicos disponíveis. Padrao: `-1`.
#' @param sby_seed Semente inteira para reproducibilidade. Padrao:
#'   `sample.int(10L^5L, 1L)`.
#'
#' @return Tibble balanceado com classe `c("tbl_df", "tbl", "data.frame")`.
#' @export
sby_nearmiss_hpc <- function(
  .data,
  formula,
  sby_k_neighbor_nearmiss = 7,
  sby_under_ratio         = 0.5,
  sby_config_max_threads  = -1,
  sby_seed                = sample.int(10L^5L, 1L)
){
  sby_adanear_check_user_interrupt()

  sby_original_column_order <- colnames(.data)

  # Nao altera variaveis de ambiente MKL/OMP dentro da chamada;
  # respeita a configuracao externa do runtime HPC.
  sby_total_threads <- sby_hpc_resolve_threads(sby_config_max_threads)

  sby_formula_data            <- sby_extract_formula_data(sby_formula = formula, sby_data = .data)
  sby_original_predictor_data <- sby_formula_data$sby_predictor_data
  sby_target_vector           <- sby_formula_data$sby_target_vector
  sby_target_name             <- sby_formula_data$sby_target_name

  # Captura os levels originais ANTES de qualquer as.factor() para preservar
  # a classe, a ordem e os labels exatos do factor de entrada.
  # c(factor, character) destruiria o factor retornando codigos numericos.
  sby_original_levels <- if (is.factor(sby_target_vector)) {
    levels(sby_target_vector)
  } else {
    unique(as.character(sby_target_vector))
  }

  sby_seed <- sby_validate_seed(sby_seed = sby_seed)
  sby_validate_sampling_inputs(sby_original_predictor_data, sby_target_vector, sby_seed = sby_seed)

  sby_x_matrix     <- sby_adanear_as_numeric_matrix(sby_original_predictor_data)
  sby_column_names <- sby_adanear_get_column_names(sby_original_predictor_data)

  # Usa os levels originais para nao reordenar alfabeticamente
  sby_target_factor <- factor(sby_target_vector, levels = sby_original_levels)
  sby_class_counts  <- sby_binary_class_counts_fast(sby_target_factor)

  # O runtime MKL/OpenMP deve ser configurado externamente pelo usuario HPC.

  # Indices 1-based das linhas da minoria no conjunto original
  sby_minority_level_int <- as.integer(sby_class_counts$sby_minority_level)
  sby_minority_idx       <- which(as.integer(sby_target_factor) == sby_minority_level_int)

  sby_k_neighbor_nearmiss <- sby_validate_positive_integer_scalar(
    sby_k_neighbor_nearmiss, "sby_k_neighbor_nearmiss"
  )

  if (!sby_adanear_hpc_available()) {
    sby_adanear_abort(
      "Motor HPC nao disponivel. Compile o pacote com suporte a MKL/AVX-512.",
      call = sys.call()
    )
  }

  set.seed(sby_seed)
  sby_hpc_result <- sby_call_native(
    "sby_nearmiss_hpc_cpp",
    sby_x_matrix,
    sby_target_factor,
    as.integer(sby_k_neighbor_nearmiss),
    as.numeric(sby_under_ratio),
    as.integer(sby_total_threads),
    sby_column_names,
    levels(sby_target_factor)
  )
  # Retorno esperado de sby_nearmiss_hpc_cpp:
  #   $sby_retained_majority_idx — IntegerVector (indices 1-based no original)
  #   $sby_scaling_info          — List(centers, scales)

  # --- Reconstrucao das 2 partes na camada R (zero aritmetica) ---

  # Parte 1: maioria remanescente — indice direto no original
  sby_maj_rows   <- sby_original_predictor_data[
    sby_hpc_result$sby_retained_majority_idx, , drop = FALSE
  ]
  sby_maj_target <- sby_target_vector[sby_hpc_result$sby_retained_majority_idx]

  # Parte 2: minoria original — integra
  sby_min_rows   <- sby_original_predictor_data[sby_minority_idx, , drop = FALSE]
  sby_min_target <- sby_target_vector[sby_minority_idx]

  sby_final_predictors <- rbind(sby_maj_rows, sby_min_rows)
  rownames(sby_final_predictors) <- NULL

  # Reconstroi vetor alvo como factor com os levels originais para preservar
  # a classe, a ordem e os labels exatos — evitando que c(factor, character)
  # retorne codigos numericos como character.
  # NearMiss nao gera sinteticas: apenas concatena os vetores originais.
  sby_final_target <- factor(
    c(
      as.character(sby_maj_target),
      as.character(sby_min_target)
    ),
    levels = sby_original_levels
  )

  sby_balanced_data <- sby_build_balanced_tibble(
    sby_predictor_data = sby_final_predictors,
    sby_target_vector  = sby_final_target
  )

  if (!identical(sby_target_name, "TARGET")) {
    names(sby_balanced_data)[names(sby_balanced_data) == "TARGET"] <- sby_target_name
  }

  sby_balanced_data <- collapse::fselect(
    .x = sby_balanced_data,
    sby_original_column_order
  )

  return(sby_balanced_data)
}
####
## Fim
#
