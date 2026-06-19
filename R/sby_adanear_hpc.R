#' Atalho HPC do pipeline combinado ADASYN e NearMiss-1
#'
#' @description
#' `sby_adanear_hpc()` e o atalho de alto desempenho do pipeline combinado de
#' balanceamento binario. Consolida ADASYN e NearMiss-1 em uma unica passada no
#' espaco padronizado, usando MKL VSL para estatisticas, cblas_sgemm para
#' distancias e interpolacao lambda com vdrnguniform. A despadronizacao das
#' sinteticas ocorre inteiramente no C++ via FMA AVX-512. A reconstrucao final
#' do tibble acontece na camada R, preservando os tipos originais das colunas.
#'
#' @details
#' Controla temporariamente apenas `MKL_NUM_THREADS`, `OMP_NUM_THREADS` e
#' `MKL_NUM_STRIPES`, restaurando os valores originais por `on.exit()` inflexivel.
#' O ambiente e configurado antes de qualquer operacao matricial para garantir
#' que o MKL leia os valores corretos desde o primeiro kernel.
#'
#' @param .data Data frame ou tibble com a coluna de desfecho e preditores
#'   numericos referenciados em `formula`.
#' @param formula Formula no formato `alvo ~ preditores`.
#' @param sby_k_neighbor_adanear Numero inteiro positivo de vizinhos da etapa
#'   ADASYN. Padrao: `3`.
#' @param sby_k_neighbor_nearmiss Numero inteiro positivo de vizinhos da etapa
#'   NearMiss-1. Padrao: `7`.
#' @param sby_config_max_threads Numero inteiro de threads do motor HPC. `-1`
#'   detecta os nucleos fisicos disponíveis. Padrao: `-1`.
#' @param sby_seed Semente inteira para o gerador de numeros pseudo-aleatorios
#'   do ADASYN. Padrao: `sample.int(10L^5L, 1L)`.
#' @param sby_over_ratio Fator de expansao da classe minoritaria pelo ADASYN.
#'   Deve ser estritamente positivo. Padrao: `0.2`.
#' @param sby_under_ratio Razao minima minoria/maioria apos o NearMiss-1.
#'   Padrao: `0.5`.
#'
#' @return Tibble balanceado com classe `c("tbl_df", "tbl", "data.frame")`.
#' @export
sby_adanear_hpc <- function(
  .data,
  formula,
  sby_k_neighbor_adanear  = 3,
  sby_k_neighbor_nearmiss = 7,
  sby_config_max_threads  = -1,
  sby_seed                = sample.int(10L^5L, 1L),
  sby_over_ratio          = 0.2,
  sby_under_ratio         = 0.5
){
  sby_adanear_check_user_interrupt()

  # Ordem das colunas do input para recompor o output na mesma sequencia
  sby_original_column_order <- colnames(.data)

  # Nao altera variaveis de ambiente MKL/OMP dentro da chamada;
  # respeita a configuracao externa do runtime HPC.
  sby_total_threads <- sby_hpc_resolve_threads(sby_config_max_threads)

  # --- Validacoes antes de qualquer operacao matricial ---
  if (!is.numeric(sby_over_ratio) || length(sby_over_ratio) != 1L ||
      is.na(sby_over_ratio) || sby_over_ratio <= 0) {
    sby_adanear_abort(
      "sby_over_ratio deve ser um numero positivo maior que zero.",
      call = sys.call()
    )
  }

  # Extrai preditores e alvo
  sby_formula_data         <- sby_extract_formula_data(sby_formula = formula, sby_data = .data)
  sby_original_predictor_data <- sby_formula_data$sby_predictor_data
  sby_target_vector        <- sby_formula_data$sby_target_vector
  sby_target_name          <- sby_formula_data$sby_target_name

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
  sby_target_factor  <- factor(sby_target_vector, levels = sby_original_levels)
  sby_class_counts   <- sby_binary_class_counts_fast(sby_target_factor)

  # O runtime MKL/OpenMP deve ser configurado externamente pelo usuario HPC.

  sby_type_info <- sby_infer_numeric_column_types(sby_original_predictor_data)

  # Indices 1-based das linhas da minoria no conjunto original
  sby_minority_level_int <- as.integer(sby_class_counts$sby_minority_level)
  sby_minority_idx       <- which(as.integer(sby_target_factor) == sby_minority_level_int)

  sby_k_neighbor_adanear  <- sby_validate_positive_integer_scalar(
    sby_k_neighbor_adanear, "sby_k_neighbor_adanear"
  )
  sby_k_neighbor_nearmiss <- sby_validate_positive_integer_scalar(
    sby_k_neighbor_nearmiss, "sby_k_neighbor_nearmiss"
  )

  # Motor HPC obrigatorio: sem fallback de engine generico
  if (!sby_adanear_hpc_available()) {
    sby_adanear_abort(
      "Motor HPC nao disponivel. Compile o pacote com suporte a MKL/AVX-512.",
      call = sys.call()
    )
  }

  set.seed(sby_seed)
  sby_hpc_result <- sby_call_native(
    "sby_adanear_hpc_result_cpp",
    sby_x_matrix,
    sby_target_factor,
    as.integer(sby_k_neighbor_adanear),
    as.integer(sby_k_neighbor_nearmiss),
    as.numeric(sby_over_ratio),
    as.numeric(sby_under_ratio),
    as.integer(sby_total_threads),
    sby_column_names,
    levels(sby_target_factor)
  )
  # Retorno esperado de sby_adanear_hpc_result_cpp:
  #   $sby_synthetic_rows        — NumericMatrix (double, despadronizado no C++)
  #   $sby_retained_majority_idx — IntegerVector (indices 1-based no original)
  #   $sby_target_synthetic      — IntegerVector (codigos de nivel das sinteticas)
  #   $sby_scaling_info          — List(centers, scales)

  # --- Reconstrucao das 3 partes na camada R ---

  # Parte 1: maioria remanescente — indice direto no original, zero aritmetica
  sby_maj_rows   <- sby_original_predictor_data[
    sby_hpc_result$sby_retained_majority_idx, , drop = FALSE
  ]
  sby_maj_target <- sby_target_vector[sby_hpc_result$sby_retained_majority_idx]

  # Parte 2: minoria original — integra, sem qualquer transformacao
  sby_min_rows   <- sby_original_predictor_data[sby_minority_idx, , drop = FALSE]
  sby_min_target <- sby_target_vector[sby_minority_idx]

  # Parte 3: sinteticas — chegam como double apos despadronizacao no C++;
  #           apenas restauro de tipos (integer -> integer, etc.)
  if (nrow(sby_hpc_result$sby_synthetic_rows) > 0L) {
    sby_syn_df <- sby_restore_numeric_column_types(
      as.data.frame(sby_hpc_result$sby_synthetic_rows, stringsAsFactors = FALSE),
      sby_type_info,
      TRUE
    )
  } else {
    sby_syn_df <- sby_original_predictor_data[0L, , drop = FALSE]
  }

  # O numero de linhas sinteticas deve ser pelo menos 1 (ceiling garantido no C++)
  # mas se vier zero (over_ratio muito pequeno), o rbind ainda e valido.
  sby_final_predictors <- rbind(sby_maj_rows, sby_min_rows, sby_syn_df)
  rownames(sby_final_predictors) <- NULL

  # Reconstroi vetor alvo: labels das sinteticas via levels originais.
  # O factor e reconstituido com os levels originais para preservar a classe,
  # a ordem e os labels exatos — evitando que c(factor, character) retorne
  # codigos numericos como character.
  sby_syn_target_labels <- sby_original_levels[
    sby_hpc_result$sby_target_synthetic
  ]
  sby_final_target <- factor(
    c(
      as.character(sby_maj_target),
      as.character(sby_min_target),
      sby_syn_target_labels
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
