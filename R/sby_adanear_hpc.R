#' Atalho HPC do pipeline combinado ADASYN e NearMiss-1
#'
#' @description
#' `sby_adanear_hpc()` e o atalho de alto desempenho do pipeline combinado de
#' balanceamento binario. Ele consolida ADASYN e NearMiss-1 em uma unica passada
#' no espaco padronizado, liquidando os gargalos classicos de normalizacao, copia
#' de memoria e geracao sintetica em infraestrutura NUMA com dois sockets Cascade
#' Lake, 48 nucleos e AVX-512. A funcao representa um atalho: as rotinas originais
#' `sby_adanear()`, `sby_adasyn()` e `sby_nearmiss()` continuam disponiveis.
#'
#' @details
#' Esta funcao substitui internamente a rota denominada "native" como caminho
#' rapido. Quando o motor HPC consolidado esta compilado e carregado, todo o
#' processamento ocorre estritamente sobre a matriz padronizada por z-score, sem
#' dupla normalizacao. As estatisticas iniciais da populacao majoritaria sao
#' calculadas pela Vector Statistics Library, a matriz de distancias usa a
#' identidade D^2 = ||A||^2 + ||B||^2 - 2 A B^T entregue ao cblas_dgemm, a
#' interpolacao lambda do ADASYN usa vdrnguniform e a reversao final do z-score
#' aproveita as unidades FMA do hardware durante a copia para os vetores da lista
#' final. O tibble e montado diretamente no C++ por zero-copy.
#'
#' Quando o motor HPC nao esta disponivel, a funcao recorre de forma transparente
#' a rota classica equivalente com `sby_knn_engine = "native"`, garantindo que
#' tudo continue funcionando.
#'
#' O ambiente de execucao e isolado: as variaveis NUMA e MKL sao capturadas com
#' `Sys.getenv(..., unset = NA)`, injetadas dinamicamente e restauradas por um
#' bloco `on.exit()` inflexivel que remove as variaveis que nao existiam.
#'
#' @param .data Data frame, tibble ou matriz com a coluna de desfecho e os
#'   preditores numericos referenciados em `formula`.
#' @param formula Formula no formato `alvo ~ preditores` que identifica a coluna
#'   de desfecho binario e os preditores numericos.
#' @param sby_k_neighbor_adanear Numero inteiro positivo de vizinhos usado pela
#'   etapa ADASYN. O padrao e `3`.
#' @param sby_k_neighbor_nearmiss Numero inteiro positivo de vizinhos usado pela
#'   etapa NearMiss-1. O padrao e `7`.
#' @param sby_config_max_threads Numero inteiro de threads do motor HPC. O padrao
#'   `-1` deixa a funcao detectar os nucleos fisicos disponiveis.
#' @param sby_over_ratio Fator relativo de expansao da classe minoritaria pela
#'   etapa ADASYN. O padrao e `0.2`.
#' @param sby_under_ratio Razao minima desejada entre minoria e maioria apos o
#'   NearMiss-1. O padrao e `0.5`.
#'
#' @return Tibble balanceado com classe `c("tbl_df", "tbl", "data.frame")`.
#' @export
sby_adanear_hpc <- function(
  .data,
  formula,
  sby_k_neighbor_adanear = 3,
  sby_k_neighbor_nearmiss = 7,
  sby_config_max_threads = -1,
  sby_over_ratio = 0.2,
  sby_under_ratio = 0.5
){
  sby_adanear_check_user_interrupt()

  # Captura o ambiente previo e instala restauro inflexivel imediatamente
  sby_previous_env <- sby_hpc_capture_env()
  on.exit(sby_hpc_restore_env(sby_previous_env), add = TRUE)

  sby_total_threads <- sby_hpc_resolve_threads(sby_config_max_threads)
  sby_hpc_apply_env(sby_total_threads = sby_total_threads, sby_affinity = "scatter")

  # Extrai preditores numericos e alvo binario pelo mesmo contrato da API publica
  sby_formula_data <- sby_extract_formula_data(sby_formula = formula, sby_data = .data)
  sby_predictor_data <- sby_formula_data$sby_predictor_data
  sby_target_vector <- sby_formula_data$sby_target_vector

  sby_validate_sampling_inputs(sby_predictor_data, sby_target_vector, sby_seed = 1L)
  sby_x_matrix <- sby_adanear_as_numeric_matrix(sby_predictor_data)
  sby_column_names <- sby_adanear_get_column_names(sby_predictor_data)
  colnames(sby_x_matrix) <- sby_column_names
  sby_target_factor <- as.factor(sby_target_vector)
  sby_target_name <- sby_formula_data$sby_target_name

  sby_k_neighbor_adanear <- sby_validate_positive_integer_scalar(
    sby_k_neighbor_adanear, "sby_k_neighbor_adanear"
  )
  sby_k_neighbor_nearmiss <- sby_validate_positive_integer_scalar(
    sby_k_neighbor_nearmiss, "sby_k_neighbor_nearmiss"
  )
  sby_compute_minority_expansion_count(sby_target_factor, sby_over_ratio)
  sby_compute_majority_retention_count(sby_target_factor, sby_under_ratio)

  # Caminho rapido: motor HPC consolidado com montagem zero-copy do tibble
  if(sby_adanear_hpc_available()){
    sby_balanced_data <- sby_call_native(
      "sby_adanear_hpc_cpp",
      sby_x_matrix,
      sby_target_factor,
      as.integer(sby_k_neighbor_adanear),
      as.integer(sby_k_neighbor_nearmiss),
      as.numeric(sby_over_ratio),
      as.numeric(sby_under_ratio),
      as.integer(sby_total_threads),
      sby_column_names,
      sby_target_name,
      levels(sby_target_factor)
    )
    return(sby_balanced_data)
  }

  # Fallback transparente: rota classica equivalente com engine nativa exata
  sby_balanced_data <- sby_adanear(
    sby_formula = formula,
    sby_data = .data,
    sby_over_ratio = sby_over_ratio,
    sby_under_ratio = sby_under_ratio,
    sby_knn_over_k = sby_k_neighbor_adanear,
    sby_knn_under_k = sby_k_neighbor_nearmiss,
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
