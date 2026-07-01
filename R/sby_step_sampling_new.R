#' Construtor interno de etapas recipes de sampling
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_subclass Subclasse recipes da etapa
#'
#' @param sby_sampling_method Metodo de sampling da etapa
#'
#' @return Objeto de etapa recipes configurado
#'
#' @noRd
sby_step_sampling_new <- function(
  sby_subclass,
  sby_sampling_method,
  sby_terms,
  sby_role,
  sby_trained,
  sby_columns,
  sby_ratio_over,
  sby_ratio_under,
  sby_knn_over_k,
  sby_knn_under_k,
  sby_seed,
  sby_audit,
  sby_restore_types,
  sby_knn_algorithm,
  sby_knn_engine,
  sby_knn_distance_metric,
  sby_knn_workers,
  sby_knn_parallel_backend,
  sby_knn_hnsw_m,
  sby_knn_hnsw_ef,
  sby_knn_query_chunk_size,
  sby_skip,
  sby_id,
  sby_config_max_threads = NA_integer_
){
  
  # Constroi objeto de etapa recipes com metadados internos de balanceamento
  return(recipes::step(
    subclass                    = sby_subclass,
    terms                       = sby_terms,
    role                        = sby_role,
    trained                     = sby_trained,
    skip                        = sby_skip,
    id                          = sby_id,
    sby_sampling_method         = sby_sampling_method,
    sby_terms                   = sby_terms,
    sby_role                    = sby_role,
    sby_trained                 = sby_trained,
    sby_columns                 = sby_columns,
    sby_ratio_over              = sby_ratio_over,
    sby_ratio_under             = sby_ratio_under,
    sby_knn_over_k              = sby_knn_over_k,
    sby_knn_under_k             = sby_knn_under_k,
    sby_seed                    = sby_seed,
    sby_audit                   = sby_audit,
    sby_restore_types           = sby_restore_types,
    sby_knn_algorithm           = sby_knn_algorithm,
    sby_knn_engine              = sby_knn_engine,
    sby_knn_distance_metric     = sby_knn_distance_metric,
    sby_knn_workers             = sby_knn_workers,
    sby_knn_parallel_backend     = sby_knn_parallel_backend,
    sby_knn_hnsw_m              = sby_knn_hnsw_m,
    sby_knn_hnsw_ef             = sby_knn_hnsw_ef,
    sby_knn_query_chunk_size  = sby_knn_query_chunk_size,
    sby_config_max_threads    = sby_config_max_threads,
    sby_skip                    = sby_skip,
    sby_id                      = sby_id
  ))
}
####
## Fim
#
