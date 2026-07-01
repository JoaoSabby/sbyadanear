#' Adicionar etapas recipes HPC de balanceamento
#'
#' @description
#' `sby_step_adasyn_hpc()`, `sby_step_nearmiss_hpc()` e
#' `sby_step_adanear_hpc()` adicionam etapas supervisionadas compatíveis com
#' `recipes` que executam, em `bake()`, as rotinas HPC correspondentes:
#' `sby_adasyn_hpc()`, `sby_nearmiss_hpc()` e `sby_adanear_hpc()`.
#'
#' @param recipe Objeto `recipe` que recebera a etapa.
#' @param ... Seletores `recipes` usados para identificar uma coluna de desfecho.
#' @param role Papel armazenado na etapa.
#' @param trained Indicador interno de treinamento da etapa.
#' @param columns Coluna de desfecho resolvida durante `prep()`.
#' @param sby_adasyn_ratio Razao de sobreamostragem ADASYN.
#' @param sby_nearmiss_ratio Razao de retencao NearMiss-1.
#' @param sby_k_adasyn Numero inteiro positivo de vizinhos para ADASYN.
#' @param sby_k_nearmiss Numero inteiro positivo de vizinhos para NearMiss-1.
#' @param sby_config_max_threads Numero inteiro de threads do motor HPC. `-1` detecta os nucleos fisicos disponiveis.
#' @param sby_seed Semente inteira.
#' @param sby_audit Indicador logico que retorna uma lista com `sby_balanced_data` quando `TRUE`.
#' @param sby_restore_types Mantido para compatibilidade com metadados da etapa.
#' @param skip Indicador logico que define se a etapa deve ser ignorada em novos dados.
#' @param id Identificador recipes da etapa.
#'
#' @return Objeto `recipe` com a etapa HPC adicionada ao pipeline.
#'
#' @export
sby_step_adasyn_hpc <- function(
  recipe,
  ...,
  role = NA,
  trained = FALSE,
  columns = NULL,
  sby_adasyn_ratio = 0.2,
  sby_k_adasyn = 3L,
  sby_config_max_threads = -1L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_restore_types = TRUE,
  skip = TRUE,
  id = recipes::rand_id("adasyn_hpc")
){
  sby_step_hpc_add(
    recipe = recipe, ..., role = role, trained = trained, columns = columns,
    sby_sampling_method = "adasyn_hpc", sby_adasyn_ratio = sby_adasyn_ratio,
    sby_nearmiss_ratio = NA_real_, sby_k_adasyn = sby_k_adasyn,
    sby_k_nearmiss = NA_integer_, sby_config_max_threads = sby_config_max_threads,
    sby_seed = sby_seed, sby_audit = sby_audit, sby_restore_types = sby_restore_types,
    skip = skip, id = id, sby_required_pkgs = required_pkgs.step_sby_step_adasyn_hpc()
  )
}

#' @rdname sby_step_adasyn_hpc
#' @export
sby_step_nearmiss_hpc <- function(
  recipe,
  ...,
  role = NA,
  trained = FALSE,
  columns = NULL,
  sby_nearmiss_ratio = 1,
  sby_k_nearmiss = 7L,
  sby_config_max_threads = -1L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_restore_types = TRUE,
  skip = TRUE,
  id = recipes::rand_id("nearmiss_hpc")
){
  sby_step_hpc_add(
    recipe = recipe, ..., role = role, trained = trained, columns = columns,
    sby_sampling_method = "nearmiss_hpc", sby_adasyn_ratio = NA_real_,
    sby_nearmiss_ratio = sby_nearmiss_ratio, sby_k_adasyn = NA_integer_,
    sby_k_nearmiss = sby_k_nearmiss, sby_config_max_threads = sby_config_max_threads,
    sby_seed = sby_seed, sby_audit = sby_audit, sby_restore_types = sby_restore_types,
    skip = skip, id = id, sby_required_pkgs = required_pkgs.step_sby_step_nearmiss_hpc()
  )
}

#' @rdname sby_step_adasyn_hpc
#' @export
sby_step_adanear_hpc <- function(
  recipe,
  ...,
  role = NA,
  trained = FALSE,
  columns = NULL,
  sby_adasyn_ratio = 0.2,
  sby_nearmiss_ratio = 1,
  sby_k_adasyn = 3L,
  sby_k_nearmiss = 7L,
  sby_config_max_threads = -1L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_restore_types = TRUE,
  skip = TRUE,
  id = recipes::rand_id("adanear_hpc")
){
  sby_step_hpc_add(
    recipe = recipe, ..., role = role, trained = trained, columns = columns,
    sby_sampling_method = "adanear_hpc", sby_adasyn_ratio = sby_adasyn_ratio,
    sby_nearmiss_ratio = sby_nearmiss_ratio, sby_k_adasyn = sby_k_adasyn,
    sby_k_nearmiss = sby_k_nearmiss, sby_config_max_threads = sby_config_max_threads,
    sby_seed = sby_seed, sby_audit = sby_audit, sby_restore_types = sby_restore_types,
    skip = skip, id = id, sby_required_pkgs = required_pkgs.step_sby_step_adanear_hpc()
  )
}

sby_step_hpc_add <- function(recipe, ..., role, trained, columns, sby_sampling_method,
                             sby_adasyn_ratio, sby_nearmiss_ratio, sby_k_adasyn,
                             sby_k_nearmiss, sby_config_max_threads, sby_seed,
                             sby_audit, sby_restore_types, skip, id, sby_required_pkgs){
  sby_adanear_check_user_interrupt()
  recipes::recipes_pkg_check(sby_required_pkgs)
  sby_terms <- rlang::enquos(...)
  sby_audit <- sby_validate_logical_scalar(sby_value = sby_audit, sby_name = "sby_audit")
  sby_restore_types <- sby_validate_logical_scalar(sby_value = sby_restore_types, sby_name = "sby_restore_types")
  skip <- sby_validate_logical_scalar(sby_value = skip, sby_name = "skip")
  sby_seed <- sby_validate_seed(sby_seed = sby_seed)
  if(!is.na(sby_k_adasyn)) sby_k_adasyn <- sby_validate_positive_integer_scalar(sby_k_adasyn, "sby_k_adasyn")
  if(!is.na(sby_k_nearmiss)) sby_k_nearmiss <- sby_validate_positive_integer_scalar(sby_k_nearmiss, "sby_k_nearmiss")
  recipes::add_step(
    rec = recipe,
    object = sby_step_sampling_new(
      sby_subclass = paste0("sby_step_", sby_sampling_method),
      sby_sampling_method = sby_sampling_method,
      sby_terms = sby_terms, sby_role = role, sby_trained = trained,
      sby_columns = columns, sby_adasyn_ratio = sby_adasyn_ratio,
      sby_nearmiss_ratio = sby_nearmiss_ratio, sby_knn_over_k = sby_k_adasyn,
      sby_knn_under_k = sby_k_nearmiss, sby_seed = sby_seed,
      sby_audit = sby_audit, sby_restore_types = sby_restore_types,
      sby_knn_algorithm = NA_character_, sby_knn_engine = NA_character_,
      sby_knn_distance_metric = NA_character_, sby_knn_workers = NA_integer_,
      sby_knn_parallel_backend = NA_character_, sby_knn_hnsw_m = NA_integer_,
      sby_knn_hnsw_ef = NA_integer_, sby_knn_query_chunk_size = NA_integer_,
      sby_skip = skip, sby_id = id, sby_config_max_threads = sby_config_max_threads
    )
  )
}
