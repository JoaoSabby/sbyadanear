#' Obter subclasse recipes primaria da etapa
#'
#' @param sby_x Objeto de etapa recipes
#'
#' @return String com subclasse primaria
#'
#' @noRd
sby_step_sampling_subclass <- function(sby_x){
  return(sub(
    pattern = "^step_",
    replacement = "",
    x = class(sby_x)[[1L]]
  ))
}

#' Preparar etapa recipes de sampling
#'
#' @param x Objeto de etapa nao treinado
#'
#' @param training Dados de treinamento
#'
#' @param info Metadados recipes
#'
#' @param sby_step_name Nome publico da etapa para mensagens
#'
#' @return Objeto treinado
#'
#' @noRd
sby_prep_step_sampling <- function(x, training, info, sby_step_name){
  
  # Normaliza argumentos S3 para nomes internos do pacote
  sby_x        <- x
  sby_training <- training
  sby_info     <- info

  # Verifica se ha solicitacao de interrupcao antes da selecao de colunas
  sby_adanear_check_user_interrupt()

  # Avalia seletores recipes para identificar coluna de desfecho
  sby_selected_columns <- recipes::recipes_eval_select(
    quos = sby_x$sby_terms,
    data = sby_training,
    info = sby_info
  )

  # Verifica se exatamente uma coluna de desfecho foi selecionada
  if(length(sby_selected_columns) != 1L){

    # Aborta quando a selecao do desfecho nao e univoca
    sby_adanear_abort(
      sby_message = paste0("'", sby_step_name, "' deve selecionar exatamente uma coluna de desfecho")
    )
  }

  # Retorna nova etapa marcada como treinada com a coluna selecionada
  return(sby_step_sampling_new(
    sby_subclass                = sby_step_sampling_subclass(sby_x),
    sby_sampling_method         = sby_x$sby_sampling_method,
    sby_terms                   = sby_x$sby_terms,
    sby_role                    = sby_x$sby_role,
    sby_trained                 = TRUE,
    sby_columns                 = names(sby_selected_columns),
    sby_over_ratio              = sby_x$sby_over_ratio,
    sby_under_ratio             = sby_x$sby_under_ratio,
    sby_knn_over_k              = sby_x$sby_knn_over_k,
    sby_knn_under_k             = sby_x$sby_knn_under_k,
    sby_seed                    = sby_x$sby_seed,
    sby_audit                   = sby_x$sby_audit,
    sby_restore_types           = sby_x$sby_restore_types,
    sby_knn_algorithm           = sby_x$sby_knn_algorithm,
    sby_knn_engine              = sby_x$sby_knn_engine,
    sby_knn_distance_metric     = sby_x$sby_knn_distance_metric,
    sby_knn_workers             = sby_x$sby_knn_workers,
    sby_knn_parallel_backend     = sby_x$sby_knn_parallel_backend,
    sby_knn_hnsw_m              = sby_x$sby_knn_hnsw_m,
    sby_knn_hnsw_ef             = sby_x$sby_knn_hnsw_ef,
    sby_knn_query_chunk_size  = sby_x$sby_knn_query_chunk_size,
    sby_skip                    = sby_x$sby_skip,
    sby_id                      = sby_x$sby_id
  ))
}

#' Aplicar etapa recipes de sampling
#'
#' @param object Objeto de etapa treinado
#'
#' @param new_data Dados novos
#'
#' @param sby_step_name Nome publico da etapa para mensagens
#'
#' @return Tibble balanceado ou lista de auditoria
#'
#' @noRd
sby_bake_step_sampling <- function(object, new_data, sby_step_name){
  
  # Normaliza argumentos S3 para nomes internos do pacote
  sby_object   <- object
  sby_new_data <- new_data

  # Verifica se ha solicitacao de interrupcao antes do balanceamento
  sby_adanear_check_user_interrupt()

  # Verifica se a etapa foi treinada antes de aplicar bake
  if(!isTRUE(sby_object$sby_trained)){

    # Aborta quando a etapa nao passou por prep
    sby_adanear_abort(
      sby_message = paste0("'", sby_step_name, "' precisa ser treinado com prep() antes de bake()")
    )
  }

  # Recupera coluna de desfecho selecionada no treinamento
  sby_target_column <- sby_object$sby_columns[[1L]]

  # Verifica se a coluna de desfecho existe nos novos dados
  if(!sby_target_column %in% names(sby_new_data)){

    # Aborta quando o desfecho selecionado esta ausente em new_data
    sby_adanear_abort(
      sby_message = paste0(
        "Coluna de desfecho nao encontrada em 'new_data': ",
        sby_target_column
      )
    )
  }

  # Define nomes de preditores excluindo a coluna de desfecho
  sby_original_names   <- names(sby_new_data)
  sby_predictor_names  <- setdiff(
    x = sby_original_names,
    y = sby_target_column
  )

  # Verifica se existe ao menos uma coluna preditora
  if(length(sby_predictor_names) < 1L){

    # Aborta quando apenas o desfecho esta disponivel
    sby_adanear_abort(
      sby_message = paste0("'", sby_step_name, "' requer ao menos uma coluna preditora")
    )
  }

  # Monta formula e dados conforme contrato publico das funcoes de sampling
  sby_formula <- stats::reformulate(
    termlabels = ".",
    response = sby_target_column
  )
  sby_data <- as.data.frame(
    x = sby_new_data[, c(sby_predictor_names, sby_target_column), drop = FALSE]
  )

  # Executa metodo de sampling configurado na etapa
  if(identical(sby_object$sby_sampling_method, "adasyn")){

    # Executa sobreamostragem ADASYN
    sby_sampling_result <- sby_adasyn(
      sby_formula             = sby_formula,
      sby_data                = sby_data,
      sby_over_ratio          = sby_object$sby_over_ratio,
      sby_knn_over_k          = sby_object$sby_knn_over_k,
      sby_seed                = sby_object$sby_seed,
      sby_audit               = TRUE,
      sby_restore_types       = sby_object$sby_restore_types,
      sby_knn_algorithm       = sby_object$sby_knn_algorithm,
      sby_knn_engine          = sby_object$sby_knn_engine,
      sby_knn_distance_metric = sby_object$sby_knn_distance_metric,
      sby_knn_workers         = sby_object$sby_knn_workers,
      sby_knn_parallel_backend = sby_object$sby_knn_parallel_backend,
      sby_knn_hnsw_m          = sby_object$sby_knn_hnsw_m,
      sby_knn_hnsw_ef         = sby_object$sby_knn_hnsw_ef,
      sby_knn_query_chunk_size  = sby_object$sby_knn_query_chunk_size
    )
  }else if(identical(sby_object$sby_sampling_method, "nearmiss")){

    # Executa subamostragem NearMiss-1
    sby_sampling_result <- sby_nearmiss(
      sby_formula             = sby_formula,
      sby_data                = sby_data,
      sby_under_ratio         = sby_object$sby_under_ratio,
      sby_knn_under_k         = sby_object$sby_knn_under_k,
      sby_seed                = sby_object$sby_seed,
      sby_audit               = TRUE,
      sby_restore_types       = sby_object$sby_restore_types,
      sby_knn_algorithm       = sby_object$sby_knn_algorithm,
      sby_knn_engine          = sby_object$sby_knn_engine,
      sby_knn_distance_metric = sby_object$sby_knn_distance_metric,
      sby_knn_workers         = sby_object$sby_knn_workers,
      sby_knn_parallel_backend = sby_object$sby_knn_parallel_backend,
      sby_knn_hnsw_m          = sby_object$sby_knn_hnsw_m,
      sby_knn_hnsw_ef         = sby_object$sby_knn_hnsw_ef,
      sby_knn_query_chunk_size  = sby_object$sby_knn_query_chunk_size
    )
  }else{

    # Executa pipeline hibrido ADASYN + NearMiss-1
    sby_sampling_result <- sby_adanear(
      sby_formula             = sby_formula,
      sby_data                = sby_data,
      sby_over_ratio          = sby_object$sby_over_ratio,
      sby_under_ratio         = sby_object$sby_under_ratio,
      sby_knn_over_k          = sby_object$sby_knn_over_k,
      sby_knn_under_k         = sby_object$sby_knn_under_k,
      sby_seed                = sby_object$sby_seed,
      sby_audit               = TRUE,
      sby_restore_types       = sby_object$sby_restore_types,
      sby_knn_algorithm       = sby_object$sby_knn_algorithm,
      sby_knn_engine          = sby_object$sby_knn_engine,
      sby_knn_distance_metric = sby_object$sby_knn_distance_metric,
      sby_knn_workers         = sby_object$sby_knn_workers,
      sby_knn_parallel_backend = sby_object$sby_knn_parallel_backend,
      sby_knn_hnsw_m          = sby_object$sby_knn_hnsw_m,
      sby_knn_hnsw_ef         = sby_object$sby_knn_hnsw_ef,
      sby_knn_query_chunk_size  = sby_object$sby_knn_query_chunk_size
    )
  }

  # Verifica se ha solicitacao de interrupcao apos o balanceamento
  sby_adanear_check_user_interrupt()

  # Retorna auditoria completa quando configurado na etapa
  if(isTRUE(sby_object$sby_audit)){

    # Entrega resultado auditavel ao chamador
    return(sby_sampling_result)
  }

  # Retorna apenas dados balanceados no fluxo recipes padrao
  return(sby_sampling_result$sby_balanced_data)
}

#' Organizar metadados de etapa recipes de sampling
#'
#' @param x Objeto de etapa recipes
#'
#' @return Data frame de metadados
#'
#' @noRd
sby_tidy_step_sampling <- function(x){
  
  # Normaliza argumento S3 para nome interno do pacote
  sby_x <- x

  # Define termos exibidos conforme estado de treinamento
  if(isTRUE(sby_x$sby_trained)){
    sby_terms <- sby_x$sby_columns
  }else{
    sby_terms <- recipes::sel2char(
      x = sby_x$sby_terms
    )
  }

  # Retorna metadados tabulares da etapa
  return(data.frame(
    sby_terms               = sby_terms,
    sby_sampling_method     = sby_x$sby_sampling_method,
    sby_over_ratio          = sby_x$sby_over_ratio,
    sby_under_ratio         = sby_x$sby_under_ratio,
    sby_knn_over_k          = sby_x$sby_knn_over_k,
    sby_knn_under_k         = sby_x$sby_knn_under_k,
    sby_knn_query_chunk_size  = sby_x$sby_knn_query_chunk_size,
    sby_audit               = sby_x$sby_audit,
    sby_id                  = sby_x$sby_id,
    stringsAsFactors        = FALSE
  ))
}

#' Imprimir etapa recipes de sampling
#'
#' @param x Objeto de etapa recipes
#'
#' @param width Largura de impressao
#'
#' @param sby_title Titulo da etapa
#'
#' @return Objeto invisivel
#'
#' @noRd
sby_print_step_sampling <- function(x, width, sby_title){
  
  # Normaliza argumentos S3 para nomes internos do pacote
  sby_x     <- x
  sby_width <- width

  # Define colunas exibidas conforme estado de treinamento
  if(isTRUE(sby_x$sby_trained)){
    sby_columns <- sby_x$sby_columns
  }else{
    sby_columns <- recipes::sel2char(
      x = sby_x$sby_terms
    )
  }

  # Imprime resumo padronizado da etapa recipes
  recipes::print_step(
    tr_obj = sby_columns,
    untr_obj = sby_x$sby_terms,
    trained = sby_x$sby_trained,
    title = sby_title,
    width = sby_width
  )

  # Retorna objeto original invisivelmente apos impressao
  return(invisible(sby_x))
}
####
## Fim
#
