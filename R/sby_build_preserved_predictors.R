#' Montar preditores preservando linhas originais observadas
#'
#' @details
#' Constroi o data frame final do wrapper tabular reaproveitando linhas
#' originais (sem reverter z-score) e inserindo apenas as linhas sinteticas
#' nas posicoes corretas. A versao vetorizada substitui um loop linha a
#' linha com `[[.data.frame` + `do.call(rbind, ...)`, que dominava o tempo
#' de execucao da API tabular (medido em ~24x mais lento que a API de
#' matriz para n = 10.000).
#'
#' @param sby_original_predictor_data Preditores originais antes do z-score
#' @param sby_retained_index Indices finais no espaco pos-ADASYN
#' @param sby_final_scaled Matriz final escalada na mesma ordem dos indices retidos
#' @param sby_scaling_info Parametros de z-score
#' @param sby_type_info Metadados de tipos numericos
#' @param sby_restore_types Indicador de restauracao das linhas sinteticas
#'
#' @return Data frame de preditores finais
#' @noRd
sby_build_preserved_predictors <- function(
  sby_original_predictor_data,
  sby_retained_index,
  sby_final_scaled,
  sby_scaling_info,
  sby_type_info,
  sby_restore_types
){
  sby_original_predictor_data <- as.data.frame(
    x = sby_original_predictor_data,
    stringsAsFactors = FALSE
  )
  sby_retained_index <- as.integer(sby_retained_index)
  sby_original_row_count <- collapse::fnrow(sby_original_predictor_data)
  sby_is_original <- sby_retained_index <= sby_original_row_count

  # Caminho rapido 1: nenhuma linha sintetica retida (NearMiss puro). Basta
  # indexar o data frame original em ordem vetorizada e zerar rownames.
  if(!any(!sby_is_original)){
    sby_final_predictors <- sby_original_predictor_data[sby_retained_index, , drop = FALSE]
    rownames(sby_final_predictors) <- NULL
    return(sby_final_predictors)
  }

  # Caminho rapido 2: todos os indices iniciais sao originais e sinteticos
  # vem em bloco no final (caso comum apos ADASYN puro). Resolve sinteticos
  # em batch e concatena com as originais por bloco.
  sby_synth_positions <- which(!sby_is_original)
  sby_orig_positions <- which(sby_is_original)
  sby_synthetic_scaled <- sby_final_scaled[sby_synth_positions, , drop = FALSE]
  sby_synthetic_original <- sby_revert_z_score_scaling_matrix(
    sby_x_matrix = sby_synthetic_scaled,
    sby_scaling_info = sby_scaling_info
  )
  sby_synthetic_predictors <- if(isTRUE(sby_restore_types)){
    sby_restore_numeric_column_types(sby_synthetic_original, sby_type_info, TRUE)
  }else{
    sby_out <- as.data.frame(sby_synthetic_original, stringsAsFactors = FALSE)
    names(sby_out) <- sby_type_info$sby_column_name
    sby_out
  }

  # Caso ADASYN puro: ordem e [originais retidas em bloco, sinteticas em bloco].
  # Usamos rbind apenas duas vezes (originais inteiras + sinteticas inteiras),
  # nao N vezes uma a uma.
  if(length(sby_orig_positions) > 0L &&
     identical(sby_orig_positions, seq_len(length(sby_orig_positions))) &&
     identical(sby_synth_positions, seq.int(length(sby_orig_positions) + 1L, length(sby_retained_index)))){
    sby_orig_block <- sby_original_predictor_data[sby_retained_index[sby_orig_positions], , drop = FALSE]
    sby_final_predictors <- data.table::rbindlist(
      l = list(sby_orig_block, sby_synthetic_predictors),
      use.names = TRUE,
      fill = FALSE
    )
    sby_final_predictors <- as.data.frame(sby_final_predictors, stringsAsFactors = FALSE)
    rownames(sby_final_predictors) <- NULL
    return(sby_final_predictors)
  }

  # Caso ADANEAR: ordem geral pode intercalar originais e sinteticas. Pre-
  # aloca o data frame final com a forma das colunas das originais (mesma
  # estrutura/typos) e preenche em duas atribuicoes vetorizadas, sem loop.
  sby_total <- length(sby_retained_index)
  # Esqueleto: aloca data frame com o tipo das colunas originais ja correto.
  sby_final_predictors <- sby_original_predictor_data[rep.int(1L, sby_total), , drop = FALSE]

  # Mapeia originais retidas: copia direta por indexacao vetorizada.
  sby_final_predictors[sby_orig_positions, ] <-
    sby_original_predictor_data[sby_retained_index[sby_orig_positions], , drop = FALSE]

  # Sinteticas: alinhamento por posicao com sby_synth_positions.
  sby_final_predictors[sby_synth_positions, ] <- sby_synthetic_predictors

  rownames(sby_final_predictors) <- NULL
  return(sby_final_predictors)
}
####
## Fim
#
