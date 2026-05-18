#' Montar preditores preservando linhas originais observadas
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
  sby_original_row_count <- nrow(sby_original_predictor_data)
  sby_is_original <- sby_retained_index <= sby_original_row_count

  sby_final_rows <- vector("list", length(sby_retained_index))

  if(any(!sby_is_original)){
    sby_synthetic_scaled <- sby_final_scaled[!sby_is_original, , drop = FALSE]
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
  }else{
    sby_synthetic_predictors <- sby_original_predictor_data[0L, , drop = FALSE]
  }

  sby_synthetic_cursor <- 1L
  for(i in seq_along(sby_retained_index)){
    if(isTRUE(sby_is_original[[i]])){
      sby_final_rows[[i]] <- sby_original_predictor_data[sby_retained_index[[i]], , drop = FALSE]
    }else{
      sby_final_rows[[i]] <- sby_synthetic_predictors[sby_synthetic_cursor, , drop = FALSE]
      sby_synthetic_cursor <- sby_synthetic_cursor + 1L
    }
  }

  sby_final_predictors <- do.call(rbind, sby_final_rows)
  rownames(sby_final_predictors) <- NULL
  return(sby_final_predictors)
}
####
## Fim
#
