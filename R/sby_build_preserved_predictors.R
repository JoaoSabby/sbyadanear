#' Reconstruir preditores preservando linhas originais quando possivel
#'
#' @noRd
sby_build_preserved_predictors <- function(
  sby_original_predictor_data,
  sby_retained_index,
  sby_final_scaled,
  sby_scaling_info,
  sby_type_info,
  sby_restore_types = TRUE
){
  sby_final_original <- sby_revert_z_score_scaling_matrix(
    sby_x_matrix = sby_final_scaled,
    sby_scaling_info = sby_scaling_info
  )
  sby_final_predictors <- as.data.frame(
    sby_final_original,
    stringsAsFactors = FALSE
  )
  names(sby_final_predictors) <- names(sby_original_predictor_data)

  if(isTRUE(sby_restore_types)){
    sby_final_predictors <- sby_restore_numeric_column_types(
      sby_x_matrix = sby_final_predictors,
      sby_type_info = sby_type_info,
      sby_as_data_frame = TRUE
    )
  }

  sby_original_n <- collapse::fnrow(sby_original_predictor_data)
  sby_original_positions <- which(sby_retained_index <= sby_original_n)
  if(length(sby_original_positions) > 0L){
    sby_final_predictors[sby_original_positions, ] <- sby_original_predictor_data[
      sby_retained_index[sby_original_positions], , drop = FALSE
    ]
  }
  rownames(sby_final_predictors) <- NULL
  sby_final_predictors
}
####
## Fim
#
