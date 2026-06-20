#' Inferir tipos numericos para restauracao posterior
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_data_frame Dados preditores numericos usados para inferencia de tipos
#'
#' @return Data frame com nomes de colunas e tipos numericos inferidos
#'
#' @noRd
sby_infer_numeric_column_types <- function(sby_data_frame){
  
  # Captura nomes de colunas associados aos preditores
  sby_column_names <- sby_adanear_get_column_names(
    sby_predictor_data = sby_data_frame
  )

  # Mantem a semantica do tipo original da coluna. Uma coluna originalmente
  # integer deve voltar como integer; uma coluna originalmente double deve
  # voltar como double mesmo quando seus valores observados parecem inteiros.
  sby_data_frame <- as.data.frame(
    x = sby_data_frame,
    stringsAsFactors = FALSE
  )

  sby_inferred_type <- vapply(
    X = sby_data_frame,
    FUN = function(sby_column_data){
      if(is.integer(sby_column_data)){
        return("integer")
      }
      return("double")
    },
    FUN.VALUE = character(1L)
  )

  sby_integer_min <- vapply(
    X = sby_data_frame,
    FUN = function(sby_column_data){
      if(is.integer(sby_column_data)){
        return(min(sby_column_data, na.rm = TRUE))
      }
      return(NA_real_)
    },
    FUN.VALUE = numeric(1L)
  )

  sby_integer_max <- vapply(
    X = sby_data_frame,
    FUN = function(sby_column_data){
      if(is.integer(sby_column_data)){
        return(max(sby_column_data, na.rm = TRUE))
      }
      return(NA_real_)
    },
    FUN.VALUE = numeric(1L)
  )

  # Retorna metadados de tipos inferidos por coluna, incluindo os limites
  # originais necessarios para truncar linhas sinteticas inteiras apos a
  # reversao do z-score.
  return(data.frame(
    sby_column_name = sby_column_names,
    sby_inferred_type = unname(sby_inferred_type),
    sby_integer_min = unname(sby_integer_min),
    sby_integer_max = unname(sby_integer_max),
    stringsAsFactors = FALSE
  ))
}
####
## Fim
#
