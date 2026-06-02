#' Inferir tipos numericos para restauracao posterior
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_data_frame Dados preditores numericos usados para inferencia de tipos
#'
#' @return Data frame com nomes de colunas e tipos numericos inferidos
#' @noRd
sby_infer_numeric_column_types <- function(sby_data_frame){
  
  # Converte dados de entrada para matriz numerica padronizada
  sby_x_matrix <- sby_adanear_as_numeric_matrix(
    sby_predictor_data = sby_data_frame
  )

  # Captura nomes de colunas associados aos preditores
  sby_column_names <- sby_adanear_get_column_names(
    sby_predictor_data = sby_data_frame
  )

  #' Inferir tipo numerico de uma coluna isolada
  #'
  #' @details
  #' A funcao local classifica valores numericos em categorias discretas ou continuas para apoiar restauracao posterior de tipos
  #'
  #' @param sby_column_data Vetor numerico correspondente a uma coluna preditora
  #'
  #' @return Marcador textual do tipo numerico inferido
  #' @noRd
  sby_infer_one <- function(sby_column_data){
    # Obtem valores unicos ordenados para classificar colunas binarias
    sby_unique_values <- sort(
      x = unique(sby_column_data)
    )

    # Identifica coluna binaria codificada numericamente
    if(length(sby_unique_values) <= 2L && all(sby_unique_values %in% c(0, 1))){

      # Retorna marcador de tipo binario
      return("binary")
    }

    # Verifica se todos os valores sao equivalentes a inteiros
    sby_is_integer_like <- all(abs(sby_column_data - round(sby_column_data)) < sqrt(.Machine$double.eps))

    # Identifica coluna numerica com semantica inteira
    if(sby_is_integer_like){

      # Retorna marcador de tipo inteiro
      return("integer")
    }

    # Retorna marcador de tipo numerico continuo
    return("double")
  }

  # Retorna metadados de tipos inferidos por coluna
  return(data.frame(
    sby_column_name = sby_column_names,
    sby_inferred_type = vapply(
      X = seq_len(collapse::fncol(sby_x_matrix)),
      FUN = function(j){
        # Infere o tipo da coluna corrente da matriz numerica
        return(sby_infer_one(
          sby_column_data = sby_x_matrix[, j]
        ))
      },
      FUN.VALUE = character(1L)
    ),
    stringsAsFactors = FALSE
  ))
}
####
## Fim
#
