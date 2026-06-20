#' Restaurar tipos numericos inferidos
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_x_matrix Matriz numerica com valores a restaurar
#'
#' @param sby_type_info Data frame com tipos numericos inferidos por coluna
#'
#' @param sby_as_data_frame Indicador logico para retornar data frame
#'
#' @return Matriz ou data frame com tipos numericos restaurados
#'
#' @noRd
sby_restore_numeric_column_types <- function(sby_x_matrix, sby_type_info, sby_as_data_frame = TRUE){
  
  # Converte dados de entrada para matriz numerica padronizada
  sby_x_matrix <- sby_adanear_as_numeric_matrix(
    sby_predictor_data = sby_x_matrix
  )

  # Verifica compatibilidade entre matriz e metadados de tipos
  if(collapse::fncol(sby_x_matrix) != collapse::fnrow(sby_type_info)){

    # Aborta quando os metadados nao correspondem as colunas da matriz
    sby_adanear_abort(
      sby_message = "Inconsistencia entre numero de colunas e sby_type_info"
    )
  }

  # Ajusta valores numericos conforme o tipo inferido de cada coluna
  for(j in seq_len(collapse::fncol(sby_x_matrix))){
    # Verifica interrupcao periodicamente durante restauracao matricial
    if(j %% 64L == 1L){

      # Executa ponto cooperativo de interrupcao
      sby_adanear_check_user_interrupt()
    }

    # Recupera o tipo inferido para a coluna corrente
    sby_inferred_type <- sby_type_info$sby_inferred_type[[j]]

    # Restaura codificacao binaria quando aplicavel
    if(identical(
      x = sby_inferred_type,
      y = "binary"
    )){

      sby_x_matrix[, j] <- kit::iif(
        test = sby_x_matrix[, j] >= 0.5,
        yes = 1,
        no = 0
      )
    }else if(identical(
      x = sby_inferred_type,
      y = "integer"
    )){

      # Arredonda valores para restaurar semantica inteira e respeita os
      # limites observados na coluna original antes da padronizacao.
      sby_x_matrix[, j] <- round(
        x = sby_x_matrix[, j]
      )

      if("sby_integer_min" %in% names(sby_type_info) &&
         "sby_integer_max" %in% names(sby_type_info) &&
         is.finite(sby_type_info$sby_integer_min[[j]]) &&
         is.finite(sby_type_info$sby_integer_max[[j]])){
        sby_x_matrix[, j] <- pmin(
          pmax(sby_x_matrix[, j], sby_type_info$sby_integer_min[[j]]),
          sby_type_info$sby_integer_max[[j]]
        )
      }
    }
  }

  # Retorna matriz quando data frame nao foi solicitado
  if(!sby_as_data_frame){

    # Entrega matriz restaurada ao chamador
    return(sby_x_matrix)
  }

  # Converte matriz restaurada em data frame nomeado
  sby_out <- as.data.frame(
    x = sby_x_matrix,
    stringsAsFactors = FALSE
  )
  names(sby_out) <- sby_type_info$sby_column_name

  # Converte colunas do data frame para os tipos numericos inferidos
  for(j in seq_len(collapse::fncol(sby_x_matrix))){
    # Verifica interrupcao periodicamente durante conversao tabular
    if(j %% 64L == 1L){

      # Executa ponto cooperativo de interrupcao
      sby_adanear_check_user_interrupt()
    }

    # Recupera o tipo inferido para a coluna corrente
    sby_inferred_type <- sby_type_info$sby_inferred_type[[j]]

    # Converte colunas discretas para inteiro e continuas para numerico
    if(identical(
      x = sby_inferred_type,
      y = "binary"
    ) || identical(
      x = sby_inferred_type,
      y = "integer"
    )){

      # Aplica tipo inteiro para colunas discretas
      sby_out[[j]] <- as.integer(sby_out[[j]])
    }else{

      # Aplica tipo numerico para colunas continuas
      sby_out[[j]] <- as.numeric(sby_out[[j]])
    }
  }

  # Retorna data frame com tipos restaurados
  return(sby_out)
}
####
## Fim
#
