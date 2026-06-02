#' Reverter z-score em matrix double
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_x_matrix Matriz numerica padronizada a ser revertida
#' @param sby_scaling_info Lista com centros e escalas de padronizacao
#'
#' @return Matriz numerica restaurada para a escala original
#' @noRd
sby_revert_z_score_scaling_matrix <- function(sby_x_matrix, sby_scaling_info){
  
  # Normaliza entrada para matriz numerica de precisao dupla
  sby_x_matrix <- sby_adanear_as_numeric_matrix(
    sby_predictor_data = sby_x_matrix
  )

  # Valida parametros de escala contra a quantidade de colunas
  sby_validate_scaling_info(
    sby_scaling_info           = sby_scaling_info,
    sby_predictor_column_count = collapse::fncol(sby_x_matrix)
  )

  # Aplica implementacao nativa quando disponivel
  if(sby_adanear_native_available()){

    # Reverte z-score por chamada nativa registrada no pacote
    sby_restored <- .Call(
      apply_z_score_c,
      sby_x_matrix,
      as.numeric(sby_scaling_info$centers),
      as.numeric(sby_scaling_info$scales),
      TRUE
    )
  }else{

    sby_unscaled <- Rfast::eachrow(
      x = sby_x_matrix,
      y = sby_scaling_info$scales,
      oper = "*"
    )
    sby_restored <- Rfast::eachrow(
      x = sby_unscaled,
      y = sby_scaling_info$centers,
      oper = "+"
    )
  }

  # Garante armazenamento numerico double apos a reversao
  storage.mode(sby_restored) <- "double"

  # Preserva nomes de colunas da matriz original
  colnames(sby_restored) <- colnames(sby_x_matrix)

  # Retorna matriz restaurada
  return(sby_restored)
}
####
## Fim
#
