#' Validar parametros de escala precomputados
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_scaling_info Lista com centros e escalas precomputados
#'
#' @param sby_predictor_column_count Quantidade esperada de colunas preditoras
#'
#' @return Retorna invisivelmente TRUE quando os parametros sao validos
#'
#' @noRd
sby_validate_scaling_info <- function(sby_scaling_info, sby_predictor_column_count){
  
  # Verifica se os componentes obrigatorios de escala estao presentes
  if(!(is.list(sby_scaling_info) && !is.null(sby_scaling_info$centers) && !is.null(sby_scaling_info$scales))){

    # Aborta quando a estrutura de escala esta incompleta
    sby_adanear_abort(
      sby_message = "'sby_precomputed_scaling' deve conter 'centers' e 'scales'"
    )
  }

  # Verifica se centros e escalas possuem comprimento compativel com os preditores
  if(length(sby_scaling_info$centers) != sby_predictor_column_count || length(sby_scaling_info$scales) != sby_predictor_column_count){

    # Aborta quando ha incompatibilidade dimensional nos parametros de escala
    sby_adanear_abort(
      sby_message = "'sby_precomputed_scaling' deve ter um centro e uma escala por coluna"
    )
  }

  # Verifica se centros e escalas sao valores numericos finitos e presentes
  if(anyNA(sby_scaling_info$centers) || anyNA(sby_scaling_info$scales) || any(!is.finite(sby_scaling_info$centers)) || any(!is.finite(sby_scaling_info$scales))){

    # Aborta quando parametros de escala contem valores invalidos
    sby_adanear_abort(
      sby_message = "'sby_precomputed_scaling' contem valores ausentes ou infinitos"
    )
  }

  # Verifica se todas as escalas sao estritamente positivas
  if(any(sby_scaling_info$scales <= 0)){

    # Aborta quando alguma escala impossibilita padronizacao valida
    sby_adanear_abort(
      sby_message = "'sby_precomputed_scaling$scales' deve conter apenas valores positivos"
    )
  }

  # Retorna sucesso invisivel apos validacao da estrutura de escala
  return(invisible(TRUE))
}
####
## Fim
#
