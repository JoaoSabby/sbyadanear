#' Validar parametro logico escalar
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_value Valor a ser validado como logico escalar
#'
#' @param sby_name Nome tecnico do parametro usado em mensagens de erro
#'
#' @return Valor logico validado sem alteracoes
#'
#' @noRd
sby_validate_logical_scalar <- function(sby_value, sby_name){
  
  # Verifica se o valor informado representa um escalar logico nao ausente
  if(!(is.logical(sby_value) && length(sby_value) == 1L && !is.na(sby_value))){

    # Aborta com mensagem descritiva para parametro logico invalido
    sby_adanear_abort(
      sby_message = paste0(
        "'",
        sby_name,
        "' deve ser logical escalar nao ausente"
      )
    )
  }

  # Retorna o valor validado sem alteracoes
  return(sby_value)
}
####
## Fim
#
