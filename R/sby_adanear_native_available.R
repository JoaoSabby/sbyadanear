#' Verificar disponibilidade das rotinas nativas
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @return Valor logico indicando se rotinas nativas podem ser chamadas
#' @noRd
sby_adanear_native_available <- function(){
  
  # Verifica se a DLL do pacote esta registrada na sessao atual
  return(is.loaded(
    symbol = "OU_ApplyZScoreC",
    PACKAGE = "sbyadanear"
  ))
}
####
## Fim
#
