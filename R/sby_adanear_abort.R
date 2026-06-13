#' Emitir erro padronizado com cli
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_message Mensagem tecnica a ser exibida no erro
#'
#' @return Interrompe a execucao com erro padronizado
#' @noRd
sby_adanear_abort <- function(sby_message){
  
  # Emite erro estruturado quando o pacote cli esta disponivel
  if(requireNamespace(
    package = "cli",
    quietly = TRUE
  )){

    # Interrompe a execucao sem incluir a chamada original
    cli::cli_abort(
      message = sby_message,
      call = NULL
    )
  }

  # Interrompe a execucao usando erro base como contingencia
  stop(
    sby_message,
    call. = FALSE
  )
}
####
## Fim
#
