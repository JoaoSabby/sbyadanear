#' Emitir mensagem informativa padronizada
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_message Mensagem tecnica a ser exibida ao usuario
#'
#' @return Retorna invisivelmente NULL apos emitir a mensagem
#'
#' @noRd
sby_adanear_inform <- function(sby_message){
  
  # Emite mensagem estruturada quando o pacote cli esta disponivel
  if(requireNamespace(
    package = "cli",
    quietly = TRUE
  )){

    # Informa a mensagem por meio do sistema cli
    cli::cli_inform(
      message = sby_message
    )

    # Retorna sem produzir saida visivel apos informar via cli
    return(invisible(NULL))
  }

  # Emite mensagem base quando cli nao esta disponivel
  message(
    sby_message
  )

  # Retorna sem produzir saida visivel apos informar via base
  return(invisible(NULL))
}
####
## Fim
#
