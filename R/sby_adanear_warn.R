#' Emitir aviso padronizado
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_message Mensagem tecnica a ser exibida no aviso
#'
#' @return Retorna invisivelmente NULL apos emitir o aviso
#'
#' @noRd
sby_adanear_warn <- function(sby_message){
  # Emite aviso estruturado quando o pacote cli esta disponivel
  if(requireNamespace(
    package = "cli",
    quietly = TRUE
  )){

    # Emite aviso por meio do sistema cli
    cli::cli_warn(
      message = sby_message
    )

    # Retorna sem produzir saida visivel apos avisar via cli
    return(invisible(NULL))
  }

  # Emite aviso base quando cli nao esta disponivel
  warning(
    sby_message,
    call. = FALSE
  )

  # Retorna sem produzir saida visivel apos avisar via base
  return(invisible(NULL))
}
####
## Fim
#
