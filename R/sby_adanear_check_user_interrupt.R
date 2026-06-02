#' Verificar interrupcao solicitada pelo usuario
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @return Retorna invisivelmente TRUE apos a verificacao de interrupcao
#' @noRd
sby_adanear_check_user_interrupt <- function(){
  
  # Utiliza rotina nativa de interrupcao quando disponivel
  if(sby_adanear_native_available()){

    # Aciona ponto de interrupcao implementado em codigo nativo
    .Call(
      check_user_interrupt_c
    )
  }else{

    # Executa pausa nula para permitir processamento cooperativo de interrupcoes
    Sys.sleep(
      time = 0
    )
  }

  # Retorna sucesso invisivel apos a checagem de interrupcao
  return(invisible(TRUE))
}
####
## Fim
#
