#' Carregar dependencias do fluxo de sampling
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @return Retorna invisivelmente TRUE apos validar dependencias
#'
#' @noRd
sby_adanear_load_packages <- function(){
  
  # Evita recarregar dependencias ja validadas na sessao
  if(isTRUE(sby_adanear_state$sby_packages_loaded)){

    # Retorna sucesso quando o estado global ja indica dependencias carregadas
    return(invisible(TRUE))
  }

  # Define pacotes obrigatorios para o fluxo de sampling.
  # Rfast e Suggests (so usado em fallback R puro quando o kernel C nao
  # esta disponivel) e por isso nao entra na lista obrigatoria.
  sby_package_names <- c(
    "cli"
  )

  # Valida disponibilidade de cada pacote obrigatorio
  for(sby_package_name in sby_package_names){
    # Verifica se o pacote corrente esta instalado
    if(!requireNamespace(
      package = sby_package_name,
      quietly = TRUE
    )){

      # Aborta informando a dependencia ausente
      sby_adanear_abort(
        sby_message = paste0(
          "Pacote necessario nao encontrado: ",
          sby_package_name
        )
      )
    }
  }

  # Registra que as dependencias foram validadas na sessao
  sby_adanear_state$sby_packages_loaded <- TRUE

  # Retorna sucesso invisivel apos carregamento das dependencias
  return(invisible(TRUE))
}
####
## Fim
#
