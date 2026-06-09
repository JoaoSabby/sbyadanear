#' Verificar disponibilidade do motor HPC consolidado
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato
#' de entrada explicito e retorno controlado. O atalho HPC depende do simbolo
#' nativo consolidado sby_adanear_hpc_result_cpp. Quando ele nao esta disponivel, o
#' chamador deve recorrer transparentemente a rota classica do pacote para que
#' tudo continue funcionando sem o motor HPC compilado.
#'
#' @return Valor logico indicando se o motor HPC consolidado pode ser chamado.
#' @noRd
sby_adanear_hpc_available <- function(){
  return(sby_native_symbol_available("sby_adanear_hpc_result_cpp"))
}
####
## Fim
#
