#' Decidir se a rota "native" deve usar o atalho HPC consolidado
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato
#' de entrada explicito e retorno controlado. O atalho HPC substitui a rota
#' "native" como caminho rapido quando o motor consolidado esta disponivel e o
#' contrato simples de saida e suficiente. As rotinas originais continuam
#' acessiveis: o atalho so e tomado quando o chamador escolhe explicitamente
#' `sby_knn_engine = "native"`, sem auditoria, sem restauro de tipos, sem escala
#' intermediaria e com metrica euclidiana. Em qualquer outro caso, o fluxo
#' classico segue inalterado.
#'
#' @return Valor logico indicando se o atalho HPC deve ser tomado.
#' @noRd
sby_should_route_native_to_hpc <- function(
  sby_knn_engine,
  sby_knn_distance_metric,
  sby_audit,
  sby_restore_types = FALSE,
  sby_return_scaled = FALSE
){
  if(!identical(sby_knn_engine, "native")){
    return(FALSE)
  }
  if(!identical(sby_knn_distance_metric, "euclidean")){
    return(FALSE)
  }
  if(isTRUE(sby_audit) || isTRUE(sby_restore_types) || isTRUE(sby_return_scaled)){
    return(FALSE)
  }
  return(sby_adanear_hpc_available())
}
####
## Fim
#
