#' Organizar metadados da etapa de balanceamento ADASYN e NearMiss
#'
#' @param x Objeto de etapa `sby_step_adanear`
#'
#' @param ... Argumentos adicionais preservados para compatibilidade S3
#'
#' @return Data frame com metadados principais da etapa
#'
#' @export
tidy.step_sby_step_adanear <- function(x, ...){
  return(sby_tidy_step_sampling(
    x = x
  ))
}
####
## Fim
#
