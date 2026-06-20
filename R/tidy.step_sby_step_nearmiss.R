#' Organizar metadados da etapa de balanceamento NearMiss
#'
#' @param x Objeto de etapa
#'
#' @param ... Argumentos adicionais preservados para compatibilidade S3
#'
#' @return Data frame com metadados principais da etapa
#'
#' @export
tidy.step_sby_step_nearmiss <- function(x, ...){
  return(sby_tidy_step_sampling(
    x = x
  ))
}
####
## Fim
#
