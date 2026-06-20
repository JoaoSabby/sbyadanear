#' Preparar etapa de balanceamento ADASYN e NearMiss
#'
#' @param x Objeto de etapa `sby_step_adanear` nao treinado
#'
#' @param training Dados de treinamento da recipe
#'
#' @param info Metadados da recipe
#'
#' @param ... Argumentos adicionais preservados para compatibilidade S3
#'
#' @return Objeto de etapa `sby_step_adanear` treinado
#'
#' @export
prep.step_sby_step_adanear <- function(x, training, info = NULL, ...){
  return(sby_prep_step_sampling(
    x = x,
    training = training,
    info = info,
    sby_step_name = "sby_step_adanear()"
  ))
}
####
## Fim
#
