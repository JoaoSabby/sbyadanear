#' Preparar etapa de balanceamento NearMiss
#'
#' @param x Objeto de etapa `sby_step_nearmiss` nao treinado
#'
#' @param training Dados de treinamento da recipe
#'
#' @param info Metadados da recipe
#'
#' @param ... Argumentos adicionais preservados para compatibilidade S3
#'
#' @return Objeto de etapa treinado
#'
#' @export
prep.step_sby_step_nearmiss <- function(x, training, info = NULL, ...){
  return(sby_prep_step_sampling(
    x = x,
    training = training,
    info = info,
    sby_step_name = "sby_step_nearmiss()"
  ))
}
####
## Fim
#
