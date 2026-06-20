#' Aplicar etapa de balanceamento ADASYN e NearMiss em novos dados
#'
#' @param object Objeto de etapa `sby_step_adanear` treinado
#'
#' @param new_data Dados novos fornecidos ao `bake()` da recipe
#'
#' @param ... Argumentos adicionais preservados para compatibilidade S3
#'
#' @return Tibble balanceado ou lista de auditoria quando configurado
#'
#' @export
bake.step_sby_step_adanear <- function(object, new_data, ...){
  return(sby_bake_step_sampling(
    object = object,
    new_data = new_data,
    sby_step_name = "sby_step_adanear()"
  ))
}
####
## Fim
#
