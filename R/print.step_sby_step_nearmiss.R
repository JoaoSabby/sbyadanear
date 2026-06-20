#' Imprimir etapa de balanceamento NearMiss
#'
#' @param x Objeto de etapa
#'
#' @param width Largura de exibicao usada pelo metodo de impressao
#'
#' @param ... Argumentos adicionais preservados para compatibilidade S3
#'
#' @return Retorna invisivelmente o objeto de etapa informado
#'
#' @export
print.step_sby_step_nearmiss <- function(x, width = max(20, options()$width - 30), ...){
  return(sby_print_step_sampling(
    x = x,
    width = width,
    sby_title = "Balanceamento NearMiss-1 usando "
  ))
}
####
## Fim
#
