#' Validar numero de workers KNN
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_knn_workers Numero de workers solicitado para o engine KNN
#'
#' @return Numero inteiro positivo de workers validado
#'
#' @noRd
sby_validate_knn_workers <- function(sby_knn_workers){
  
  # Verifica se a quantidade de workers e um escalar numerico positivo finito
  if(!(is.numeric(sby_knn_workers) && length(sby_knn_workers) == 1L && !is.na(sby_knn_workers) && is.finite(sby_knn_workers) && sby_knn_workers >= 1L && sby_knn_workers == floor(sby_knn_workers) && sby_knn_workers <= .Machine$integer.max)){

    # Aborta quando a configuracao de workers nao e valida
    sby_adanear_abort(
      sby_message = "'sby_knn_workers' deve ser inteiro positivo"
    )
  }

  # Retorna a quantidade de workers como inteiro
  return(as.integer(sby_knn_workers))
}
####
## Fim
#
