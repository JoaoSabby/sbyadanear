#' Validar parametros do engine RcppHNSW
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_knn_hnsw_m Conectividade solicitada para o indice HNSW
#'
#' @param sby_knn_hnsw_ef Tamanho da lista dinamica HNSW solicitada
#'
#' @return Lista com parametros HNSW convertidos para inteiros
#'
#' @noRd
sby_validate_hnsw_params <- function(sby_knn_hnsw_m, sby_knn_hnsw_ef){
  
  # Verifica se a conectividade HNSW e um escalar numerico finito valido
  if(!(is.numeric(sby_knn_hnsw_m) && length(sby_knn_hnsw_m) == 1L && !is.na(sby_knn_hnsw_m) && is.finite(sby_knn_hnsw_m) && sby_knn_hnsw_m >= 2L && sby_knn_hnsw_m == floor(sby_knn_hnsw_m) && sby_knn_hnsw_m <= .Machine$integer.max)){

    # Aborta quando a conectividade HNSW esta fora do dominio permitido
    sby_adanear_abort(
      sby_message = "'sby_knn_hnsw_m' deve ser inteiro >= 2"
    )
  }

  # Verifica se a lista dinamica HNSW e um escalar numerico positivo finito
  if(!(is.numeric(sby_knn_hnsw_ef) && length(sby_knn_hnsw_ef) == 1L && !is.na(sby_knn_hnsw_ef) && is.finite(sby_knn_hnsw_ef) && sby_knn_hnsw_ef >= 1L && sby_knn_hnsw_ef == floor(sby_knn_hnsw_ef) && sby_knn_hnsw_ef <= .Machine$integer.max)){

    # Aborta quando o parametro de busca HNSW esta fora do dominio permitido
    sby_adanear_abort(
      sby_message = "'sby_knn_hnsw_ef' deve ser inteiro positivo"
    )
  }

  # Retorna parametros HNSW normalizados como inteiros
  return(list(
    sby_knn_hnsw_m  = as.integer(sby_knn_hnsw_m),
    sby_knn_hnsw_ef = as.integer(sby_knn_hnsw_ef)
  ))
}
####
## Fim
#
