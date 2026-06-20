#' Validar tamanho de bloco para consultas KNN
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_knn_query_chunk_size Tamanho de bloco solicitado para consultas KNN
#'
#' @return Tamanho de bloco inteiro positivo validado
#'
#' @noRd
sby_validate_knn_query_chunk_size <- function(sby_knn_query_chunk_size){
  
  # Define tamanho padrao quando o valor informado e nulo
  if(is.null(sby_knn_query_chunk_size)){

    # Usa tamanho padrao conservador para consultas em blocos
    sby_knn_query_chunk_size <- 1000L
  }

  # Verifica se o tamanho de bloco e um escalar numerico dentro do limite inteiro
  if(!(is.numeric(sby_knn_query_chunk_size) && length(sby_knn_query_chunk_size) == 1L && !is.na(sby_knn_query_chunk_size) && is.finite(sby_knn_query_chunk_size) && sby_knn_query_chunk_size >= 1L && sby_knn_query_chunk_size == floor(sby_knn_query_chunk_size) && sby_knn_query_chunk_size <= .Machine$integer.max)){

    # Aborta quando o tamanho de bloco KNN e invalido
    sby_adanear_abort(
      sby_message = "'sby_knn_query_chunk_size' deve ser inteiro positivo"
    )
  }

  # Retorna tamanho de bloco normalizado como inteiro
  return(as.integer(floor(sby_knn_query_chunk_size)))
}
####
## Fim
#
