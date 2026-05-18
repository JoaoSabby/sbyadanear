#' Verificar orcamento de memoria para matrizes densas
#'
#' @param sby_row_count Numero de linhas esperado
#' @param sby_column_count Numero de colunas esperado
#' @param sby_extra_matrix_count Multiplicador de matrizes densas simultaneas
#' @param sby_max_dense_gb Limite maximo em GiB
#' @param sby_context Contexto usado na mensagem
#'
#' @return Invisivelmente TRUE quando dentro do orcamento
#' @noRd
sby_check_dense_memory_budget <- function(
  sby_row_count,
  sby_column_count,
  sby_extra_matrix_count = 1L,
  sby_max_dense_gb = Inf,
  sby_context = "balanceamento"
){
  if(!is.finite(sby_max_dense_gb)){
    return(invisible(TRUE))
  }

  sby_estimated_bytes <- as.numeric(sby_row_count) * as.numeric(sby_column_count) * 8 * as.numeric(sby_extra_matrix_count)
  sby_limit_bytes <- as.numeric(sby_max_dense_gb) * 1024^3

  if(is.na(sby_estimated_bytes) || sby_estimated_bytes > sby_limit_bytes){
    sby_adanear_abort(
      sby_message = paste0(
        "Orcamento de memoria densa excedido em ", sby_context,
        ": estimativa ", format(round(sby_estimated_bytes / 1024^3, 3), scientific = FALSE),
        " GiB > limite ", format(round(sby_max_dense_gb, 3), scientific = FALSE), " GiB"
      )
    )
  }

  return(invisible(TRUE))
}
