#' Validar alinhamento de objeto matricial sem densificar
#'
#' @param sby_x Objeto com dimensao de linhas
#' @param sby_y_vector Vetor alvo associado
#'
#' @return Invisivelmente TRUE quando alinhado
#' @noRd
sby_validate_matrix_like_row_count <- function(sby_x, sby_y_vector){
  sby_row_count <- tryCatch(
    collapse::fnrow(sby_x),
    error = function(e) NULL
  )
  if(is.null(sby_row_count) || length(sby_row_count) != 1L || is.na(sby_row_count)){
    sby_adanear_abort("'sby_x_matrix' deve possuir numero de linhas valido")
  }
  if(length(sby_y_vector) != sby_row_count){
    sby_adanear_abort("'sby_y_vector' deve ter comprimento igual ao numero de linhas de 'sby_x_matrix'")
  }
  return(invisible(TRUE))
}
####
## Fim
#
