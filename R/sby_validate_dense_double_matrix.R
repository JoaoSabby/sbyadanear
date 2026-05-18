#' Validar matriz densa double
#'
#' @param sby_x_matrix Objeto de entrada esperado como matrix double densa
#' @param sby_name Nome usado em mensagens de erro
#' @param sby_allow_integer Permite converter matrix integer para double uma unica vez
#' @param sby_allow_na Permite valores ausentes
#' @param sby_allow_infinite Permite valores infinitos
#'
#' @return A propria matriz quando ja estiver no contrato, ou uma copia double quando integer for permitido
#' @noRd
sby_validate_dense_double_matrix <- function(
  sby_x_matrix,
  sby_name = "sby_x_matrix",
  sby_allow_integer = FALSE,
  sby_allow_na = FALSE,
  sby_allow_infinite = FALSE
){
  if(sby_is_sparse_matrix(sby_x_matrix)){
    sby_adanear_abort(
      sby_message = paste0("'", sby_name, "' deve ser matrix double densa; matrizes esparsas devem ser tratadas antes da chamada")
    )
  }

  if(!is.matrix(sby_x_matrix)){
    sby_adanear_abort(
      sby_message = paste0("'", sby_name, "' deve ser uma matrix double densa")
    )
  }

  if(nrow(sby_x_matrix) < 2L){
    sby_adanear_abort(
      sby_message = paste0("'", sby_name, "' deve conter ao menos duas linhas")
    )
  }

  if(ncol(sby_x_matrix) < 1L){
    sby_adanear_abort(
      sby_message = paste0("'", sby_name, "' deve conter ao menos uma coluna")
    )
  }

  sby_matrix_type <- typeof(sby_x_matrix)
  if(!identical(sby_matrix_type, "double")){
    if(isTRUE(sby_allow_integer) && identical(sby_matrix_type, "integer")){
      sby_x_matrix <- matrix(
        data = as.double(sby_x_matrix),
        nrow = nrow(sby_x_matrix),
        ncol = ncol(sby_x_matrix),
        dimnames = dimnames(sby_x_matrix)
      )
    }else{
      sby_adanear_abort(
        sby_message = paste0("'", sby_name, "' deve ter storage mode double")
      )
    }
  }

  if(!isTRUE(sby_allow_na) && anyNA(sby_x_matrix)){
    sby_adanear_abort(
      sby_message = paste0("'", sby_name, "' nao pode conter NA")
    )
  }

  if(!isTRUE(sby_allow_infinite) && any(!is.finite(sby_x_matrix))){
    sby_adanear_abort(
      sby_message = paste0("'", sby_name, "' nao pode conter Inf ou -Inf")
    )
  }

  return(sby_x_matrix)
}
