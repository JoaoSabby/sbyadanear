#' Normalizar componente de resultado KNN para matriz
#'
#' @param sby_component Componente de indices ou distancias retornado por engine KNN
#'
#' @param sby_k Numero de vizinhos esperado
#'
#' @param sby_storage Tipo de armazenamento de saida
#'
#' @return Matriz com uma linha por consulta e `sby_k` colunas
#'
#' @noRd
sby_as_knn_matrix <- function(sby_component, sby_k, sby_storage = c("integer", "double")){
  sby_storage <- match.arg(sby_storage)

  if(is.matrix(sby_component)){
    sby_out <- sby_component
  }else if(is.data.frame(sby_component)){
    sby_out <- as.matrix(sby_component)
  }else if(is.list(sby_component)){
    sby_out <- matrix(
      data = unlist(sby_component, use.names = FALSE),
      ncol = sby_k,
      byrow = TRUE
    )
  }else if(!is.null(dim(sby_component))){
    sby_out <- as.matrix(sby_component)
  }else{
    sby_out <- matrix(
      data = sby_component,
      ncol = sby_k
    )
  }

  if(!identical(collapse::fncol(sby_out), as.integer(sby_k))){
    sby_adanear_abort(
      sby_message = "Resultado KNN retornou numero inesperado de colunas"
    )
  }

  if(identical(sby_storage, "integer")){
    storage.mode(sby_out) <- "integer"
  }else{
    storage.mode(sby_out) <- "double"
  }
  return(sby_out)
}
