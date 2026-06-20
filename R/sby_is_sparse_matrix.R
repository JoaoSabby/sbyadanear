#' Identificar matrizes esparsas do pacote Matrix
#'
#' @details
#' A funcao usa apenas inspecao de classe para detectar objetos esparsos sem
#' exigir carregamento do pacote Matrix. Essa verificacao evita conversoes
#' densas silenciosas em rotinas KNN e nativas que esperam matrix double.
#'
#' @param sby_x Objeto a ser inspecionado
#'
#' @return Valor logico indicando se `sby_x` herda de uma classe esparsa Matrix
#'
#' @noRd
sby_is_sparse_matrix <- function(sby_x){
  # Detecta classes esparsas formais do pacote Matrix sem adicionar dependencia forte
  return(inherits(
    x = sby_x,
    what = "sparseMatrix"
  ))
}
####
## Fim
#
