#' Calcular matriz de covariância com rotina otimizada
#'
#' @param sby_x_matrix Matriz numérica de entrada.
#'
#' @return Matriz de covariância calculada por `coop::covar()`.
#' @noRd
sby_fast_covariance_matrix <- function(sby_x_matrix){
  sby_x_matrix <- sby_adanear_as_numeric_matrix(
    sby_predictor_data = sby_x_matrix
  )
  coop::covar(sby_x_matrix)
}
