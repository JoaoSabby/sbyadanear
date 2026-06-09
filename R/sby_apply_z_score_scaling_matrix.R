#' Aplicar z-score em matrix double
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_x_matrix Matriz numerica a ser padronizada
#' @param sby_scaling_info Lista com centros e escalas de padronizacao
#'
#' @return Matriz numerica padronizada por z-score
#' @noRd
sby_apply_z_score_scaling_matrix <- function(
  sby_x_matrix, 
  sby_scaling_info,
  sby_engine = c("auto", "native", "FNN", "RcppHNSW", "KernelKnn", "bigKNN")
){
  sby_engine <- match.arg(sby_engine)

  # Normaliza entrada para matriz numerica de precisao dupla
  sby_x_matrix <- sby_adanear_as_numeric_matrix(
    sby_predictor_data = sby_x_matrix
  )

  # Valida parametros de escala contra a quantidade de colunas
  sby_validate_scaling_info(
    sby_scaling_info           = sby_scaling_info,
    sby_predictor_column_count = collapse::fncol(sby_x_matrix)
  )

  # Aplica implementacao nativa Fortran se solicitada
  if(sby_engine == "native" && sby_native_symbol_available("apply_zscore_fortran_c")){
    sby_configure_blas_threads(sby_workers = 1L)
    sby_scaled <- .Call(
      apply_zscore_fortran_c,
      sby_x_matrix,
      as.numeric(sby_scaling_info$centers),
      as.numeric(sby_scaling_info$scales)
    )
  } else if(sby_adanear_native_available()){
    # Aplica implementacao nativa C
    sby_scaled <- .Call(
      apply_z_score_c,
      sby_x_matrix,
      as.numeric(sby_scaling_info$centers),
      as.numeric(sby_scaling_info$scales),
      FALSE
    )
  } else {

    sby_centered <- Rfast::eachrow(
      x = sby_x_matrix,
      y = sby_scaling_info$centers,
      oper = "-"
    )
    sby_scaled <- Rfast::eachrow(
      x = sby_centered,
      y = sby_scaling_info$scales,
      oper = "/"
    )
  }

  # Garante armazenamento numerico double apos o calculo
  storage.mode(sby_scaled) <- "double"

  # Preserva nomes de colunas da matriz original
  colnames(sby_scaled) <- colnames(sby_x_matrix)

  # Retorna matriz padronizada
  return(sby_scaled)
}
####
## Fim
#
