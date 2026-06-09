#' Calcular parametros de z-score
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_x_matrix Matriz numerica usada para estimar centro e escala
#'
#' @return Lista com vetores numericos `centers` e `scales`
#' @noRd
sby_compute_z_score_params <- function(
  sby_x_matrix,
  sby_engine = c("auto", "native", "FNN", "RcppHNSW", "KernelKnn", "bigKNN")
){
  sby_engine <- match.arg(sby_engine)
  
  # Normaliza entrada para matriz numerica de precisao dupla
  sby_x_matrix <- sby_adanear_as_numeric_matrix(
    sby_predictor_data = sby_x_matrix
  )

  # Calcula parametros pela engine nativa Fortran
  if(sby_engine == "native" && sby_native_symbol_available("compute_zscore_population_fortran_c")){
    sby_configure_blas_threads(sby_workers = 1L)
    sby_params_f <- sby_call_native(
      "compute_zscore_population_fortran_c",
      sby_x_matrix
    )
    sby_params <- list(
      centers = sby_params_f$means,
      scales = sby_params_f$sds
    )
  } else if(sby_adanear_native_available()){
    # Estima centros e escalas por chamada nativa em C
    sby_params <- sby_call_native(
      "compute_z_score_params_c",
      sby_x_matrix
    )
  } else {
    sby_params <- list(
      centers = Rfast::colmeans(sby_x_matrix),
      scales = Rfast::colVars(
        x = sby_x_matrix,
        std = TRUE
      )
    )
  }

  # Identifica escalas ausentes, infinitas ou nao positivas
  sby_invalid <- is.na(sby_params$scales) | !is.finite(sby_params$scales) | sby_params$scales <= 0

  # Verifica se todas as colunas possuem desvio padrao valido
  if(any(sby_invalid)){

    # Captura nomes de colunas para mensagem de erro diagnostica
    sby_column_names <- colnames(sby_x_matrix)

    # Gera nomes padronizados quando a matriz nao possui nomes
    if(is.null(sby_column_names)){

      # Cria nomes sequenciais para as colunas invalidas
      sby_column_names <- paste0(
        "V",
        seq_len(collapse::fncol(sby_x_matrix))
      )
    }

    # Aborta informando colunas com escala indefinida
    sby_adanear_abort(
      sby_message = paste0(
        "Colunas com desvio padrao zero ou indefinido: ",
        paste(
          sby_column_names[sby_invalid],
          collapse = ", "
        )
      )
    )
  }

  # Retorna parametros de escala em formato numerico simples
  return(list(
    centers = as.numeric(sby_params$centers),
    scales  = as.numeric(sby_params$scales)
  ))
}
####
## Fim
#
