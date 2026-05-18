#' Calcular quantidade sintetica da minoria
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_target_factor Fator binario com as classes observadas
#' @param sby_over_ratio Fator relativo de expansao da classe minoritaria
#'
#' @return Quantidade inteira de linhas sinteticas a gerar
#' @noRd
sby_compute_minority_expansion_count <- function(sby_target_factor, sby_over_ratio){
  
  # Verifica se o fator de expansao e um escalar numerico positivo
  if(!(is.numeric(sby_over_ratio) && length(sby_over_ratio) == 1L && !is.na(sby_over_ratio) && sby_over_ratio > 0)){

    # Aborta quando o fator de oversampling e invalido
    sby_adanear_abort(
      sby_message = "'sby_over_ratio' deve ser escalar numerico positivo"
    )
  }

  # Identifica papeis de classe por contagem leve baseada em tabulate
  sby_class_roles <- sby_binary_class_counts_fast(
    sby_y_vector = sby_target_factor
  )
  sby_minority_count <- sby_class_roles$sby_minority_count
  sby_synthetic_count <- max(
    1L,
    as.integer(floor(sby_minority_count * sby_over_ratio))
  )

  # Retorna quantidade de amostras sinteticas a gerar
  return(sby_synthetic_count)
}
####
## Fim
#
