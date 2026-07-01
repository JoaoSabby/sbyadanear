#' Calcular quantidade retida da maioria
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_target_factor Fator binario com as classes observadas
#'
#' @param sby_ratio_under Multiplicador positivo da quantidade minoritaria para definir a maioria retida
#'
#' @param sby_minority_label Rotulo fixo opcional da classe minoritaria
#'
#' @param sby_majority_label Rotulo fixo opcional da classe majoritaria
#'
#' @return Quantidade inteira de linhas majoritarias a manter
#'
#' @noRd
sby_compute_majority_retention_count <- function(
  sby_target_factor,
  sby_ratio_under,
  sby_minority_label = NULL,
  sby_majority_label = NULL
){
  
  # Verifica se a razao alvo e positiva
  if(!(is.numeric(sby_ratio_under) && length(sby_ratio_under) == 1L && !is.na(sby_ratio_under) && sby_ratio_under > 0)){

    # Aborta quando a razao de undersampling e invalida
    sby_adanear_abort(
      sby_message = "'sby_ratio_under' deve ser escalar numerico maior que zero"
    )
  }

  # Identifica papeis de classe para calcular os tamanhos atuais
  if(is.null(sby_minority_label) && is.null(sby_majority_label)){
    sby_class_roles <- sby_binary_class_counts_fast(
      sby_y_vector = sby_target_factor
    )
  }else{
    sby_class_roles <- sby_get_binary_class_roles(
      sby_target_factor    = sby_target_factor,
      sby_minority_label   = sby_minority_label,
      sby_majority_label   = sby_majority_label
    )
    sby_class_roles$sby_minority_count <- as.integer(sby_class_roles$sby_class_counts[sby_class_roles$sby_minority_label])
    sby_class_roles$sby_majority_count <- as.integer(sby_class_roles$sby_class_counts[sby_class_roles$sby_majority_label])
  }
  sby_minority_count <- sby_class_roles$sby_minority_count
  sby_majority_count <- sby_class_roles$sby_majority_count

  # Interpreta sby_ratio_under como multiplicador da quantidade minoritaria
  # final disponivel para o NearMiss-1. Exemplo: 1.0 iguala maioria a
  # minoria; 0.5 retem metade da minoria; 2.0 retem ate duas vezes a
  # minoria, limitado a maioria disponivel.
  sby_target_majority_count <- floor(sby_minority_count * sby_ratio_under)

  # Limita a retencao ao total majoritario disponivel para evitar indices ausentes
  sby_retained_count <- min(
    sby_majority_count,
    sby_target_majority_count
  )

  # Verifica se a configuracao reteve ao menos uma linha majoritaria
  if(sby_retained_count < 1L){

    # Aborta quando a configuracao removeria toda a classe majoritaria
    sby_adanear_abort(
      sby_message = "'sby_ratio_under' reteve zero linhas"
    )
  }

  # Retorna quantidade majoritaria retida pelo criterio configurado
  return(as.integer(sby_retained_count))
}
####
## Fim
#
