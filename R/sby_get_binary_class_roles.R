#' Identificar classes minoritaria e majoritaria
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_target_factor Fator binario com as classes observadas
#' @param sby_minority_label Rotulo fixo opcional da classe minoritaria
#' @param sby_majority_label Rotulo fixo opcional da classe majoritaria
#'
#' @return Lista com contagens e rotulos das classes minoritaria e majoritaria
#' @noRd
sby_get_binary_class_roles <- function(
  sby_target_factor,
  sby_minority_label = NULL,
  sby_majority_label = NULL
){
  
  # Calcula distribuicao de frequencias por classe. Usa tabulate ao inves de
  # table porque tabulate e cerca de 5x-10x mais rapido em fatores grandes
  # (evita hash de strings, formatacao e ordenacao alfabetica de levels).
  if(!is.factor(sby_target_factor)){
    sby_target_factor <- as.factor(sby_target_factor)
  }
  sby_levels <- levels(sby_target_factor)
  sby_counts_int <- tabulate(
    bin = as.integer(sby_target_factor),
    nbins = length(sby_levels)
  )
  sby_class_counts <- array(
    data = sby_counts_int,
    dim = length(sby_counts_int),
    dimnames = structure(list(sby_levels), names = "sby_target_factor")
  )
  class(sby_class_counts) <- "table"

  # Verifica se o alvo possui exatamente duas classes
  if(length(sby_class_counts) != 2L){

    # Aborta quando o alvo nao e binario
    sby_adanear_abort(
      sby_message = "'sby_target_vector' deve ser binario"
    )
  }

  # Usa papeis fixos quando o chamador precisa preservar a classe rara original
  if(!is.null(sby_minority_label) || !is.null(sby_majority_label)){

    # Verifica se ambos os rotulos fixos foram informados em conjunto
    if(is.null(sby_minority_label) || is.null(sby_majority_label)){

      # Aborta quando apenas um papel fixo foi informado
      sby_adanear_abort(
        sby_message = "'sby_minority_label' e 'sby_majority_label' devem ser informados em conjunto"
      )
    }

    # Verifica se os rotulos fixos sao escalares validos
    if(length(sby_minority_label) != 1L || length(sby_majority_label) != 1L || is.na(sby_minority_label) || is.na(sby_majority_label)){

      # Aborta quando algum papel fixo nao e escalar ou esta ausente
      sby_adanear_abort(
        sby_message = "Rotulos fixos de classe devem ser escalares nao ausentes"
      )
    }

    # Normaliza rotulos fixos para comparacao com os nomes da tabela
    sby_minority_label <- as.character(sby_minority_label)
    sby_majority_label <- as.character(sby_majority_label)

    # Verifica se os papeis fixos sao distintos
    if(identical(sby_minority_label, sby_majority_label)){

      # Aborta quando os dois papeis apontam para a mesma classe
      sby_adanear_abort(
        sby_message = "'sby_minority_label' e 'sby_majority_label' devem ser distintos"
      )
    }

    # Verifica se os rotulos fixos existem no alvo atual
    if(!all(c(sby_minority_label, sby_majority_label) %in% names(sby_class_counts))){

      # Aborta quando algum papel fixo nao esta presente nos dados atuais
      sby_adanear_abort(
        sby_message = "Rotulos fixos de classe devem existir em 'sby_target_vector'"
      )
    }

    # Retorna papeis fixos com as contagens atuais de cada classe
    return(list(
      sby_class_counts   = sby_class_counts,
      sby_minority_label = sby_minority_label,
      sby_majority_label = sby_majority_label
    ))
  }

  # Verifica se ha desbalanceamento entre as classes
  if(sby_class_counts[[1L]] == sby_class_counts[[2L]]){

    # Aborta quando a rotina de balanceamento nao tem classe majoritaria definida
    sby_adanear_abort(
      sby_message = "As rotinas de sampling requerem classes desbalanceadas"
    )
  }

  # Retorna papeis de classe calculados a partir das frequencias
  return(list(
    sby_class_counts   = sby_class_counts,
    sby_minority_label = names(sby_class_counts)[which.min(sby_class_counts)],
    sby_majority_label = names(sby_class_counts)[which.max(sby_class_counts)]
  ))
}
####
## Fim
#
