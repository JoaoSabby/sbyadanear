#' Validar entradas de sampling
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_predictor_data Dados preditores em data frame ou matriz
#'
#' @param sby_target_vector Vetor alvo binario associado aos preditores
#'
#' @param sby_seed Semente numerica usada para reprodutibilidade
#'
#' @return Retorna invisivelmente TRUE quando as entradas sao validas
#'
#' @noRd
sby_validate_sampling_inputs <- function(
  sby_predictor_data,
  sby_target_vector,
  sby_seed
){
  
  # Carrega dependencias necessarias para as rotinas de sampling
  sby_adanear_load_packages()

  # Verifica se preditores possuem classe tabular suportada
  if(!(is.data.frame(sby_predictor_data) || is.matrix(sby_predictor_data))){

    # Aborta quando os preditores nao sao data frame nem matriz
    sby_adanear_abort(
      sby_message = "'sby_predictor_data' deve ser data.frame ou matrix"
    )
  }

  # Verifica se os preditores possuem ao menos uma linha
  if(collapse::fnrow(sby_predictor_data) == 0L){

    # Aborta quando nao ha observacoes para balanceamento
    sby_adanear_abort(
      sby_message = "'sby_predictor_data' deve conter ao menos uma linha"
    )
  }

  # Verifica se alvo e preditores possuem o mesmo numero de linhas
  if(length(sby_target_vector) != collapse::fnrow(sby_predictor_data)){

    # Aborta quando o alvo nao esta alinhado aos preditores
    sby_adanear_abort(
      sby_message = "'sby_target_vector' deve ter o mesmo numero de linhas de 'sby_predictor_data'"
    )
  }

  # Verifica ausencia de valores faltantes nos preditores
  if(anyNA(sby_predictor_data)){

    # Aborta quando preditores contem valores ausentes
    sby_adanear_abort(
      sby_message = "'sby_predictor_data' nao pode conter NA ou NaN"
    )
  }

  # Verifica ausencia de valores faltantes no alvo
  if(anyNA(sby_target_vector)){

    # Aborta quando alvo contem valores ausentes
    sby_adanear_abort(
      sby_message = "'sby_target_vector' nao pode conter NA"
    )
  }

  # Define objeto auxiliar para verificacao de tipos numericos
  if(is.matrix(sby_predictor_data)){

    # Reutiliza matriz de entrada para validacao de tipo
    sby_x_check <- sby_predictor_data
  }else{

    # Reutiliza data frame de entrada para validacao por coluna
    sby_x_check <- sby_predictor_data
  }

  # Avalia se todos os preditores sao numericos
  sby_is_numeric_column <- if(is.matrix(sby_x_check)){

    # Verifica tipo numerico de matriz completa
    is.numeric(sby_x_check)
  }else{

    # Verifica tipo numerico para cada coluna do data frame
    vapply(
      X = sby_x_check,
      FUN = is.numeric,
      FUN.VALUE = logical(1L)
    )
  }

  # Verifica se todos os preditores passaram na validacao numerica
  if(!all(sby_is_numeric_column)){

    # Aborta quando ha preditor nao numerico
    sby_adanear_abort(
      sby_message = "Todos os preditores devem ser numericos"
    )
  }

  # Verifica ausencia de valores infinitos nos preditores numericos
  if(any(!is.finite(as.matrix(sby_x_check)))){

    # Aborta quando preditores contem Inf ou -Inf, que invalidam escala e KNN
    sby_adanear_abort(
      sby_message = "'sby_predictor_data' nao pode conter Inf ou -Inf"
    )
  }

  # Converte alvo para fator para validacao de classes
  sby_target_factor <- as.factor(
    x = sby_target_vector
  )

  # Verifica se o alvo possui exatamente duas classes
  if(nlevels(sby_target_factor) != 2L){

    # Aborta quando o alvo nao e binario
    sby_adanear_abort(
      sby_message = "'sby_target_vector' deve ser binario"
    )
  }

  # Calcula distribuicao de classes do alvo
  sby_class_counts <- table(
    sby_target_factor
  )

  # Verifica se cada classe possui observacoes suficientes
  if(any(sby_class_counts < 2L)){

    # Aborta quando alguma classe tem menos de duas observacoes
    sby_adanear_abort(
      sby_message = "Cada classe deve ter ao menos 2 observacoes"
    )
  }

  # Valida a semente e evita truncamentos silenciosos em set.seed()
  sby_validate_seed(
    sby_seed = sby_seed
  )

  # Retorna sucesso invisivel apos validacao das entradas
  return(invisible(TRUE))
}
####
## Fim
#
