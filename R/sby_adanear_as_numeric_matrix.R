#' Converter preditores para matriz numerica
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_predictor_data Dados preditores em data frame ou matriz
#'
#' @return Matriz numerica com nomes de colunas preservados
#'
#' @noRd
sby_adanear_as_numeric_matrix <- function(sby_predictor_data){
  # Bloqueia Matrix esparsa antes de qualquer conversao densa implicita. As
  # etapas atuais de z-score, KNN e geracao nativa operam sobre matrix double
  # densa; aceitar sparse aqui poderia materializar bases grandes na memoria e
  # parecer travamento do processamento.
  if(sby_is_sparse_matrix(
    sby_x = sby_predictor_data
  )){

    # Aborta com diagnostico explicito em vez de chamar as.matrix silenciosamente
    sby_adanear_abort(
      sby_message = paste0(
        "Matrizes esparsas do pacote Matrix ainda nao sao suportadas por ",
        "sby_adanear/sby_adasyn/sby_nearmiss. Converta conscientemente para ",
        "matrix densa somente se houver memoria suficiente, ou use preditores ",
        "densos ja materializados."
      )
    )
  }

  # Converte a entrada para matriz quando os preditores estao em data frame
  if(is.data.frame(sby_predictor_data)){

    sby_numeric_column <- vapply(
      X = sby_predictor_data,
      FUN = is.numeric,
      FUN.VALUE = logical(1L)
    )
    if(all(sby_numeric_column) && exists("qM", envir = asNamespace("collapse"), mode = "function")){
      sby_x_matrix <- collapse::qM(
        sby_predictor_data
      )
    }else{
      sby_x_matrix <- data.matrix(
        frame = sby_predictor_data
      )
    }
  }else{

    # Reutiliza matriz de entrada para conversao de armazenamento
    sby_x_matrix <- as.matrix(
      x = sby_predictor_data
    )
  }

  # Garante armazenamento double para calculos numericos posteriores
  storage.mode(sby_x_matrix) <- "double"

  # Preserva nomes de colunas existentes ou gerados para a matriz numerica
  colnames(sby_x_matrix) <- sby_adanear_get_column_names(
    sby_predictor_data = sby_predictor_data
  )

  # Retorna matriz numerica normalizada
  return(sby_x_matrix)
}
####
## Fim
#
