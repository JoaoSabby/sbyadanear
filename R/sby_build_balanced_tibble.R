#' Construir tibble balanceado com alvo padronizado
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_predictor_data Dados preditores ja balanceados
#'
#' @param sby_target_vector Vetor alvo correspondente aos preditores
#'
#' @return Tibble com preditores e coluna TARGET
#'
#' @noRd
sby_build_balanced_tibble <- function(sby_predictor_data, sby_target_vector){
  
  # Verifica se nome reservado de alvo ja existe nos preditores
  if("TARGET" %in% names(sby_predictor_data)){

    # Aborta para evitar sobrescrita ambigua da coluna de alvo. A coluna
    # 'TARGET' e reservada pelo pacote para o alvo no tibble de saida; se
    # 'TARGET' aparece como preditor, o usuario deve renomea-la antes da
    # chamada para evitar colisao silenciosa.
    sby_adanear_abort(
      sby_message = paste0(
        "A coluna 'TARGET' e reservada pelo pacote para o alvo no tibble ",
        "retornado. Renomeie a coluna preditora chamada 'TARGET' antes de ",
        "chamar as funcoes de balanceamento do sbyadanear."
      )
    )
  }

  # Converte preditores para data frame antes da composicao tabular
  sby_out <- as.data.frame(
    x = sby_predictor_data,
    stringsAsFactors = FALSE
  )

  # Anexa coluna alvo com fator preservado
  sby_out$TARGET <- sby_target_vector

  # Retorna tibble balanceado para consumo externo
  return(tibble::as_tibble(
    x = sby_out
  ))
}
####
## Fim
#
