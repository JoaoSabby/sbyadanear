#' Preparar etapa de balanceamento ADASYN e NearMiss
#'
#' Detecta a configuração inerte antes de calcular estatísticas de treinamento. Quando
#' ambas as proporções são zero, devolve a etapa treinada e marcada para ser ignorada
#' por `bake()`, evitando qualquer trabalho de amostragem.
#'
#' @param x Objeto não treinado da etapa `sby_step_adanear`.
#' @param training Dados de treinamento da receita.
#' @param info Metadados da receita, ou `NULL`.
#' @param ... Argumentos adicionais preservados para compatibilidade com o método S3.
#'
#' @return Objeto treinado da etapa `sby_step_adanear`.
#'
#' @export
prep.step_sby_step_adanear <- function(x, training, info = NULL, ...){
  # Extrai a proporcao de sobreamostragem armazenada na etapa.
  sby_adasyn_ratio <- x$sby_adasyn_ratio

  # Extrai a proporcao de subamostragem armazenada na etapa.
  sby_nearmiss_ratio <- x$sby_nearmiss_ratio

  # Identifica a inercia somente quando ambas as proporcoes sao estritamente zero.
  ratios_are_zero <- identical(as.numeric(sby_adasyn_ratio), 0) &&
    identical(as.numeric(sby_nearmiss_ratio), 0)

  # Interrompe a preparacao estatistica quando nenhuma amostragem foi solicitada.
  if(ratios_are_zero){
    # Registra o estado treinado exigido pelo contrato nativo de recipes.
    x$trained <- TRUE

    # Sincroniza o estado treinado mantido internamente pelo pacote.
    x$sby_trained <- TRUE

    # Sinaliza a recipes que a etapa deve ser ignorada em aplicacoes futuras.
    x$skip <- TRUE

    # Sincroniza a sinalizacao de bypass mantida internamente pelo pacote.
    x$sby_skip <- TRUE

    # Retorna a etapa inerte sem avaliar seletores ou estatisticas de treinamento.
    return(x)
  }

  # Executa a preparacao completa quando ao menos uma proporcao e positiva.
  return(sby_prep_step_sampling(
    x = x,
    training = training,
    info = info,
    sby_step_name = "sby_step_adanear()"
  ))
}
####
## Fim
#
