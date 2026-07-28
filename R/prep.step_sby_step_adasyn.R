#' Preparar etapa de balanceamento ADASYN
#'
#' Detecta a configuração inerte antes de calcular estatísticas de treinamento. Quando
#' a proporção é zero, devolve a etapa treinada e marcada para ser ignorada por
#' `bake()`, evitando qualquer trabalho de amostragem.
#'
#' @param x Objeto não treinado da etapa `sby_step_adasyn`.
#' @param training Dados de treinamento da receita.
#' @param info Metadados da receita, ou `NULL`.
#' @param ... Argumentos adicionais preservados para compatibilidade com o método S3.
#'
#' @return Objeto treinado da etapa `sby_step_adasyn`.
#'
#' @export
prep.step_sby_step_adasyn <- function(x, training, info = NULL, ...){
  # Extrai a proporcao de sobreamostragem armazenada na etapa.
  sby_adasyn_ratio <- x$sby_adasyn_ratio

  # Identifica a inercia somente quando a proporcao e estritamente zero.
  ratio_is_zero <- identical(as.numeric(sby_adasyn_ratio), 0)

  # Interrompe a preparacao estatistica quando nenhuma amostragem foi solicitada.
  if(ratio_is_zero){
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

  # Executa a preparacao completa quando a proporcao e positiva.
  return(sby_prep_step_sampling(
    x = x,
    training = training,
    info = info,
    sby_step_name = "sby_step_adasyn()"
  ))
}
####
## Fim
#
