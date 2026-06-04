#' Resolver engine KNN automatico quando solicitado
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado.
#' A heuristica auto e conservadora: para ADASYN e NearMiss, busca aproximada so e selecionada quando explicitamente permitida.
#'
#' @param sby_knn_engine Engine KNN informado pelo chamador
#' @param sby_knn_workers Numero de workers validado para consulta KNN
#' @param sby_knn_distance_metric Metrica de distancia configurada
#' @param sby_row_count Numero total de linhas dos dados (opcional)
#' @param sby_predictor_column_count Numero de colunas preditoras (opcional)
#' @param sby_knn_allow_approx Permite que a heuristica automatica selecione engines aproximados
#'
#' @return Nome do engine KNN resolvido
#' @noRd
sby_resolve_knn_engine <- function(
  sby_knn_engine,
  sby_knn_workers,
  sby_knn_distance_metric = "euclidean",
  sby_row_count = NA_integer_,
  sby_predictor_column_count = NA_integer_,
  sby_knn_allow_approx = getOption("sbyadanear.sby_knn_allow_approx", FALSE)
){

  # Retorna engine explicito quando modo automatico nao foi solicitado
  if(!identical(
    x = sby_knn_engine,
    y = "auto"
  )){

    # Mantem a escolha explicita do chamador
    return(sby_knn_engine)
  }

  sby_knn_allow_approx <- isTRUE(sby_knn_allow_approx)

  # Metricas nao euclidianas exigem HNSW na integracao atual. Como ADASYN e
  # NearMiss sao sensiveis a vizinhos, a rota aproximada so e escolhida quando
  # o usuario permite explicitamente.
  if(!identical(sby_knn_distance_metric, "euclidean")){
    if(isTRUE(sby_knn_allow_approx) && requireNamespace("RcppHNSW", quietly = TRUE)){
      sby_adanear_inform(
        sby_message = paste0(
          "KNN automatico: sby_knn_engine = \"RcppHNSW\". ",
          "Justificativa: a metrica \"", sby_knn_distance_metric,
          "\" nao e suportada pelas rotas exatas atuais e busca aproximada foi permitida."
        )
      )
      return("RcppHNSW")
    }

    sby_adanear_abort(
      sby_message = paste0(
        "KNN automatico exato suporta apenas 'sby_knn_distance_metric = euclidean'. ",
        "Para usar '", sby_knn_distance_metric,
        "' automaticamente, defina a opcao sbyadanear.sby_knn_allow_approx = TRUE ou escolha sby_knn_engine = 'RcppHNSW' explicitamente."
      )
    )
  }

  # Preferencia conservadora: usa a engine nativa exata quando disponivel.
  if(sby_adanear_native_available()){
    sby_adanear_inform(
      sby_message = paste0(
        "KNN automatico: sby_knn_engine = \"native\". ",
        "Justificativa: a engine nativa executa KNN euclidiano exato denso ",
        "com contrato nn.index/nn.dist e evita selecionar busca aproximada sem permissao explicita."
      )
    )
    return("native")
  }

  # Fallback exato externo quando o kernel nativo nao estiver carregado.
  if(requireNamespace("FNN", quietly = TRUE)){
    sby_adanear_inform(
      sby_message = paste0(
        "KNN automatico: sby_knn_engine = \"FNN\". ",
        "Justificativa: a engine nativa nao esta disponivel e FNN fornece busca euclidiana exata."
      )
    )
    return("FNN")
  }

  # Ultimo recurso aproximado, apenas com permissao explicita.
  if(isTRUE(sby_knn_allow_approx) && requireNamespace("RcppHNSW", quietly = TRUE)){
    sby_adanear_inform(
      sby_message = paste0(
        "KNN automatico: sby_knn_engine = \"RcppHNSW\". ",
        "Justificativa: nenhuma engine exata esta disponivel e busca aproximada foi permitida."
      )
    )
    return("RcppHNSW")
  }

  sby_adanear_abort(
    sby_message = "Nenhuma engine KNN exata compativel esta disponivel. Instale FNN ou carregue as rotinas nativas do pacote."
  )
}
####
## Fim
#
