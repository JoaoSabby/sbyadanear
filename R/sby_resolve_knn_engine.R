#' Resolver engine KNN automatico quando solicitado
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado.
#' A heuristica auto considera dimensionalidade, tamanho da base e metrica para escolher entre FNN (exato) e RcppHNSW (aproximado).
#'
#' @param sby_knn_engine Engine KNN informado pelo chamador
#' @param sby_knn_workers Numero de workers validado para consulta KNN
#' @param sby_knn_distance_metric Metrica de distancia configurada
#' @param sby_row_count Numero total de linhas dos dados (opcional)
#' @param sby_predictor_column_count Numero de colunas preditoras (opcional)
#'
#' @return Nome do engine KNN resolvido
#' @noRd
sby_resolve_knn_engine <- function(
  sby_knn_engine,
  sby_knn_workers,
  sby_knn_distance_metric = "euclidean",
  sby_row_count = NA_integer_,
  sby_predictor_column_count = NA_integer_
){

  # Retorna engine explicito quando modo automatico nao foi solicitado
  if(!identical(
    x = sby_knn_engine,
    y = "auto"
  )){

    # Mantem a escolha explicita do chamador
    return(sby_knn_engine)
  }

  # Heuristica 1: metricas nao euclidianas exigem RcppHNSW (FNN nao suporta)
  if(!identical(sby_knn_distance_metric, "euclidean")){
    sby_adanear_inform(
      sby_message = paste0(
        "KNN automatico: sby_knn_engine = \"RcppHNSW\". ",
        "Justificativa: a metrica \"", sby_knn_distance_metric,
        "\" so e suportada pelo engine RcppHNSW; FNN aceita apenas euclidean."
      )
    )
    return("RcppHNSW")
  }

  # Heuristica 2: bases grandes e de alta dimensionalidade tendem a se beneficiar
  # de HNSW (aproximado) em vez de FNN exato. Os limiares sao conservadores e podem
  # ser sobrescritos por opcao instenginer.sby_auto_engine_hnsw_min_cells.
  sby_min_cells <- as.numeric(getOption(
    x = "instenginer.sby_auto_engine_hnsw_min_cells",
    default = 5e6
  ))
  sby_have_dims <- !is.na(sby_row_count) && !is.na(sby_predictor_column_count) &&
    is.finite(sby_row_count) && is.finite(sby_predictor_column_count)

  if(sby_have_dims){
    sby_cells <- as.numeric(sby_row_count) * as.numeric(sby_predictor_column_count)
    if(sby_cells >= sby_min_cells && sby_predictor_column_count >= 50L){
      sby_adanear_inform(
        sby_message = paste0(
          "KNN automatico: sby_knn_engine = \"RcppHNSW\". ",
          "Justificativa: dados grandes (n=", sby_row_count,
          ", p=", sby_predictor_column_count, ") favorecem busca aproximada HNSW. ",
          "Para forcar FNN exato, defina sby_knn_engine = \"FNN\"."
        )
      )
      return("RcppHNSW")
    }
  }

  # Heuristica 3: caso padrao - FNN exato, sequencial ou paralelo por blocos
  sby_adanear_inform(
    sby_message = paste0(
      "KNN automatico: sby_knn_engine = \"FNN\". ",
      "Justificativa: FNN e o engine exato padrao; quando sby_knn_workers e ",
      "maior que 1, o pacote paraleliza a consulta por blocos sem depender ",
      "de dependencias adicionais."
    )
  )
  return("FNN")
}
####
## Fim
#
