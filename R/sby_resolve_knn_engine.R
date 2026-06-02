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

  sby_min_cells <- 5e6
  sby_have_dims <- !is.na(sby_row_count) && !is.na(sby_predictor_column_count) &&
    is.finite(sby_row_count) && is.finite(sby_predictor_column_count)
  sby_cells <- 0
  if(sby_have_dims){
    sby_cells <- as.numeric(sby_row_count) * as.numeric(sby_predictor_column_count)
  }

  sby_route_code <- kit::nif(
    !identical(sby_knn_distance_metric, "euclidean"),
    1L,
    sby_have_dims && sby_cells >= sby_min_cells && sby_predictor_column_count >= 50L,
    2L,
    default = 3L
  )
  sby_route_code <- as.integer(sby_route_code[[1L]])

  if(identical(sby_route_code, 1L)){
    sby_adanear_inform(
      sby_message = paste0(
        "KNN automático: sby_knn_engine = \"RcppHNSW\". ",
        "Justificativa: a métrica \"", sby_knn_distance_metric,
        "\" só é suportada pelo engine RcppHNSW; FNN aceita apenas euclidean."
      )
    )
    return("RcppHNSW")
  }

  if(identical(sby_route_code, 2L)){
    sby_adanear_inform(
      sby_message = paste0(
        "KNN automático: sby_knn_engine = \"RcppHNSW\". ",
        "Justificativa: dados grandes (n=", sby_row_count,
        ", p=", sby_predictor_column_count, ") favorecem busca aproximada HNSW. ",
        "Para forçar FNN exato, defina sby_knn_engine = \"FNN\"."
      )
    )
    return("RcppHNSW")
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
