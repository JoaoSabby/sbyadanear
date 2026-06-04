#' Resolver algoritmo KNN automatico quando solicitado
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado.
#' A resolucao considera o engine efetivo para que `auto` produza algoritmos compativeis com FNN, RcppHNSW, KernelKnn ou bigKNN.
#'
#' @param sby_knn_algorithm Algoritmo KNN informado pelo chamador
#' @param sby_predictor_column_count Quantidade de colunas preditoras
#' @param sby_knn_engine Engine KNN resolvido
#'
#' @return Nome do algoritmo KNN resolvido
#' @noRd
sby_resolve_knn_algorithm <- function(sby_knn_algorithm, sby_predictor_column_count, sby_knn_engine){

  # Retorna algoritmo explicito quando modo automatico nao foi solicitado
  if(!identical(
    x = sby_knn_algorithm,
    y = "auto"
  )){

    # Mantem a escolha explicita do chamador
    return(sby_knn_algorithm)
  }

  # Engines externos gerenciam seu algoritmo internamente
  if(sby_knn_engine %in% c("native", "RcppHNSW", "KernelKnn", "bigKNN")){

    # Informa ao usuario que o engine selecionado gerencia o algoritmo internamente
    sby_adanear_inform(
      sby_message = paste0(
        "KNN automatico: sby_knn_algorithm = \"auto\" com sby_knn_engine = \"",
        sby_knn_engine, "\". Justificativa: o engine selecionado controla ",
        "internamente a estrategia de busca, entao algoritmos externos nao sao aplicados nessa rota."
      )
    )

    # Mantem marcador automatico porque o engine ignora algoritmos externos
    return("auto")
  }

  # Seleciona algoritmo FNN por dimensionalidade quando o modo automatico e usado
  if(identical(
    x = sby_knn_engine,
    y = "FNN"
  )){

    # Usa busca bruta para dimensionalidade mais alta
    if(sby_predictor_column_count > 15L){

      # Informa ao usuario a selecao automatica e a justificativa da decisao
      sby_adanear_inform(
        sby_message = paste0(
          "KNN automatico: sby_knn_algorithm = \"brute\" com sby_knn_engine = \"FNN\". ",
          "Justificativa: os dados tem ", sby_predictor_column_count,
          " colunas preditoras, acima do limite de 15; nessa dimensionalidade, ",
          "a busca bruta tende a ser mais estavel que estruturas de arvore."
        )
      )

      return("brute")
    }

    # Informa ao usuario a selecao automatica e a justificativa da decisao
    sby_adanear_inform(
      sby_message = paste0(
        "KNN automatico: sby_knn_algorithm = \"kd_tree\" com sby_knn_engine = \"FNN\". ",
        "Justificativa: os dados tem ", sby_predictor_column_count,
        " colunas preditoras, ate o limite de 15; nessa dimensionalidade, ",
        "kd_tree costuma acelerar consultas exatas de vizinhos."
      )
    )

    # Usa kd-tree para dimensionalidade mais baixa
    return("kd_tree")
  }

  # Aborta engines desconhecidos para evitar combinacoes obsoletas ou removidas
  sby_adanear_abort(
    sby_message = "'sby_knn_engine' deve ser um de 'auto', 'native', 'FNN', 'RcppHNSW', 'KernelKnn' ou 'bigKNN'"
  )

}
####
## Fim
#
