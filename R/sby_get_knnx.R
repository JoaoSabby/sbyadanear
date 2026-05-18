#' Executar consulta KNN usando FNN ou RcppHNSW
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_data Matriz de referencia para busca KNN
#' @param sby_query Matriz de consulta para busca KNN
#' @param sby_k Numero de vizinhos solicitados
#' @param sby_knn_algorithm Algoritmo KNN configurado
#' @param sby_knn_engine Engine KNN configurado
#' @param sby_knn_workers Numero de workers KNN configurado
#' @param sby_knn_hnsw_m Conectividade HNSW configurada
#' @param sby_knn_hnsw_ef Lista dinamica HNSW configurada
#' @param sby_knn_query_chunk_size Tamanho de bloco para consultas KNN
#' @param sby_query_is_data Indicador de que consulta e referencia sao a mesma matriz
#'
#' @return Lista com matrizes `nn.index` e `nn.dist`
#' @noRd
sby_get_knnx <- function(
  sby_data,
  sby_query,
  sby_k,
  sby_knn_algorithm,
  sby_knn_engine,
  sby_knn_distance_metric,
  sby_knn_workers,
  sby_knn_hnsw_m,
  sby_knn_hnsw_ef,
  sby_knn_query_chunk_size = getOption("instenginer.sby_knn_query_chunk_size", 1000L),
  sby_query_is_data = FALSE,
  sby_knn_return = c("both", "index", "dist")
){
  
  # Verifica se ha solicitacao de interrupcao antes da consulta KNN
  sby_adanear_check_user_interrupt()

  # Resolve quais componentes KNN devem permanecer no objeto retornado
  sby_knn_return <- match.arg(sby_knn_return)

  # Valida tamanho de bloco para consultas KNN interrompiveis
  sby_knn_query_chunk_size <- sby_validate_knn_query_chunk_size(
    sby_knn_query_chunk_size = sby_knn_query_chunk_size
  )

  # Aplica normalizacao L2 obrigatoria para metricas angulares ou produto interno
  if(!identical(
    x = sby_knn_distance_metric,
    y = "euclidean"
  )){

    # Normaliza referencia depois do z-score e antes do engine KNN
    sby_data <- sby_normalize_l2(
      sby_x_matrix = sby_data
    )

    # Reutiliza a mesma matriz normalizada em consultas self-KNN para evitar
    # uma segunda passagem completa por dados grandes de minoria.
    sby_query <- if(isTRUE(sby_query_is_data)){
      sby_data
    }else{
      sby_normalize_l2(
        sby_x_matrix = sby_query
      )
    }
  }

  # Executa consulta pelo engine FNN quando selecionado
  if(identical(
    x = sby_knn_engine,
    y = "FNN"
  )){

    # Bloqueia metricas nao euclidianas porque FNN implementa apenas distancia euclidiana neste pacote
    if(!identical(
      x = sby_knn_distance_metric,
      y = "euclidean"
    )){

      # Aborta combinacao incompatvel para evitar regressao silenciosa de metrica
      sby_adanear_abort(
        sby_message = "'sby_knn_engine = FNN' suporta apenas 'sby_knn_distance_metric = euclidean'"
      )
    }

    # Verifica disponibilidade do pacote FNN
    if(!requireNamespace(
      package = "FNN",
      quietly = TRUE
    )){

      # Aborta quando o engine FNN nao esta instalado
      sby_adanear_abort(
        sby_message = "'sby_knn_engine = FNN' requer o pacote FNN"
      )
    }

    # Bloqueia algoritmos exclusivos dos outros engines no FNN
    if(!(sby_knn_algorithm %in% c("auto", "kd_tree", "cover_tree", "brute"))){

      # Aborta algoritmo incompatvel com FNN
      sby_adanear_abort(
        sby_message = "'sby_knn_algorithm' deve ser um de 'auto', 'kd_tree', 'cover_tree' ou 'brute' quando 'sby_knn_engine = FNN'"
      )
    }

    # Consulta vizinhos FNN em blocos interrompiveis, sequenciais ou paralelos
    sby_knn_result <- sby_query_knn_in_chunks(
      sby_query = sby_query,
      sby_k = sby_k,
      sby_knn_query_chunk_size = sby_knn_query_chunk_size,
      sby_knn_workers = sby_knn_workers,
      sby_query_fun = function(sby_query_chunk){
        # Executa consulta FNN para o bloco corrente
        return(FNN::get.knnx(
          data = sby_data,
          query = sby_query_chunk,
          k = sby_k,
          algorithm = sby_knn_algorithm
        ))
      }
    )

    # Verifica se ha solicitacao de interrupcao apos consulta FNN
    sby_adanear_check_user_interrupt()

    # Retorna resultado KNN produzido pelo engine FNN
    return(sby_trim_knn_result(sby_knn_result, sby_knn_return))
  }

  # Executa consulta pelo engine RcppHNSW quando selecionado
  if(identical(
    x = sby_knn_engine,
    y = "RcppHNSW"
  )){

    # Verifica disponibilidade do pacote RcppHNSW
    if(!requireNamespace(
      package = "RcppHNSW",
      quietly = TRUE
    )){

      # Aborta quando o engine RcppHNSW nao esta instalado
      sby_adanear_abort(
        sby_message = "'sby_knn_engine = RcppHNSW' requer o pacote RcppHNSW. Instale-o com install.packages('RcppHNSW')."
      )
    }

    # Define parametro efetivo de busca HNSW limitado pelo tamanho dos dados
    sby_effective_ef <- min(
      max(
        as.integer(sby_knn_hnsw_ef),
        as.integer(sby_k)
      ),
      nrow(sby_data)
    )

    # Define rotina HNSW completa em funcao local. O caminho padrao executa
    # diretamente para evitar custo alto de fork/serializacao em bases grandes;
    # o fork permanece disponivel por opcao para investigacoes interativas nas
    # quais a capacidade de matar chamadas nativas bloqueantes seja prioritaria.
    sby_run_hnsw_query <- function(){
      # Constroi indice HNSW para a matriz de referencia
      sby_hnsw_index <- RcppHNSW::hnsw_build(
        X = sby_data,
        distance = sby_knn_distance_metric,
        M = as.integer(sby_knn_hnsw_m),
        ef = sby_effective_ef,
        verbose = FALSE,
        progress = "none",
        n_threads = sby_knn_workers,
        byrow = TRUE
      )

      # Verifica se ha solicitacao de interrupcao apos construcao do indice HNSW
      sby_adanear_check_user_interrupt()

      # Valida tamanho de bloco especifico para consultas HNSW. O valor
      # padrao acompanha o tamanho geral de blocos KNN para evitar milhares
      # de chamadas nativas pequenas, que podem transformar execucoes de
      # minutos em horas em bases grandes.
      sby_hnsw_query_chunk_size <- sby_validate_knn_query_chunk_size(
        sby_knn_query_chunk_size = getOption(
          "instenginer.sby_hnsw_query_chunk_size",
          sby_knn_query_chunk_size
        )
      )

      # Consulta vizinhos HNSW em blocos interrompiveis
      sby_knn_result <- sby_query_knn_in_chunks(
        sby_query = sby_query,
        sby_k = sby_k,
        sby_knn_query_chunk_size = sby_hnsw_query_chunk_size,
        sby_knn_workers = 1L,
        sby_query_fun = function(sby_query_chunk){
          # Executa busca HNSW para o bloco corrente
          sby_hnsw_result <- RcppHNSW::hnsw_search(
            X = sby_query_chunk,
            ann = sby_hnsw_index,
            k = sby_k,
            ef = sby_effective_ef,
            verbose = FALSE,
            progress = "none",
            n_threads = sby_knn_workers,
            byrow = TRUE
          )

          # Retorna indices e distancias no contrato comum de KNN
          return(list(
            nn.index = sby_hnsw_result$idx,
            nn.dist  = sby_hnsw_result$dist
          ))
        }
      )

      # Verifica se ha solicitacao de interrupcao apos consulta HNSW
      sby_adanear_check_user_interrupt()

      # Retorna resultado KNN produzido pelo engine HNSW
      return(sby_trim_knn_result(sby_knn_result, sby_knn_return))
    }

    # Executa HNSW em fork somente quando solicitado explicitamente. Em bases
    # grandes, retornar matrizes de vizinhos pelo pipe do processo filho pode
    # dominar o tempo de execucao; por isso o padrao privilegia desempenho e
    # mantem a interrupcao cooperativa entre blocos de consulta.
    if(isTRUE(getOption(
      x = "instenginer.sby_hnsw_interruptible_fork",
      default = FALSE
    ))){
      return(sby_run_interruptible_fork(
        sby_run_hnsw_query()
      ))
    }

    # Mantem caminho direto opcional para ambientes que desativarem isolamento por fork
    return(sby_run_hnsw_query())
  }

  # Aborta engines desconhecidos para manter o contrato publico atual
  sby_adanear_abort(
    sby_message = "'sby_knn_engine' deve ser um de 'auto', 'FNN' ou 'RcppHNSW'. Use FNN com sby_knn_workers > 1L para paralelismo exato."
  )
}

#' Executar consultas KNN em blocos interrompiveis
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_query Matriz de consulta para busca KNN
#' @param sby_k Numero de vizinhos solicitados
#' @param sby_knn_query_chunk_size Tamanho de bloco para consultas KNN
#' @param sby_knn_workers Numero de workers usados para processar blocos de consulta
#' @param sby_query_fun Funcao que executa consulta KNN em um bloco
#'
#' @return Lista com matrizes `nn.index` e `nn.dist`
#' @noRd
sby_query_knn_in_chunks <- function(sby_query, sby_k, sby_knn_query_chunk_size, sby_knn_workers, sby_query_fun){
  # Calcula numero de linhas da matriz de consulta
  sby_query_rows <- nrow(sby_query)

  # Executa consulta diretamente quando os dados cabem em um unico bloco
  if(sby_query_rows <= sby_knn_query_chunk_size){

    # Retorna resultado da funcao de consulta sem particionamento
    return(sby_query_fun(
      sby_query
    ))
  }

  # Cria descritores de blocos preservando a ordem original das linhas
  sby_chunk_starts <- seq.int(
    from = 1L,
    to = sby_query_rows,
    by = sby_knn_query_chunk_size
  )
  sby_chunk_ranges <- lapply(
    X = sby_chunk_starts,
    FUN = function(sby_chunk_start){
      # Define intervalo de linhas do bloco corrente
      sby_chunk_end <- min(
        sby_chunk_start + sby_knn_query_chunk_size - 1L,
        sby_query_rows
      )

      # Retorna indices globais do bloco
      return(seq.int(
        from = sby_chunk_start,
        to = sby_chunk_end
      ))
    }
  )

  # Funcao local que executa um bloco e preserva seus indices globais
  sby_run_chunk <- function(sby_chunk_index){
    # Executa consulta KNN para o bloco corrente
    sby_chunk_result <- sby_query_fun(
      sby_query[sby_chunk_index, , drop = FALSE]
    )

    # Retorna resultado associado ao intervalo original para consolidacao ordenada
    return(list(
      sby_chunk_index = sby_chunk_index,
      nn.index = sby_chunk_result$nn.index,
      nn.dist = sby_chunk_result$nn.dist
    ))
  }

  # Limita workers efetivos ao numero real de blocos disponiveis
  sby_effective_workers <- min(
    as.integer(sby_knn_workers),
    length(sby_chunk_ranges)
  )

  # Processa consultas KNN em blocos sequenciais quando apenas um worker e usado
  if(sby_effective_workers <= 1L){
    sby_chunk_results <- vector(
      mode = "list",
      length = length(sby_chunk_ranges)
    )

    # Processa cada bloco com verificacoes cooperativas de interrupcao
    for(sby_chunk_position in seq_along(sby_chunk_ranges)){
      # Verifica interrupcao antes de cada bloco de consulta
      sby_adanear_check_user_interrupt()

      # Executa consulta KNN para o bloco corrente
      sby_chunk_results[[sby_chunk_position]] <- sby_run_chunk(
        sby_chunk_ranges[[sby_chunk_position]]
      )

      # Verifica interrupcao apos cada bloco de consulta
      sby_adanear_check_user_interrupt()
    }
  }else if(identical(
    x = .Platform$OS.type,
    y = "windows"
  )){
    # Em Windows usa cluster PSOCK; os blocos sao independentes e a ordem do
    # parLapply acompanha a ordem de entrada, preservando exatidao e alinhamento.
    sby_cluster <- parallel::makeCluster(
      spec = sby_effective_workers,
      type = "PSOCK"
    )
    on.exit(
      expr = parallel::stopCluster(sby_cluster),
      add = TRUE
    )

    # Executa blocos em paralelo por sockets
    sby_chunk_results <- parallel::parLapply(
      cl = sby_cluster,
      X = sby_chunk_ranges,
      fun = sby_run_chunk
    )
  }else{
    # Em sistemas POSIX usa fork para reduzir copias fisicas da matriz de referencia
    # e manter a consolidacao ordenada dos resultados independentes por bloco.
    sby_chunk_results <- parallel::mclapply(
      X = sby_chunk_ranges,
      FUN = sby_run_chunk,
      mc.cores = sby_effective_workers,
      mc.preschedule = TRUE
    )
  }

  # Inicializa matrizes de saida
  sby_nn_index <- matrix(
    data = NA_integer_,
    nrow = sby_query_rows,
    ncol = sby_k
  )
  sby_nn_dist <- matrix(
    data = NA_real_,
    nrow = sby_query_rows,
    ncol = sby_k
  )

  # Consolida resultados por indice global, sem depender da ordem de conclusao
  for(sby_chunk_result in sby_chunk_results){
    # Copia indices e distancias do bloco para as matrizes completas
    sby_nn_index[sby_chunk_result$sby_chunk_index, ] <- sby_chunk_result$nn.index
    sby_nn_dist[sby_chunk_result$sby_chunk_index, ]  <- sby_chunk_result$nn.dist
  }

  # Verifica se ha solicitacao de interrupcao apos consolidacao dos blocos
  sby_adanear_check_user_interrupt()

  # Retorna matrizes completas de indices e distancias
  return(list(
    nn.index = sby_nn_index,
    nn.dist  = sby_nn_dist
  ))
}
####
## Fim
#


#' Descartar componentes KNN nao solicitados
#'
#' @param sby_knn_result Resultado KNN com `nn.index` e/ou `nn.dist`
#' @param sby_knn_return Componentes a manter
#'
#' @return Lista KNN reduzida
#' @noRd
sby_trim_knn_result <- function(sby_knn_result, sby_knn_return){
  if(identical(sby_knn_return, "index")){
    sby_knn_result$nn.dist <- NULL
  }else if(identical(sby_knn_return, "dist")){
    sby_knn_result$nn.index <- NULL
  }
  return(sby_knn_result)
}
