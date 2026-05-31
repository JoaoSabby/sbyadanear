#' Gerar amostras sinteticas ADASYN em matriz ja escalada
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_x_scaled Matriz de preditores ja padronizada
#' @param sby_target_factor Fator binario associado aos preditores
#' @param sby_synthetic_count Numero de amostras sinteticas a gerar
#' @param sby_knn_over_k Numero de vizinhos usado pelo ADASYN
#' @param sby_knn_algorithm Algoritmo KNN configurado
#' @param sby_knn_engine Engine KNN configurado
#' @param sby_knn_distance_metric Metrica de distancia KNN configurada
#' @param sby_knn_workers Numero de workers KNN configurado
#' @param sby_knn_hnsw_m Conectividade HNSW configurada
#' @param sby_knn_hnsw_ef Lista dinamica HNSW configurada
#'
#' @return Lista com matriz expandida `x` e fator expandido `y`
#' @noRd
sby_generate_adasyn_samples <- function(
  sby_x_scaled,
  sby_target_factor,
  sby_synthetic_count,
  sby_knn_over_k,
  sby_knn_algorithm,
  sby_knn_engine,
  sby_knn_distance_metric,
  sby_knn_workers,
  sby_knn_hnsw_m,
  sby_knn_hnsw_ef,
  sby_knn_query_chunk_size
){
  
  # Identifica papeis de classe para localizar a minoria
  sby_class_roles <- sby_get_binary_class_roles(
    sby_target_factor = sby_target_factor
  )

  # Extrai indices e matriz da classe minoritaria
  sby_minority_index  <- which(
    x = sby_target_factor == sby_class_roles$sby_minority_label
  )
  sby_minority_matrix <- sby_x_scaled[sby_minority_index, , drop = FALSE]

  # Calcula vizinhos de cada linha minoritaria contra todo o conjunto escalado
  sby_effective_all_k     <- min(
    as.integer(sby_knn_over_k) + 1L,
    nrow(sby_x_scaled)
  )
  sby_all_neighbor_result <- sby_get_knnx(
    sby_data                    = sby_x_scaled,
    sby_query                   = sby_minority_matrix,
    sby_k                       = sby_effective_all_k,
    sby_knn_algorithm           = sby_knn_algorithm,
    sby_knn_engine             = sby_knn_engine,
    sby_knn_distance_metric         = sby_knn_distance_metric,
    sby_knn_workers             = sby_knn_workers,
    sby_knn_hnsw_m                  = sby_knn_hnsw_m,
    sby_knn_hnsw_ef                 = sby_knn_hnsw_ef,
    sby_knn_query_chunk_size      = sby_knn_query_chunk_size,
    sby_knn_return                = "index"
  )

  # Verifica se ha solicitacao de interrupcao apos a primeira consulta KNN
  sby_adanear_check_user_interrupt()

  # Remove o proprio ponto dos vizinhos quando ha vizinhos suficientes
  sby_neighbor_index <- sby_all_neighbor_result$nn.index
  sby_desired_all_k  <- min(
    as.integer(sby_knn_over_k),
    nrow(sby_x_scaled) - 1L
  )

  # Verifica se a consulta retornou mais de um vizinho candidato
  if(sby_effective_all_k > 1L){

    # Remove autorreferencias da matriz de vizinhos globais
    sby_neighbor_index <- sby_drop_self_neighbor_index(
      sby_neighbor_index = sby_neighbor_index,
      sby_self_index     = sby_minority_index,
      sby_desired_k      = sby_desired_all_k
    )
  }

  # Calcula mascara de vizinhos pertencentes a classe majoritaria
  sby_majority_mask <- sby_target_factor[as.vector(sby_neighbor_index)] == sby_class_roles$sby_majority_label

  # Reaproveita o vetor logico como matriz para evitar copia extra em bases grandes
  dim(sby_majority_mask) <- dim(sby_neighbor_index)

  # Verifica se ha solicitacao de interrupcao antes do calculo de razoes
  sby_adanear_check_user_interrupt()

  # Calcula proporcao de vizinhos majoritarios por linha minoritaria
  sby_majority_ratio <- rowMeans(
    x = sby_majority_mask
  )

  # Verifica se ha solicitacao de interrupcao antes da ponderacao de geracao
  sby_adanear_check_user_interrupt()

  # Define pesos de geracao a partir da dificuldade local da minoria
  if(sum(sby_majority_ratio) <= 0){

    # Distribui geracao uniformemente quando nao ha vizinhos majoritarios
    sby_generation_weights <- rep.int(
      x = 1 / length(sby_minority_index),
      times = length(sby_minority_index)
    )
  }else{

    # Normaliza razoes de maioria para probabilidades de geracao
    sby_generation_weights <- sby_majority_ratio / sum(sby_majority_ratio)
  }

  # Converte pesos de geracao em contagens inteiras por linha minoritaria
  sby_raw_counts        <- sby_synthetic_count * sby_generation_weights
  sby_synthetic_per_row <- floor(sby_raw_counts)
  sby_remaining         <- sby_synthetic_count - sum(sby_synthetic_per_row)

  # Distribui amostras restantes pelas maiores partes fracionarias
  if(sby_remaining > 0L){

    # Ordena linhas minoritarias por fracao residual decrescente
    sby_fractional_order <- order(
      sby_raw_counts - sby_synthetic_per_row,
      decreasing = TRUE
    )

    # Acrescenta uma amostra as linhas com maior residuo fracionario
    sby_synthetic_per_row[sby_fractional_order[seq_len(sby_remaining)]] <- sby_synthetic_per_row[sby_fractional_order[seq_len(sby_remaining)]] + 1L
  }

  # Calcula vizinhos minoritarios usados para interpolacao sintetica
  sby_effective_minority_k     <- min(
    as.integer(sby_knn_over_k) + 1L,
    nrow(sby_minority_matrix)
  )
  sby_minority_neighbor_result <- sby_get_knnx(
    sby_data                    = sby_minority_matrix,
    sby_query                   = sby_minority_matrix,
    sby_k                       = sby_effective_minority_k,
    sby_knn_algorithm           = sby_knn_algorithm,
    sby_knn_engine             = sby_knn_engine,
    sby_knn_distance_metric         = sby_knn_distance_metric,
    sby_knn_workers             = sby_knn_workers,
    sby_knn_hnsw_m                  = sby_knn_hnsw_m,
    sby_knn_hnsw_ef                 = sby_knn_hnsw_ef,
    sby_knn_query_chunk_size      = sby_knn_query_chunk_size,
    sby_query_is_data           = TRUE,
    sby_knn_return                = "index"
  )

  # Verifica se ha solicitacao de interrupcao apos a consulta KNN minoritaria
  sby_adanear_check_user_interrupt()

  # Remove autorreferencias dos vizinhos minoritarios quando possivel
  sby_minority_neighbor_index <- sby_minority_neighbor_result$nn.index
  sby_desired_minority_k      <- min(
    as.integer(sby_knn_over_k),
    nrow(sby_minority_matrix) - 1L
  )

  # Verifica se a consulta minoritaria retornou mais de um vizinho candidato
  if(sby_effective_minority_k > 1L){

    # Remove autorreferencias da matriz de vizinhos minoritarios
    sby_minority_neighbor_index <- sby_drop_self_neighbor_index(
      sby_neighbor_index = sby_minority_neighbor_index,
      sby_self_index     = seq_len(nrow(sby_minority_matrix)),
      sby_desired_k      = sby_desired_minority_k
    )
  }

  # Gera amostras sinteticas por rotina nativa quando disponivel
  if(sby_adanear_native_available()){

    # Garante modos de armazenamento esperados pela rotina nativa
    storage.mode(sby_minority_matrix)         <- "double"
    storage.mode(sby_minority_neighbor_index) <- "integer"

    # Escolha de kernel ADASYN nativo.
    #
    # Default = "row": kernel original linha-a-linha (OU_GenerateSyntheticAdasynC).
    # Vence em testes empiricos com p moderado (p <= 200) e n_synthetic moderado
    # porque tem overhead minimo: nao pre-aloca vetores temporarios e gera +
    # escreve cada sintetico em uma unica passada.
    #
    # Opcao "col": kernel column-friendly (OU_GenerateSyntheticAdasynColC), que
    # pre-resolve (baseRow, nbrRow, weight) em vetores contiguos e percorre as
    # colunas no laco externo. Em teoria favorece locality em dimensionalidade
    # muito alta (p >> 200) e n_synthetic alto (> 10^5). Em casos pequenos a
    # pre-alocacao adicional anula o ganho.
    #
    # Automatic selection: prefer row kernel as default stable path.
    sby_adasyn_kernel <- "row"
    if(identical(sby_adasyn_kernel, "col")){
      sby_synthetic_matrix <- .Call(
        OU_GenerateSyntheticAdasynColC,
        sby_minority_matrix,
        sby_minority_neighbor_index,
        as.integer(sby_synthetic_per_row)
      )
    }else{
      sby_synthetic_matrix <- .Call(
        OU_GenerateSyntheticAdasynC,
        sby_minority_matrix,
        sby_minority_neighbor_index,
        as.integer(sby_synthetic_per_row)
      )
    }

    # Verifica se ha solicitacao de interrupcao apos geracao nativa
    sby_adanear_check_user_interrupt()
  }else{

    # Inicializa matriz sintetica e ponteiros de escrita para geracao em R
    sby_synthetic_matrix <- matrix(
      data = 0,
      nrow = sby_synthetic_count,
      ncol = NCOL(sby_x_scaled)
    )
    sby_write_start   <- 1L
    sby_positive_rows <- which(
      x = sby_synthetic_per_row > 0L
    )

    # Gera amostras sinteticas por interpolacao nas linhas positivas
    for(i in sby_positive_rows){
      # Verifica interrupcao a cada linha minoritaria sintetizada
      sby_adanear_check_user_interrupt()

      # Calcula intervalo de escrita para a linha minoritaria corrente
      sby_row_count   <- sby_synthetic_per_row[[i]]
      sby_write_end   <- sby_write_start + sby_row_count - 1L
      sby_base_rows   <- matrix(
        data = sby_minority_matrix[i, ],
        nrow = sby_row_count,
        ncol = NCOL(sby_x_scaled),
        byrow = TRUE
      )

      # Amostra vizinhos minoritarios para interpolacao local
      sby_selected_neighbor_rows <- sby_minority_neighbor_index[i, sample.int(
        n = ncol(sby_minority_neighbor_index),
        size = sby_row_count,
        replace = TRUE
      )]
      sby_neighbor_rows <- sby_minority_matrix[sby_selected_neighbor_rows, , drop = FALSE]

      # Preenche bloco sintetico por interpolacao aleatoria entre base e vizinho
      sby_synthetic_matrix[sby_write_start:sby_write_end, ] <- sby_base_rows + stats::runif(
        n = sby_row_count
      ) * (sby_neighbor_rows - sby_base_rows)

      # Avanca ponteiro de escrita para o proximo bloco sintetico
      sby_write_start <- sby_write_end + 1L
    }
  }

  # Preserva nomes de colunas na matriz sintetica
  colnames(sby_synthetic_matrix) <- colnames(sby_x_scaled)

  # Retorna matriz expandida e alvo expandido com niveis preservados.
  # Evita o ciclo factor -> character -> factor (que aloca vetores de strings
  # do tamanho do conjunto inteiro): trabalhamos diretamente em codigos inteiros
  # e reconstruimos o fator no final.
  sby_minority_level_code <- match(
    sby_class_roles$sby_minority_label,
    levels(sby_target_factor)
  )
  sby_y_codes <- c(
    unclass(sby_target_factor),
    rep.int(sby_minority_level_code, sby_synthetic_count)
  )
  attributes(sby_y_codes) <- NULL
  storage.mode(sby_y_codes) <- "integer"
  sby_y_factor <- structure(
    sby_y_codes,
    levels = levels(sby_target_factor),
    class = "factor"
  )

  # Consolida a matriz expandida por kernel nativo quando possivel.
  # base::rbind() e generico e faz validacoes/despacho desnecessarios aqui;
  # neste ponto ambas as entradas ja sao matrizes double com o mesmo numero
  # de colunas. O caminho nativo aloca a saida uma unica vez e copia cada
  # coluna em blocos contiguos.
  if(sby_adanear_native_available()){
    storage.mode(sby_x_scaled) <- "double"
    storage.mode(sby_synthetic_matrix) <- "double"
    sby_expanded_x <- .Call(
      OU_RbindDoubleMatrixC,
      sby_x_scaled,
      sby_synthetic_matrix
    )
  }else{
    sby_expanded_x <- rbind(
      sby_x_scaled,
      sby_synthetic_matrix
    )
  }

  return(list(
    x = sby_expanded_x,
    y = sby_y_factor
  ))
}
####
## Fim
#
