#' Remover o proprio ponto de uma matriz de vizinhos de forma robusta
#'
#' @details
#' A entrada e uma matriz de indices retornada por uma consulta KNN que pode
#' incluir a propria linha (quando query = data). A funcao retorna uma matriz
#' com `sby_desired_k` colunas em que a auto-referencia foi removida. A
#' implementacao e vetorizada: percorre apenas as linhas em que o self
#' efetivamente aparece, evitando o custo de um loop em R por linha. Para
#' bases com poucos empates de self isso reduz drasticamente o overhead em
#' comparacao a implementacao linha-a-linha anterior.
#'
#' @param sby_neighbor_index Matriz de indices de vizinhos candidatos
#' @param sby_self_index Vetor com indices das proprias linhas
#' @param sby_desired_k Numero desejado de vizinhos apos a remocao
#'
#' @return Matriz de indices de vizinhos sem o proprio ponto
#' @noRd
sby_drop_self_neighbor_index <- function(sby_neighbor_index, sby_self_index, sby_desired_k){
  sby_n <- nrow(sby_neighbor_index)
  sby_k_plus <- ncol(sby_neighbor_index)
  sby_desired_k <- as.integer(sby_desired_k)
  sby_self_index <- as.integer(sby_self_index)

  if(sby_desired_k < 1L){
    return(matrix(integer(0), nrow = sby_n, ncol = 0L))
  }

  # Caso degenerado: nao ha colunas suficientes para garantir k vizinhos
  # mesmo se nao houvesse self - retorna NA matriz para que o caller falhe
  # de forma controlada.
  if(sby_k_plus < sby_desired_k){
    sby_adanear_abort(
      sby_message = "Nao foi possivel remover o proprio ponto mantendo vizinhos suficientes"
    )
  }

  # Caminho nativo (rapido): kernel C com varredura linear por linha em vez
  # do baseline R + repair loop. Em bases grandes (n_minority >> 10^4) reduz
  # o overhead R consideravelmente.
  if(sby_adanear_native_available()){
    storage.mode(sby_neighbor_index) <- "integer"
    return(.Call(
      OU_DropSelfNeighborC,
      sby_neighbor_index,
      sby_self_index,
      sby_desired_k
    ))
  }

  # Caminho rapido vetorizado: assume que a coluna 1 do KNN, quando query=data,
  # contem o proprio ponto na maioria das linhas (caso tipico de FNN exato).
  # Estrategia em duas fases:
  #   1) Caso geral: dropa a ultima coluna -> baseline com k colunas.
  #   2) Repara apenas as linhas em que o self ainda esta presente nas k
  #      primeiras colunas.
  sby_out <- sby_neighbor_index[, seq_len(sby_desired_k), drop = FALSE]

  # Identifica linhas que ainda contem self dentro do bloco baseline.
  sby_baseline_self <- sby_out == sby_self_index
  sby_baseline_self[is.na(sby_baseline_self)] <- FALSE
  sby_rows_to_fix <- which(.rowSums(sby_baseline_self, m = sby_n, n = sby_desired_k) > 0L)

  if(length(sby_rows_to_fix) == 0L){
    # Mas tambem precisamos confirmar que o self nao esta na ultima coluna
    # (descartada): se estiver, baseline ja esta correto. Se nao estava em
    # nenhuma das k+1 colunas, baseline ainda esta correto.
    return(sby_out)
  }

  # Para as linhas que precisam de conserto, reconstroi vizinhos validos por
  # linha. Loop limitado as linhas problematicas, nao a base inteira.
  for(i in sby_rows_to_fix){
    sby_candidates <- sby_neighbor_index[i, ]
    sby_candidates <- sby_candidates[!is.na(sby_candidates) & sby_candidates != sby_self_index[[i]]]
    if(length(sby_candidates) < sby_desired_k){
      sby_adanear_abort(
        sby_message = "Nao foi possivel remover o proprio ponto mantendo vizinhos suficientes"
      )
    }
    sby_out[i, ] <- sby_candidates[seq_len(sby_desired_k)]
  }

  return(sby_out)
}
####
## Fim
#
