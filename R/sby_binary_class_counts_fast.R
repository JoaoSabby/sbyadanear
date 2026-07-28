#' Contar classes binarias rapidamente
#'
#' @param sby_y_vector Fator binario sem ausentes
#'
#' @return Lista com contagens nomeadas, papeis de classe, codigos inteiros de
#' nivel (`sby_minority_level`/`sby_majority_level`) e razao maioria/minoria
#'
#' @noRd
sby_binary_class_counts_fast <- function(sby_y_vector){
  if(!is.factor(sby_y_vector)){
    sby_adanear_abort(
      sby_message = "'sby_y_vector' deve ser factor"
    )
  }

  if(anyNA(sby_y_vector)){
    sby_adanear_abort(
      sby_message = "'sby_y_vector' nao pode conter NA"
    )
  }

  sby_levels <- levels(sby_y_vector)
  if(length(sby_levels) != 2L){
    sby_adanear_abort(
      sby_message = "'sby_y_vector' deve possuir exatamente dois niveis"
    )
  }

  sby_counts <- tabulate(
    bin = as.integer(sby_y_vector),
    nbins = length(sby_levels)
  )
  names(sby_counts) <- sby_levels

  if(any(sby_counts < 1L)){
    sby_adanear_abort(
      sby_message = "'sby_y_vector' deve conter observacoes em ambos os niveis"
    )
  }

  # Com contagens empatadas, which.min() e which.max() devolvem o mesmo nivel e
  # os papeis de classe colapsam em silencio. Aborta com a mesma mensagem de
  # sby_get_binary_class_roles() para que as duas rotas tenham o mesmo contrato.
  if(sby_counts[[1L]] == sby_counts[[2L]]){
    sby_adanear_abort(
      sby_message = "As rotinas de sampling requerem classes desbalanceadas"
    )
  }

  sby_minority_position <- which.min(sby_counts)
  sby_majority_position <- which.max(sby_counts)

  return(list(
    sby_class_counts = sby_counts,
    sby_minority_label = sby_levels[[sby_minority_position]],
    sby_majority_label = sby_levels[[sby_majority_position]],
    # Codigos inteiros de nivel: as rotas HPC indexam o factor por codigo, do
    # mesmo modo que sby_resolve_minority_role() no motor C++.
    sby_minority_level = as.integer(sby_minority_position),
    sby_majority_level = as.integer(sby_majority_position),
    sby_minority_count = as.integer(sby_counts[[sby_minority_position]]),
    sby_majority_count = as.integer(sby_counts[[sby_majority_position]]),
    sby_class_ratio = as.numeric(sby_counts[[sby_majority_position]] / sby_counts[[sby_minority_position]])
  ))
}
