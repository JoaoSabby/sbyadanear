#' Validar semente pseudoaleatoria
#'
#' @details
#' A funcao concentra a validacao de sementes usadas em `set.seed()` para evitar
#' truncamentos silenciosos, valores infinitos e entradas fora do intervalo de
#' inteiros aceito pela geracao pseudoaleatoria do R.
#'
#' @param sby_seed Valor informado como semente
#'
#' @return Inteiro escalar validado para uso em `set.seed()`
#'
#' @noRd
sby_validate_seed <- function(sby_seed){
  # Verifica tipo, tamanho e finitude antes de qualquer conversao para inteiro
  if(!(is.numeric(sby_seed) && length(sby_seed) == 1L && !is.na(sby_seed) && is.finite(sby_seed))){

    # Aborta quando a semente nao atende ao contrato esperado
    sby_adanear_abort(
      sby_message = "'sby_seed' deve ser escalar numerico finito"
    )
  }

  # Evita truncamento silencioso de valores decimais por set.seed/as.integer
  if(sby_seed != floor(sby_seed)){

    # Aborta quando a semente nao e representavel como inteiro
    sby_adanear_abort(
      sby_message = "'sby_seed' deve ser numero inteiro"
    )
  }

  # Garante faixa segura para conversao ao tipo integer do R
  if(sby_seed < 0 || sby_seed > .Machine$integer.max){

    # Aborta quando a semente esta fora da faixa aceita
    sby_adanear_abort(
      sby_message = "'sby_seed' deve estar entre 0 e .Machine$integer.max"
    )
  }

  # Retorna semente normalizada como inteiro escalar
  return(as.integer(sby_seed))
}
####
## Fim
#
