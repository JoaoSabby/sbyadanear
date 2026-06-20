#' Executar expressao com semente local preservando o RNG do chamador
#'
#' @details
#' Salva `RNGkind()` e `.Random.seed` antes de inicializar a semente local.
#' Ao sair, restaura o tipo de RNG e o estado aleatorio global exatamente como
#' estavam antes da chamada, inclusive quando `.Random.seed` nao existia.
#'
#' @param sby_seed Semente inteira validada para uso em `set.seed()`.
#'
#' @param sby_expr Expressao a executar com RNG local.
#'
#' @return Valor produzido por `sby_expr`.
#'
#' @keywords internal
sby_with_seed <- function(sby_seed, sby_expr){
  sby_rng_kind <- RNGkind()
  sby_had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  sby_random_seed <- if(sby_had_random_seed){
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }else{
    NULL
  }

  on.exit({
    do.call(RNGkind, as.list(sby_rng_kind))
    if(isTRUE(sby_had_random_seed)){
      assign(".Random.seed", sby_random_seed, envir = .GlobalEnv)
    }else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(sby_seed)
  force(sby_expr)
}
