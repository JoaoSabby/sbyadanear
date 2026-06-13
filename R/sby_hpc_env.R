#' Capturar, injetar e restaurar threads MKL e OpenMP para a rota HPC
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato
#' de entrada explicito e retorno controlado. As rotinas abaixo isolam somente
#' `MKL_NUM_THREADS` e `OMP_NUM_THREADS` antes de acionar o motor HPC, e fornecem
#' um plano de restauro inflexivel para devolver o ambiente ao estado anterior,
#' removendo as variaveis que nao existiam.
#'
#' @return Lista com o estado anterior das variaveis modificadas.
#' @noRd

# Lista canonica das variaveis de ambiente controladas pela rota HPC
sby_hpc_env_keys <- function(){
  c(
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OMP_PROC_BIND",
    "OMP_PLACES"
  )
}

# Resolve o numero de threads efetivo para o motor HPC
sby_hpc_resolve_threads <- function(sby_config_max_threads = -1L){
  sby_config_max_threads <- suppressWarnings(as.integer(sby_config_max_threads))
  if(length(sby_config_max_threads) != 1L || is.na(sby_config_max_threads)){
    sby_config_max_threads <- -1L
  }

  sby_detected <- tryCatch(
    as.integer(parallel::detectCores(logical = FALSE)),
    error = function(sby_error) NA_integer_
  )
  if(is.na(sby_detected) || sby_detected < 1L){
    sby_detected <- 1L
  }

  if(sby_config_max_threads > 0L){
    return(min(sby_config_max_threads, sby_detected))
  }
  return(sby_detected)
}

# Captura o estado anterior das variaveis controladas usando unset = NA
sby_hpc_capture_env <- function(){
  sby_keys <- sby_hpc_env_keys()
  sby_previous <- stats::setNames(
    lapply(sby_keys, function(sby_key) Sys.getenv(sby_key, unset = NA)),
    sby_keys
  )
  return(sby_previous)
}

# Injeta somente o numero de threads MKL e OpenMP para a rotina atual
sby_hpc_apply_env <- function(
  sby_total_threads
){
  sby_total_threads <- max(1L, as.integer(sby_total_threads))

  Sys.setenv(
    MKL_NUM_THREADS = as.character(sby_total_threads),
    OMP_NUM_THREADS = as.character(sby_total_threads),
    OMP_PROC_BIND = "spread",
    OMP_PLACES = "cores"
  )

  invisible(TRUE)
}

# Restaura o ambiente original, removendo as variaveis que nao existiam antes
sby_hpc_restore_env <- function(sby_previous){
  if(is.null(sby_previous)){
    return(invisible(TRUE))
  }
  for(sby_key in names(sby_previous)){
    sby_value <- sby_previous[[sby_key]]
    if(is.na(sby_value)){
      Sys.unsetenv(sby_key)
    }else{
      sby_args <- stats::setNames(list(sby_value), sby_key)
      do.call(Sys.setenv, sby_args)
    }
  }
  invisible(TRUE)
}
####
## Fim
#
