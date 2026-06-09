#' Capturar, injetar e restaurar o ambiente NUMA e MKL para a rota HPC
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato
#' de entrada explicito e retorno controlado. As rotinas abaixo isolam as variaveis
#' de ambiente que controlam a distribuicao NUMA e o comportamento do MKL antes de
#' acionar o motor HPC, e fornecem um plano de restauro inflexivel para devolver o
#' ambiente ao estado anterior, removendo o que nao existia.
#'
#' @return Lista com o estado anterior das variaveis modificadas.
#' @noRd

# Lista canonica das variaveis de ambiente controladas pela rota HPC
sby_hpc_env_keys <- function(){
  c(
    "KMP_AFFINITY",
    "MKL_NUM_STRIPES",
    "MKL_DISABLE_FAST_MM",
    "MKL_DYNAMIC",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS"
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

# Injeta as variaveis NUMA e MKL para distribuir o uso das FMA de forma simetrica
sby_hpc_apply_env <- function(
  sby_total_threads,
  sby_affinity = c("scatter", "balanced")
){
  sby_affinity <- match.arg(sby_affinity)
  sby_total_threads <- max(1L, as.integer(sby_total_threads))

  # MKL_NUM_STRIPES segue a heuristica numero total de threads dividido por 4
  sby_num_stripes <- max(1L, as.integer(floor(sby_total_threads / 4L)))

  Sys.setenv(
    KMP_AFFINITY = paste0("granularity=fine,", sby_affinity),
    MKL_NUM_STRIPES = as.character(sby_num_stripes),
    MKL_DISABLE_FAST_MM = "1",
    MKL_DYNAMIC = "FALSE",
    MKL_NUM_THREADS = as.character(sby_total_threads),
    OMP_NUM_THREADS = as.character(sby_total_threads)
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
