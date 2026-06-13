#' Capturar, injetar e restaurar threads MKL e OpenMP para a rota HPC
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato
#' de entrada explicito e retorno controlado. As rotinas abaixo isolam as
#' variaveis de paralelismo antes de acionar o motor HPC, e fornecem um plano de
#' restauro inflexivel para devolver o ambiente ao estado anterior, removendo as
#' variaveis que nao existiam.
#'
#' @return Lista com o estado anterior das variaveis modificadas.
#' @noRd

# Lista canonica das variaveis de ambiente controladas pela rota HPC
sby_hpc_env_keys <- function(){
  c(
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_STRIPES"
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

# Resolve uma sugestao conservadora para MKL_NUM_STRIPES em GEMM column-major.
sby_hpc_resolve_mkl_num_stripes <- function(
  sby_total_threads,
  sby_majority_count = NA_integer_,
  sby_minority_count = NA_integer_
){
  sby_total_threads <- max(1L, suppressWarnings(as.integer(sby_total_threads)))
  if(sby_total_threads < 2L){
    return(1L)
  }

  sby_majority_count <- suppressWarnings(as.integer(sby_majority_count))
  sby_minority_count <- suppressWarnings(as.integer(sby_minority_count))

  if(length(sby_majority_count) == 1L && length(sby_minority_count) == 1L &&
     !is.na(sby_majority_count) && !is.na(sby_minority_count) &&
     sby_majority_count > 0L && sby_minority_count > 0L){
    sby_shape_ratio <- sby_majority_count / sby_minority_count
    if(sby_shape_ratio >= 2){
      return(max(1L, as.integer(ceiling(sby_total_threads / 2))))
    }
    if(sby_shape_ratio <= 0.5){
      return(1L)
    }
  }

  sby_2d_limit <- max(1L, as.integer(floor((sby_total_threads - 1L) / 2L)))
  max(1L, min(sby_2d_limit, as.integer(ceiling(sqrt(sby_total_threads)))))
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

# Injeta numero de threads MKL/OpenMP e uma sugestao de stripes para GEMM
sby_hpc_apply_env <- function(
  sby_total_threads,
  sby_majority_count = NA_integer_,
  sby_minority_count = NA_integer_
){
  sby_total_threads <- max(1L, as.integer(sby_total_threads))
  sby_num_stripes <- sby_hpc_resolve_mkl_num_stripes(
    sby_total_threads = sby_total_threads,
    sby_majority_count = sby_majority_count,
    sby_minority_count = sby_minority_count
  )

  sby_temporary_env <- c(
    MKL_NUM_THREADS = as.character(sby_total_threads),
    OMP_NUM_THREADS = as.character(sby_total_threads),
    MKL_NUM_STRIPES = as.character(sby_num_stripes)
  )

  do.call(Sys.setenv, as.list(sby_temporary_env))

  message(
    "sbyadanear HPC: variaveis de ambiente temporarias: ",
    paste(names(sby_temporary_env), sby_temporary_env, sep = "=", collapse = ", ")
  )

  invisible(sby_temporary_env)
}

# Restaura o ambiente original, removendo as variaveis que nao existiam antes
sby_hpc_restore_env <- function(sby_previous){
  if(is.null(sby_previous)){
    return(invisible(TRUE))
  }
  sby_restored_env <- character()
  for(sby_key in names(sby_previous)){
    sby_value <- sby_previous[[sby_key]]
    if(is.na(sby_value)){
      Sys.unsetenv(sby_key)
      sby_restored_env[[sby_key]] <- "<unset>"
    }else{
      sby_args <- stats::setNames(list(sby_value), sby_key)
      do.call(Sys.setenv, sby_args)
      sby_restored_env[[sby_key]] <- sby_value
    }
  }
  message(
    "sbyadanear HPC: variaveis de ambiente restauradas: ",
    paste(names(sby_restored_env), sby_restored_env, sep = "=", collapse = ", ")
  )
  invisible(sby_restored_env)
}
####
## Fim
#
