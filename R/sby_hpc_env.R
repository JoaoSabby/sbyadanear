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
#'
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

# Resolve uma sugestao adaptativa para MKL_NUM_STRIPES em GEMM column-major.
sby_hpc_resolve_mkl_num_stripes <- function(
  sby_total_threads,
  sby_majority_count = NA_integer_,
  sby_minority_count = NA_integer_,
  sby_column_count = NA_integer_
){
  sby_total_threads <- suppressWarnings(as.integer(sby_total_threads)[1L])
  if(is.na(sby_total_threads) || sby_total_threads < 1L){
    sby_total_threads <- 1L
  }
  if(sby_total_threads < 2L){
    return(1L)
  }

  sby_majority_count <- suppressWarnings(as.integer(sby_majority_count)[1L])
  sby_minority_count <- suppressWarnings(as.integer(sby_minority_count)[1L])
  sby_column_count <- suppressWarnings(as.integer(sby_column_count)[1L])

  # Limite 2D seguro para matrizes sem informacao de forma: distribui o
  # paralelismo sem criar stripes demais em cargas pequenas ou balanceadas.
  sby_default_stripes <- max(1L, as.integer(ceiling(sqrt(sby_total_threads))))

  if(!is.na(sby_majority_count) && !is.na(sby_minority_count) &&
     sby_majority_count > 0L && sby_minority_count > 0L){
    sby_shape_ratio <- sby_majority_count / sby_minority_count

    # Heuristica conservadora para detectar matrizes FP32 que tendem a exceder
    # uma cache L3 compartilhada aproximada de 35 MB em topologias duplas.
    sby_is_large_matrix <- FALSE
    if(!is.na(sby_column_count) && sby_column_count > 0L){
      sby_matrix_bytes <- as.double(sby_majority_count) * as.double(sby_column_count) * 4
      sby_is_large_matrix <- sby_matrix_bytes > 35000000
    }

    if(sby_shape_ratio >= 2){
      if(sby_is_large_matrix || sby_shape_ratio > 100){
        return(sby_total_threads)
      }
      return(max(1L, as.integer(ceiling(sby_total_threads / 2))))
    }
    if(sby_shape_ratio <= 0.5){
      return(1L)
    }
  }

  sby_default_stripes
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
  sby_minority_count = NA_integer_,
  sby_column_count = NA_integer_
){
  sby_total_threads <- suppressWarnings(as.integer(sby_total_threads)[1L])
  if(is.na(sby_total_threads) || sby_total_threads < 1L){
    sby_total_threads <- 1L
  }
  sby_num_stripes <- sby_hpc_resolve_mkl_num_stripes(
    sby_total_threads = sby_total_threads,
    sby_majority_count = sby_majority_count,
    sby_minority_count = sby_minority_count,
    sby_column_count = sby_column_count
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
