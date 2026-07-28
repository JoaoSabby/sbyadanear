#' Capturar e restaurar threads MKL e OpenMP para a rota HPC
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
    "OMP_NUM_THREADS"
  )
}

# Teto de CPUs imposto pelo cgroup corrente, em NA quando nao ha limite.
# parallel::detectCores() enxerga a maquina inteira e ignora containers e slices
# do systemd, o que faz o motor abrir dezenas de threads dentro de uma cota de
# poucos nucleos e degradar o desempenho por oversubscription.
sby_hpc_cgroup_cpu_quota <- function(){
  sby_read_first_line <- function(sby_path){
    if(!file.exists(sby_path)){
      return(NA_character_)
    }
    tryCatch(
      readLines(sby_path, n = 1L, warn = FALSE)[1L],
      error = function(sby_error) NA_character_
    )
  }

  # cgroup v2: "<quota|max> <period>"
  sby_v2 <- sby_read_first_line("/sys/fs/cgroup/cpu.max")
  if(!is.na(sby_v2)){
    sby_parts <- strsplit(trimws(sby_v2), "[[:space:]]+")[[1L]]
    if(length(sby_parts) == 2L && !identical(sby_parts[[1L]], "max")){
      sby_quota  <- suppressWarnings(as.numeric(sby_parts[[1L]]))
      sby_period <- suppressWarnings(as.numeric(sby_parts[[2L]]))
      if(!is.na(sby_quota) && !is.na(sby_period) && sby_quota > 0 && sby_period > 0){
        return(max(1L, as.integer(floor(sby_quota / sby_period))))
      }
    }
  }

  # cgroup v1: quota e period em arquivos separados, quota -1 quando ilimitada.
  sby_quota  <- suppressWarnings(as.numeric(
    sby_read_first_line("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
  ))
  sby_period <- suppressWarnings(as.numeric(
    sby_read_first_line("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
  ))
  if(!is.na(sby_quota) && !is.na(sby_period) && sby_quota > 0 && sby_period > 0){
    return(max(1L, as.integer(floor(sby_quota / sby_period))))
  }

  NA_integer_
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

  sby_quota <- sby_hpc_cgroup_cpu_quota()
  if(!is.na(sby_quota) && sby_quota >= 1L){
    sby_detected <- min(sby_detected, sby_quota)
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
