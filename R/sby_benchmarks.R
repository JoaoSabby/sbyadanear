# Rotinas de benchmark para a engine nativa do sbyadanear
#
# Funcoes de benchmark para medir tempo de execucao e pico de
# memoria residente (RSS) das rotinas criticas do pacote.
# Os resultados sao retornados como tibbles para facil analise.
#
# @keywords internal

# Le o pico de memoria residente (RSS) do processo atual via procfs (Linux).
# Retorna o valor em megabytes. Em sistemas nao-Linux, retorna NA.
sby_benchmark_peak_rss_mb <- function() {
  status_file <- "/proc/self/status"

  if (!file.exists(status_file)) {
    return(NA_real_)
  }

  lines <- readLines(status_file, warn = FALSE)
  vmrss_line <- grep("^VmRSS:", lines, value = TRUE)

  if (length(vmrss_line) == 0L) {
    return(NA_real_)
  }

  # Extrai valor numerico em kB e converte para MB
  rss_kb <- as.numeric(gsub("[^0-9]", "", vmrss_line[1L]))
  return(rss_kb / 1024)
}


# Executa benchmark de uma expressao N vezes e retorna tibble com estatisticas.
#
# @param expr Expressao a ser avaliada (nao-avaliada via substitute).
# @param n_reps integer. Numero de repeticoes.
# @param label character. Rotulo descritivo do benchmark.
#
# @return tibble com label, min_time, max_time, median_time, peak_rss_mb.
sby_benchmark_run <- function(expr, n_reps = 5L, label = "benchmark") {
  expr_q <- substitute(expr)
  times <- numeric(n_reps)

  for (i in seq_len(n_reps)) {
    gc(verbose = FALSE, full = TRUE)
    start_time <- proc.time()["elapsed"]
    eval(expr_q, envir = parent.frame())
    end_time <- proc.time()["elapsed"]
    times[i] <- end_time - start_time
  }

  peak_rss <- sby_benchmark_peak_rss_mb()

  tibble::tibble(
    label      = label,
    n_reps     = n_reps,
    min_time   = min(times),
    max_time   = max(times),
    median_time = stats::median(times),
    mean_time  = collapse::fmean(times),
    peak_rss_mb = peak_rss
  )
}


# Benchmark de z-score populacional nativo.
#
# @param sby_formula Formula alvo ~ preditores.
# @param sby_data Data frame ou tibble.
# @param n_reps integer. Numero de repeticoes.
#
# @return tibble com metricas de tempo e memoria.
sby_mutate_zscore_benchmark <- function(
    sby_formula,
    sby_data,
    n_reps = 5L
){
  sby_benchmark_run(
    sby_mutate_zscore(sby_formula, sby_data, sby_engine = "native"),
    n_reps = n_reps,
    label = "sby_mutate_zscore_native"
  )
}


# Benchmark de ADASYN completo.
#
# @param sby_formula Formula alvo ~ preditores.
# @param sby_data Data frame ou tibble.
# @param n_reps integer. Numero de repeticoes.
# @param ... Argumentos extras para sby_adasyn().
#
# @return tibble com metricas de tempo e memoria.
sby_adasyn_benchmark <- function(
    sby_formula,
    sby_data,
    n_reps = 5L,
    ...
){
  sby_benchmark_run(
    sby_adasyn(sby_formula, sby_data, ...),
    n_reps = n_reps,
    label = "sby_adasyn"
  )
}


# Benchmark de NearMiss completo.
#
# @param sby_formula Formula alvo ~ preditores.
# @param sby_data Data frame ou tibble.
# @param n_reps integer. Numero de repeticoes.
# @param ... Argumentos extras para sby_nearmiss().
#
# @return tibble com metricas de tempo e memoria.
sby_nearmiss_benchmark <- function(
    sby_formula,
    sby_data,
    n_reps = 5L,
    ...
){
  sby_benchmark_run(
    sby_nearmiss(sby_formula, sby_data, ...),
    n_reps = n_reps,
    label = "sby_nearmiss"
  )
}


# Benchmark de pipeline hibrido adanear completo.
#
# @param sby_formula Formula alvo ~ preditores.
# @param sby_data Data frame ou tibble.
# @param n_reps integer. Numero de repeticoes.
# @param ... Argumentos extras para sby_adanear().
#
# @return tibble com metricas de tempo e memoria.
sby_adanear_benchmark <- function(
    sby_formula,
    sby_data,
    n_reps = 5L,
    ...
){
  sby_benchmark_run(
    sby_adanear(sby_formula, sby_data, ...),
    n_reps = n_reps,
    label = "sby_adanear"
  )
}


# Benchmark comparativo multi-engine
#
# Executa o mesmo pipeline com todas as engines KNN disponiveis
# e retorna um tibble unico para comparacao direta. Cada engine
# e testada isoladamente com gc() completo entre rodadas.
#
# @param sby_formula Formula alvo ~ preditores.
# @param sby_data Data frame ou tibble.
# @param sby_pipeline character. Pipeline a testar: "adasyn", "nearmiss" ou "adanear".
# @param n_reps integer. Numero de repeticoes por engine.
# @param sby_engines character vector. Engines a testar. O padrao testa todas.
# @param ... Argumentos extras passados para a funcao de pipeline.
#
# @return tibble com colunas: engine, pipeline, n_rows, n_cols, n_reps,
#   min_sec, median_sec, mean_sec, max_sec, peak_rss_mb.
# @export
sby_benchmark_engines <- function(
    sby_formula,
    sby_data,
    sby_pipeline = c("adasyn", "nearmiss", "adanear"),
    n_reps = 3L,
    sby_engines = c("native", "FNN", "RcppHNSW", "KernelKnn", "bigKNN"),
    ...
){
  sby_pipeline <- match.arg(sby_pipeline)

  # Resolve funcao de pipeline
  pipeline_fn <- switch(sby_pipeline,
    "adasyn"  = sby_adasyn,
    "nearmiss" = sby_nearmiss,
    "adanear" = sby_adanear
  )

  n_rows <- collapse::fnrow(sby_data)
  n_cols <- collapse::fncol(sby_data) - 1L  # exclui target

  # Detecta engines disponiveis no ambiente atual
  available_engines <- character(0L)
  for (eng in sby_engines) {
    is_available <- tryCatch({
      if (eng == "native") {
        # Engine nativa sempre disponivel se o pacote compilou
        is.loaded("brute_force_knn_c")
      } else if (eng == "FNN") {
        requireNamespace("FNN", quietly = TRUE)
      } else if (eng == "RcppHNSW") {
        requireNamespace("RcppHNSW", quietly = TRUE)
      } else if (eng == "KernelKnn") {
        requireNamespace("KernelKnn", quietly = TRUE)
      } else if (eng == "bigKNN") {
        requireNamespace("bigKNN", quietly = TRUE) &&
          requireNamespace("bigmemory", quietly = TRUE)
      } else {
        FALSE
      }
    }, error = function(e) FALSE)

    if (isTRUE(is_available)) {
      available_engines <- c(available_engines, eng)
    } else {
      message("Engine '", eng, "' nao disponivel, pulando.")
    }
  }

  if (length(available_engines) == 0L) {
    stop("Nenhuma engine KNN disponivel para benchmark.")
  }

  results_list <- vector("list", length(available_engines))

  for (idx in seq_along(available_engines)) {
    eng <- available_engines[idx]
    message("Benchmarking engine: ", eng, " (", idx, "/", length(available_engines), ")")

    times <- numeric(n_reps)
    rss_values <- numeric(n_reps)

    for (rep_i in seq_len(n_reps)) {
      gc(verbose = FALSE, full = TRUE)
      rss_before <- sby_benchmark_peak_rss_mb()

      start_time <- proc.time()["elapsed"]

      tryCatch({
        result <- pipeline_fn(
          sby_formula,
          sby_data,
          sby_knn_engine = eng,
          ...
        )
      }, error = function(e) {
        message("  Erro na engine '", eng, "' rep ", rep_i, ": ", conditionMessage(e))
        times[rep_i] <<- NA_real_
        rss_values[rep_i] <<- NA_real_
        return(NULL)
      })

      end_time <- proc.time()["elapsed"]
      times[rep_i] <- end_time - start_time
      rss_values[rep_i] <- sby_benchmark_peak_rss_mb()
    }

    # Remove NAs de runs com erro
    valid_times <- times[!is.na(times)]
    valid_rss <- rss_values[!is.na(rss_values)]

    results_list[[idx]] <- tibble::tibble(
      engine      = eng,
      pipeline    = sby_pipeline,
      n_rows      = n_rows,
      n_cols      = n_cols,
      n_reps      = length(valid_times),
      min_sec     = if (length(valid_times) > 0L) min(valid_times) else NA_real_,
      median_sec  = if (length(valid_times) > 0L) stats::median(valid_times) else NA_real_,
      mean_sec    = if (length(valid_times) > 0L) collapse::fmean(valid_times) else NA_real_,
      max_sec     = if (length(valid_times) > 0L) max(valid_times) else NA_real_,
      peak_rss_mb = if (length(valid_rss) > 0L) max(valid_rss) else NA_real_
    )
  }

  # Empilha todos os resultados em um unico tibble
  benchmark_table <- do.call(rbind, results_list)
  benchmark_table <- tibble::as_tibble(benchmark_table)

  # Ordena por tempo mediano crescente
  benchmark_table <- benchmark_table[order(benchmark_table$median_sec), ]

  message("\nBenchmark concluido. ", length(available_engines), " engines testadas.")
  return(benchmark_table)
}

