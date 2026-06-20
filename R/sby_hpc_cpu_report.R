#' Relatorio de CPU e compilacao do motor HPC
#'
#' @description
#' `sby_hpc_cpu_report()` consolida evidencias objetivas para confirmar se o
#' servidor e o pacote estao prontos para usar o caminho Intel Cascade Lake com
#' AVX-512, AVX2, FMA e OpenMP.
#'
#' @details
#' A confirmacao fica dividida em tres camadas:
#' * `runtime_flags`: instrucoes anunciadas pelo processador em `/proc/cpuinfo`.
#' * `compile_report`: macros gravadas no binario nativo durante a compilacao.
#' * `binary_scan`: busca opcional por mnemonicos AVX-512/FMA no `.so` com
#'   `objdump`, quando a ferramenta esta instalada no sistema.
#'
#' Para a meta do cliente, considere OK quando `runtime_ok`, `compile_ok` e
#' `openmp_ok` forem `TRUE`. `binary_scan$has_zmm_or_avx512` e uma evidencia
#' adicional forte de que o objeto contem instrucoes AVX-512, mas pode ficar
#' `NA` em ambientes sem `objdump`.
#'
#' @return Lista com indicadores logicos e evidencias de runtime, compilacao,
#' threads MKL/OpenMP e varredura opcional do binario nativo.
#'
#' @export
sby_hpc_cpu_report <- function(){
  sby_flags <- character()
  if(file.exists("/proc/cpuinfo")){
    sby_cpuinfo <- readLines("/proc/cpuinfo", warn = FALSE)
    sby_flag_lines <- grep("^(flags|Features)[[:space:]]*:", sby_cpuinfo, value = TRUE)
    sby_flags <- unique(unlist(strsplit(tolower(paste(sby_flag_lines, collapse = " ")), "[[:space:]]+")))
    sby_flags <- setdiff(sby_flags, c("flags", "features", ":"))
  }

  sby_required_runtime <- c("sse4_2", "avx", "avx2", "avx512f", "avx512cd", "avx512bw", "avx512dq", "avx512vl", "fma")
  sby_runtime_flags <- stats::setNames(sby_required_runtime %in% sby_flags, sby_required_runtime)

  sby_compile_report <- tryCatch(
    sby_call_native("sby_hpc_compile_report_cpp"),
    error = function(sby_error) list(error = conditionMessage(sby_error))
  )

  sby_compile_required <- c("cascade_lake_native", "avx512f", "avx512cd", "avx512bw", "avx512dq", "avx512vl", "avx2", "fma")
  sby_compile_ok <- if(is.list(sby_compile_report) && is.null(sby_compile_report$error)){
    all(vapply(sby_compile_required, function(sby_key) isTRUE(sby_compile_report[[sby_key]]), logical(1)))
  }else{
    FALSE
  }

  sby_openmp_ok <- is.list(sby_compile_report) && isTRUE(sby_compile_report$openmp)

  sby_lib_path <- tryCatch(system.file("libs", package = "sbyadanear"), error = function(sby_error) "")
  sby_shared_object <- character()
  if(nzchar(sby_lib_path) && dir.exists(sby_lib_path)){
    sby_shared_object <- list.files(sby_lib_path, pattern = paste0(.Platform$dynlib.ext, "$"), full.names = TRUE)
  }

  sby_objdump <- Sys.which("objdump")
  sby_binary_scan <- list(
    checked = FALSE,
    shared_object = if(length(sby_shared_object)) sby_shared_object[[1]] else NA_character_,
    has_zmm_or_avx512 = NA,
    has_vfmadd = NA,
    note = "objdump indisponivel ou biblioteca nativa nao localizada"
  )
  if(nzchar(sby_objdump) && length(sby_shared_object)){
    sby_disasm <- tryCatch(
      system2(sby_objdump, c("-d", sby_shared_object[[1]]), stdout = TRUE, stderr = FALSE),
      error = function(sby_error) character()
    )
    sby_binary_text <- tolower(paste(sby_disasm, collapse = "\n"))
    sby_binary_scan <- list(
      checked = length(sby_disasm) > 0L,
      shared_object = sby_shared_object[[1]],
      has_zmm_or_avx512 = grepl("\\bzmm[0-9]+\\b|avx512", sby_binary_text),
      has_vfmadd = grepl("vfmadd", sby_binary_text),
      note = if(length(sby_disasm) > 0L) "varredura objdump concluida" else "objdump nao retornou desassembly"
    )
  }

  sby_env_keys <- c("MKL_NUM_THREADS", "OMP_NUM_THREADS")

  list(
    runtime_ok = all(sby_runtime_flags),
    compile_ok = sby_compile_ok,
    openmp_ok = sby_openmp_ok,
    runtime_flags = sby_runtime_flags,
    compile_report = sby_compile_report,
    binary_scan = sby_binary_scan,
    hpc_environment = Sys.getenv(sby_env_keys, unset = NA),
    guidance = c(
      "runtime_ok confirma que a CPU anuncia AVX-512/AVX2/FMA.",
      "compile_ok confirma que o pacote foi compilado com -march=cascadelake e macros AVX-512.",
      "openmp_ok confirma que os loops Fortran/C++ foram ligados com OpenMP.",
      "binary_scan com zmm/vfmadd TRUE e a prova adicional no binario nativo."
    )
  )
}
