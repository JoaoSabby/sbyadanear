# Padronizacao z-score populacional via engine nativa Fortran
#
# Calcula e aplica z-score populacional sobre preditores numericos
# usando o motor nativo em Fortran com paralelismo OpenMP. Opera
# sobre a matriz transposta (P x N) para maximizar localidade de
# cache e vetorizacao AVX-512.
#
# @param sby_formula Formula no formato alvo ~ preditores.
# @param sby_data Data frame ou tibble com os dados.
# @param sby_engine character. Engine a usar: "native" (Fortran) ou "r" (base).
#
# @return Lista com matriz padronizada, vetores de media e desvio,
#   e nomes das features.
# @keywords internal
sby_mutate_zscore <- function(
    sby_formula,
    sby_data,
    sby_engine = c("native", "r")
){
  sby_engine <- match.arg(sby_engine)

  # Extrai preditores e target usando a infraestrutura existente
  sby_formula_data <- sby_extract_formula_data(
    sby_formula = sby_formula,
    sby_data = sby_data
  )

  sby_predictor_data <- sby_formula_data$sby_predictor_data
  feature_names <- colnames(sby_predictor_data)
  n <- collapse::fnrow(sby_predictor_data)
  p <- collapse::fncol(sby_predictor_data)

  # Converte para matriz numerica densa
  x_matrix <- as.matrix(sby_predictor_data)
  storage.mode(x_matrix) <- "double"

  if (sby_engine == "native") {
    sby_configure_blas_threads(sby_workers = 1L)
    # Usa layout R n x p diretamente para evitar copias por transposicao
    zscore_params <- sby_call_native(
      "compute_zscore_population_fortran_c",
      x_matrix
    )

    # Aplica z-score via Fortran
    x_scaled <- sby_call_native(
      "apply_zscore_fortran_c",
      x_matrix,
      zscore_params$means,
      zscore_params$sds
    )

    return(list(
      x_scaled = x_scaled,
      x_original = x_matrix,
      means = zscore_params$means,
      sds = zscore_params$sds,
      feature_names = feature_names,
      n = n,
      p = p,
      target_vector = sby_formula_data$sby_target_vector
    ))
  }

  # Fallback R puro: usa as rotinas C existentes
  zscore_params <- sby_call_native(
    "compute_z_score_params_c",
    x_matrix
  )
  x_scaled <- sby_call_native(
    "apply_z_score_c",
    x_matrix,
    zscore_params$centers,
    zscore_params$scales,
    FALSE
  )

  return(list(
    x_scaled = x_scaled,
    x_original = x_matrix,
    means = zscore_params$centers,
    sds = zscore_params$scales,
    feature_names = feature_names,
    n = n,
    p = p,
    target_vector = sby_formula_data$sby_target_vector
  ))
}
