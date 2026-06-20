#' @title Preparar preditores com z-score populacional
#'
#' @usage sby_mutate_zscore(sby_formula, sby_data, sby_engine = c("native", "r"))
#'
#' @description
#' Extrai os preditores numéricos definidos por uma fórmula, calcula parâmetros
#' de padronização populacional e devolve a matriz transformada junto com os
#' insumos necessários para rotinas de balanceamento.
#'
#' @details
#' A função é interna e preserva uma cópia da matriz original para que registros
#' sintéticos possam ser revertidos à escala observada quando necessário. Na
#' engine `"native"`, a rotina delega cálculo e aplicação do z-score ao núcleo
#' Fortran registrado via `.Call`; na engine `"r"`, usa os símbolos nativos já
#' existentes para manter a mesma convenção de escala.
#'
#' A padronização usa denominador populacional, isto é:
#'
#' $$z_{ij} = \frac{x_{ij} - \mu_j}{\sigma_j}, \quad
#' \sigma_j = \sqrt{\frac{1}{n}\sum_{i = 1}^{n}(x_{ij} - \mu_j)^2}$$
#'
#' **Fluxo interno:**
#'
#' ```mermaid
#' flowchart LR
#'   A[Fórmula e dados] --> B[Extração de alvo e preditores]
#'   B --> C[Matriz double densa]
#'   C --> D{Engine native?}
#'   D -->|Sim| E[Fortran registrado]
#'   D -->|Não| F[Símbolos nativos auxiliares]
#'   E --> G[Lista técnica de z-score]
#'   F --> G
#' ```
#'
#' > **Nota:** colunas com desvio zero devem ser validadas pelas camadas
#' > chamadoras, pois a função documenta e executa apenas a etapa de preparo.
#'
#' @param sby_formula Fórmula no formato `alvo ~ preditores`.
#'
#' @param sby_data Data frame ou tibble com a coluna alvo e os preditores.
#'
#' @param sby_engine Engine de padronização, aceita `"native"` ou `"r"`.
#'
#' @return Lista com `x_scaled`, `x_original`, médias, desvios, nomes de
#'   features, dimensões e vetor alvo extraído.
#'
#' @section Pré-condições:
#' Os preditores extraídos devem ser conversíveis para matriz `double` densa.
#'
#' @section Pós-condições:
#' A matriz escalonada mantém a mesma ordem de linhas e colunas da matriz
#' original usada pela fórmula.
#'
#' @seealso sby_extract_formula_data, sby_call_native, sby_configure_blas_threads
#'
#' @references He, H., Bai, Y., Garcia, E. A., & Li, S. (2008). ADASYN.
#'
#' @examples
#' \dontrun{
#' sbyadanear:::sby_mutate_zscore(target ~ ., dados, sby_engine = "native")
#' }
#'
#' @keywords internal
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
    sby_previous_blas_env <- sby_configure_blas_threads(sby_workers = 1L)
    on.exit(sby_restore_blas_threads(sby_previous_blas_env), add = TRUE)
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
