#' Listar dependencias comuns das etapas recipes de sampling
#'
#' @details
#' Centraliza a declaracao de pacotes usados pelas etapas para evitar divergencia
#' entre ADASYN, NearMiss e ADANEAR. Dependencias opcionais de engines sugeridos
#' continuam validadas no momento em que o engine e selecionado. Rfast e Suggests
#' porque so e usado no fallback R puro; quando o kernel C nativo esta disponivel
#' (caso padrao), Rfast nao e acionado.
#'
#' @return Vetor de nomes de pacotes requeridos pelas etapas recipes
#' @noRd
sby_required_sampling_pkgs <- function(){
  # Retorna dependencias necessarias para construir e executar as etapas padrao
  return(c(
    "sbyadanear",
    "recipes",
    "generics",
    "rlang",
    "cli",
    "tibble",
    "FNN",
    "RcppHNSW",
    "RcppParallel"
  ))
}
####
## Fim
#
