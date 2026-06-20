#' Validar inteiro positivo escalar
#'
#' @details
#' A funcao padroniza validacoes de hiperparametros KNN que representam contagens
#' discretas. Valores decimais sao rejeitados para evitar truncamento silencioso
#' por chamadas posteriores a `as.integer()`.
#'
#' @param sby_value Valor a validar
#'
#' @param sby_name Nome do argumento usado na mensagem de erro
#'
#' @return Inteiro positivo escalar
#'
#' @noRd
sby_validate_positive_integer_scalar <- function(sby_value, sby_name){
  # Verifica contrato numerico antes de converter o valor
  if(!(is.numeric(sby_value) && length(sby_value) == 1L && !is.na(sby_value) && is.finite(sby_value))){

    # Aborta quando a entrada nao e escalar numerico finito
    sby_adanear_abort(
      sby_message = paste0("'", sby_name, "' deve ser escalar inteiro positivo")
    )
  }

  # Rejeita valores nao inteiros para impedir truncamento silencioso
  if(sby_value != floor(sby_value) || sby_value < 1L || sby_value > .Machine$integer.max){

    # Aborta quando o valor nao representa contagem positiva valida
    sby_adanear_abort(
      sby_message = paste0("'", sby_name, "' deve ser escalar inteiro positivo")
    )
  }

  # Retorna valor em armazenamento integer para chamadas KNN posteriores
  return(as.integer(sby_value))
}
####
## Fim
#
