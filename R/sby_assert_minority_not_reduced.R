#' Garantir que a classe rara original nao foi reduzida
#'
#' @details
#' Esta pos-condicao centraliza a regra de seguranca das rotas ADASYN: todas as
#' observacoes da classe rara de origem devem permanecer no resultado. O rotulo
#' raro e determinado exclusivamente na entrada, mesmo que a sobreamostragem
#' altere qual classe seria considerada minoritaria na saida.
#'
#' @param sby_input_target Vetor alvo binario antes do ADASYN.
#' @param sby_output_target Vetor alvo depois do ADASYN ou ADANEAR.
#' @param sby_context Nome da rotina usado na mensagem de erro.
#' @param sby_minority_label Rotulo da classe rara identificado na entrada. Se
#' `NULL`, o rotulo e inferido apenas para compatibilidade com chamadas internas.
#' @param sby_input_count Contagem rara original ja calculada. Se `NULL`, a
#' contagem e obtida de `sby_input_target`.
#'
#' @return Invisivelmente, a quantidade final da classe rara original.
#'
#' @noRd
sby_assert_minority_not_reduced <- function(sby_input_target, sby_output_target,
                                            sby_context = "ADASYN",
                                            sby_minority_label = NULL,
                                            sby_input_count = NULL){
  sby_input_factor <- if(is.factor(sby_input_target)){
    sby_input_target
  }else{
    factor(sby_input_target, levels = unique(as.character(sby_input_target)))
  }
  if(is.null(sby_minority_label)){
    sby_input_roles <- sby_binary_class_counts_fast(sby_input_factor)
    sby_minority_label <- sby_input_roles$sby_minority_label
  }
  sby_minority_label <- as.character(sby_minority_label)
  if(is.null(sby_input_count)){
    sby_input_count <- sum(as.character(sby_input_factor) == sby_minority_label)
  }
  sby_output_count <- sum(as.character(sby_output_target) == sby_minority_label)

  if(sby_output_count < sby_input_count){
    sby_adanear_abort(paste0(
      "Falha de seguranca em ", sby_context,
      ": a classe rara original nunca pode ser reduzida (entrada = ",
      sby_input_count, ", saida = ", sby_output_count, ")."
    ))
  }

  invisible(sby_output_count)
}
####
## Fim
#
