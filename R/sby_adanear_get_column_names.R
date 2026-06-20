#' Obter nomes de colunas dos preditores
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_predictor_data Dados preditores em data frame ou matriz
#'
#' @return Vetor de nomes de colunas validado ou gerado
#'
#' @noRd
sby_adanear_get_column_names <- function(sby_predictor_data){
  
  # Captura nomes de colunas existentes no objeto de entrada
  sby_names <- colnames(sby_predictor_data)

  # Gera nomes padronizados quando a entrada nao possui nomes de colunas
  if(is.null(sby_names)){

    # Cria nomes sequenciais para todas as colunas preditoras
    sby_names <- paste0(
      "V",
      seq_len(collapse::fncol(sby_predictor_data))
    )
  }

  # Retorna nomes de colunas disponiveis para uso posterior
  return(sby_names)
}
####
## Fim
#
