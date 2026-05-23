#' Extrair preditores e alvo de formula e dados
#'
#' @details
#' A funcao implementa uma unidade interna do fluxo de balanceamento com contrato de entrada explicito e retorno controlado
#' A documentacao descreve a intencao operacional para apoiar manutencao, auditoria e revisao tecnica do pacote
#'
#' @param sby_formula Formula no formato alvo ~ preditores
#' @param sby_data Dados tabulares contendo alvo e preditores
#'
#' @return Lista com preditores, alvo e nomes resolvidos
#' @noRd
sby_extract_formula_data <- function(sby_formula, sby_data){
  
  # Verifica se a especificacao recebida e uma formula supervisionada
  if(!(inherits(sby_formula, "formula") && length(sby_formula) == 3L)){
    
    # Aborta quando a formula nao contem lado esquerdo e lado direito
    sby_adanear_abort(
      sby_message = "'sby_formula' deve ser formula no formato alvo ~ preditores"
    )
  }
  
  # Bloqueia matrizes esparsas antes da conversao para data frame. O pipeline
  # numerico atual e baseado em matrizes densas para z-score, KNN e rotinas C.
  if(sby_is_sparse_matrix(
    sby_x = sby_data
  )){

    # Aborta cedo para evitar densificacao acidental de objetos Matrix grandes
    sby_adanear_abort(
      sby_message = paste0(
        "'sby_data' recebeu uma matriz esparsa do pacote Matrix, mas as ",
        "rotinas atuais usam matrix double densa. Forneca data.frame/matrix ",
        "densa ou reduza/materialize os preditores antes do balanceamento."
      )
    )
  }

  # Verifica se os dados possuem classe tabular suportada
  if(!(is.data.frame(sby_data) || is.matrix(sby_data))){
    
    # Aborta quando os dados nao sao data frame nem matriz
    sby_adanear_abort(
      sby_message = "'sby_data' deve ser tibble, data.frame ou matrix"
    )
  }
  
  # Converte matriz para data frame para selecao por nomes de colunas
  sby_data_frame <- as.data.frame(
    x = sby_data,
    stringsAsFactors = FALSE
  )
  
  # Verifica se os dados possuem nomes de colunas para resolver a formula
  if(is.null(names(sby_data_frame)) || any(!nzchar(names(sby_data_frame)))){
    
    # Aborta quando nao ha nomes suficientes para mapear formula aos dados
    sby_adanear_abort(
      sby_message = "'sby_data' deve possuir nomes de colunas"
    )
  }
  
  # Resolve o nome da coluna de desfecho no lado esquerdo da formula
  sby_target_name <- all.vars(
    expr = sby_formula[[2L]]
  )
  
  # Verifica se exatamente um desfecho foi informado
  if(length(sby_target_name) != 1L){
    
    # Aborta quando a formula nao aponta para um unico alvo
    sby_adanear_abort(
      sby_message = "'sby_formula' deve possuir exatamente uma coluna de desfecho"
    )
  }
  
  # Verifica se a coluna de desfecho existe nos dados
  if(!sby_target_name %in% names(sby_data_frame)){
    
    # Aborta quando o alvo referenciado nao esta disponivel
    sby_adanear_abort(
      sby_message = paste0(
        "Coluna de desfecho nao encontrada em 'sby_data': ",
        sby_target_name
      )
    )
  }
  
  # Expande termos do lado direito, incluindo atalhos como ponto
  sby_terms <- stats::terms(
    x = sby_formula,
    data = sby_data_frame
  )
  sby_predictor_names <- attr(
    x = sby_terms,
    which = "term.labels"
  )
  
  # Remove o alvo da lista de preditores por seguranca
  sby_predictor_names <- setdiff(
    x = sby_predictor_names,
    y = sby_target_name
  )
  
  # Verifica se a formula selecionou ao menos um preditor
  if(length(sby_predictor_names) < 1L){
    
    # Aborta quando nao ha colunas preditoras para balancear
    sby_adanear_abort(
      sby_message = "'sby_formula' deve selecionar ao menos uma coluna preditora"
    )
  }
  
  # Restringe a interface a selecao direta de colunas existentes
  if(!all(sby_predictor_names %in% colnames(sby_data_frame))){
    
    predictors_data_frame <- colnames(sby_data_frame)
    
    # Aborta formulas com transformacoes ou interacoes nao materializadas como colunas
    sby_adanear_abort(
      sby_message = paste0(
        "'sby_formula' deve referenciar apenas colunas existentes em 'sby_data'. ",
        "Transformacoes, interacoes e offsets devem ser materializados em colunas ",
        "antes de chamar sbyadanear. Preditores solicitados: ",
        paste(sby_predictor_names, collapse = ", "),
        "; colunas disponiveis: ",
        paste(predictors_data_frame, collapse = ", ")
      )
    )
  }
  
  # Retorna preditores e alvo resolvidos para as rotinas numericas internas
  return(list(
    sby_predictor_data  = sby_data_frame[, sby_predictor_names, drop = FALSE],
    sby_target_vector   = sby_data_frame[[sby_target_name]],
    sby_predictor_names = sby_predictor_names,
    sby_target_name     = sby_target_name
  ))
}
####
## Fim
#
