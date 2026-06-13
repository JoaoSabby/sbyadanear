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

  # Verifica se a especificacao recebida e uma formula supervisionada two-sided.
  # formula.tools::is.two.sided() cobre formula, call e expression sem depender
  # de length(formula) == 3L, que falha para formulas com environments alterados.
  if(!inherits(sby_formula, "formula") || !formula.tools::is.two.sided(sby_formula)){

    # Aborta quando a formula nao contem lado esquerdo e lado direito
    sby_adanear_abort(
      sby_message = "'sby_formula' deve ser formula no formato alvo ~ preditores"
    )
  }

  # Bloqueia matrizes esparsas antes de qualquer operacao. O pipeline
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

  # Verifica se os dados possuem nomes de colunas para resolver a formula
  if(is.null(colnames(sby_data)) || any(!nzchar(colnames(sby_data)))){

    # Aborta quando nao ha nomes suficientes para mapear formula aos dados
    sby_adanear_abort(
      sby_message = "'sby_data' deve possuir nomes de colunas"
    )
  }

  # Resolve o nome da coluna de desfecho no lado esquerdo da formula.
  # formula.tools::lhs.vars() extrai os nomes de variaveis do LHS sem
  # depender de stats::terms() nem de conversao para data.frame.
  sby_target_name <- formula.tools::lhs.vars(sby_formula)

  # Verifica se exatamente um desfecho foi informado
  if(length(sby_target_name) != 1L){

    # Aborta quando a formula nao aponta para um unico alvo
    sby_adanear_abort(
      sby_message = "'sby_formula' deve possuir exatamente uma coluna de desfecho"
    )
  }

  # Verifica se a coluna de desfecho existe nos dados
  if(!sby_target_name %in% colnames(sby_data)){

    # Aborta quando o alvo referenciado nao esta disponivel
    sby_adanear_abort(
      sby_message = paste0(
        "Coluna de desfecho nao encontrada em 'sby_data': ",
        sby_target_name
      )
    )
  }

  # Extrai o TARGET diretamente do objeto original, sem nenhuma conversao
  # de tipo. Isso preserva factor (com levels e ordered), character, integer
  # e qualquer outra classe exatamente como entrou — em especial evita que
  # as.data.frame(..., stringsAsFactors = FALSE) converta factor para integer.
  sby_target_vector <- sby_data[[sby_target_name]]

  # Resolve os nomes dos preditores no lado direito da formula.
  # formula.tools::rhs.vars(data = sby_data) expande o atalho '.' para todos
  # os nomes de colunas de sby_data, excluindo o alvo automaticamente, sem
  # precisar converter sby_data para data.frame nem chamar stats::terms().
  sby_predictor_names <- formula.tools::rhs.vars(sby_formula, data = sby_data)

  # Remove o alvo da lista de preditores por seguranca (caso o usuario
  # use '.' e o alvo esteja entre as colunas)
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

  # Restringe a interface a selecao direta de colunas existentes.
  # Transformacoes (log(x), x^2) e interacoes (x:z) geram nomes que nao
  # existem como colunas e devem ser materializados antes.
  if(!all(sby_predictor_names %in% colnames(sby_data))){

    # Aborta formulas com transformacoes ou interacoes nao materializadas
    sby_adanear_abort(
      sby_message = paste0(
        "'sby_formula' deve referenciar apenas colunas existentes em 'sby_data'. ",
        "Transformacoes, interacoes e offsets devem ser materializados em colunas ",
        "antes de chamar sbyadanear. Preditores solicitados: ",
        paste(sby_predictor_names, collapse = ", "),
        "; colunas disponiveis: ",
        paste(colnames(sby_data), collapse = ", ")
      )
    )
  }

  # Seleciona preditores via dplyr::select() — preserva os tipos originais de
  # todas as colunas (integer permanece integer, double permanece double, etc.)
  # sem nenhuma coercao implicita. A operacao e leve (selecao de colunas por
  # nome) e compativel com tibble, data.frame e tbl_df.
  sby_predictor_data <- dplyr::select(
    .data = sby_data,
    dplyr::all_of(sby_predictor_names)
  )

  # Retorna preditores e alvo resolvidos para as rotinas numericas internas.
  # Contrato de saida:
  #   sby_predictor_data  — tibble/data.frame identico ao input em tipos de colunas
  #   sby_target_vector   — vetor identico ao original: factor, character ou integer
  #   sby_predictor_names — character vector com os nomes resolvidos
  #   sby_target_name     — scalar character com o nome do desfecho
  return(list(
    sby_predictor_data  = sby_predictor_data,
    sby_target_vector   = sby_target_vector,
    sby_predictor_names = sby_predictor_names,
    sby_target_name     = sby_target_name
  ))
}
####
## Fim
#
