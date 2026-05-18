#' Executar balanceamento combinado ADASYN e NearMiss-1
#'
#' @description
#' `sby_adanear()` executa um pipeline híbrido de balanceamento para problemas
#' binários: primeiro aplica ADASYN para expandir adaptativamente a classe
#' minoritária e, em seguida, aplica NearMiss-1 para reduzir a classe majoritária.
#' O objetivo é combinar geração sintética em regiões difíceis com retenção
#' seletiva de exemplos majoritários próximos à fronteira de decisão.
#'
#' @usage
#' sby_adanear(
#'   sby_formula,
#'   sby_data,
#'   sby_over_ratio = 0.2,
#'   sby_under_ratio = 0.5,
#'   sby_knn_over_k = 5L,
#'   sby_knn_under_k = 5L,
#'   sby_seed = sample.int(10L^5L, 1L),
#'   sby_audit = FALSE,
#'   sby_restore_types = TRUE,
#'   sby_knn_algorithm = c(
#'     "auto", "kd_tree", "cover_tree", "brute"
#'   ),
#'   sby_knn_engine = c(
#'     "auto", "FNN", "RcppHNSW"
#'   ),
#'   sby_knn_distance_metric = c(
#'     "euclidean", "ip", "cosine"
#'   ),
#'   sby_knn_workers = 1L,
#'   sby_knn_hnsw_m = 16L,
#'   sby_knn_hnsw_ef = 200L
#' )
#'
#' @param sby_formula Fórmula no formato `alvo ~ preditores` usada para identificar uma única coluna de desfecho binário e as colunas preditoras numéricas em `sby_data`. O lado direito deve referenciar apenas colunas ja existentes; transformacoes, interacoes e offsets precisam ser materializados antes da chamada. Não possui valor padrão; use `alvo ~ .` para selecionar todos os demais campos como preditores.
#' @param sby_data Data frame, tibble ou matriz com a coluna de desfecho e as variáveis preditoras numéricas referenciadas em `sby_formula`. Não possui valor padrão. Esse objeto define o espaço comum no qual serão calculadas tanto a dificuldade adaptativa do ADASYN quanto a proximidade NearMiss-1.
#' @param sby_over_ratio Valor numérico escalar que controla a intensidade da sobreamostragem ADASYN antes da etapa NearMiss. O padrão é `0.2`; em bases pequenas, qualquer valor positivo gera ao menos uma linha sintética antes da subamostragem. Valores maiores podem melhorar cobertura da minoria, mas também propagam ruído em regiões ambíguas.
#' @param sby_under_ratio Valor numérico escalar no intervalo `(0, 1]` que controla a razão mínima desejada entre minoria e maioria após o NearMiss-1. O padrão `0.5` permite reter até duas vezes a quantidade minoritária; `1` reduz a maioria até igualar a minoria quando houver registros suficientes.
#' @param sby_knn_over_k Número inteiro positivo de vizinhos usados pela etapa ADASYN para estimar dificuldade local e gerar amostras sintéticas. O padrão é `5L`. Essa escolha influencia onde a expansão minoritária será concentrada e o grau de suavização da avaliação local.
#' @param sby_knn_under_k Número inteiro positivo de vizinhos minoritários usados pela etapa NearMiss-1 para ranquear observações majoritárias. O padrão é `5L`. Essa escolha controla a sensibilidade da retenção majoritária à fronteira local criada após a sobreamostragem.
#' @param sby_seed Valor numérico inteiro utilizado para inicializar o gerador de números pseudoaleatórios nas duas etapas do pipeline. O padrão é `sample.int(10L^5L, 1L)`, gerando uma semente inteira aleatória quando o usuário não informa valor. Informe uma semente fixa para garantir reprodutibilidade conjunta das amostras sintéticas, desempates, índices retidos e diagnósticos finais.
#' @param sby_audit Indicador lógico escalar que controla se o retorno deve incluir resultados intermediários das etapas ADASYN e NearMiss. O padrão é `FALSE`, retornando somente os dados balanceados finais. Quando `TRUE`, o retorno passa a ser uma lista com auditoria completa, útil para inspeção de contagens e validação experimental.
#' @param sby_restore_types Indicador lógico escalar que define se tipos numéricos originais devem ser restaurados ao final do pipeline. O padrão é `TRUE`. Mantê-lo ativado facilita reintegração com bases tabulares originais; desativá-lo mantém a saída no formato numérico produzido pelas transformações internas.
#' @param sby_knn_algorithm String escalar que escolhe a estratégia de busca KNN: `"auto"`, `"kd_tree"`, `"cover_tree"` ou `"brute"`. Use `"auto"` para deixar o pacote escolher uma opção compatível com o engine e a dimensionalidade; informe uma alternativa explícita quando quiser controlar o compromisso entre exatidão, velocidade de execução, consumo de memória e suporte a métricas. Consulte os detalhes para recomendações por algoritmo.
#' @param sby_knn_engine String escalar que escolhe a biblioteca usada para executar a busca KNN: `"auto"`, `"FNN"` ou `"RcppHNSW"`. Na maioria dos casos, mantenha `"auto"`; informe explicitamente apenas quando precisar de uma implementação específica, paralelismo exato via `FNN` ou busca aproximada HNSW por `RcppHNSW`. Consulte os detalhes para saber quando o engine precisa ser declarado.
#' @param sby_knn_distance_metric String escalar que define a geometria da vizinhança: `"euclidean"`, `"cosine"` ou `"ip"`. A escolha muda o significado de proximidade e também restringe engines e algoritmos disponíveis; `"euclidean"` é a opção mais geral, `"cosine"` privilegia direção angular e `"ip"` usa produto interno via `RcppHNSW`. Consulte os detalhes para recomendações.
#' @param sby_knn_workers Número inteiro positivo de workers usados nas consultas KNN quando o engine efetivo oferece suporte. O padrão é `1L`. Aumentar esse número pode acelerar o pipeline completo, mas eleva uso de recursos e deve ser compatível com o backend escolhido.
#' @param sby_knn_hnsw_m Número inteiro positivo usado apenas quando o engine efetivo é `"RcppHNSW"`. Controla a conectividade máxima do grafo (`M`): valores maiores aumentam a chance de recuperar vizinhos melhores e tornam o índice mais robusto, mas consomem mais memória e tempo de construção. O padrão `16L` costuma ser um bom equilíbrio; aumente em bases grandes, ruidosas ou de alta dimensionalidade quando recall for mais importante que memória.
#' @param sby_knn_hnsw_ef Número inteiro positivo usado apenas quando o engine efetivo é `"RcppHNSW"`. Controla a largura da lista dinâmica de candidatos (`ef`) durante construção/consulta: valores maiores aproximam a busca do resultado exato e estabilizam ADASYN/NearMiss, mas deixam as consultas mais lentas. O padrão `200L` prioriza qualidade; reduza para velocidade ou aumente quando a vizinhança aproximada precisar de mais fidelidade.
#'
#' @details
#' A função utiliza uma arquitetura KNN configurável sobre preditores numéricos
#' padronizados por Z-score. Três escolhas trabalham em conjunto:
#' `sby_knn_algorithm` define a estrutura/estratégia de busca,
#' `sby_knn_engine` define a implementação computacional e
#' `sby_knn_distance_metric` define a noção de proximidade.
#'
#' ## Preciso informar `sby_knn_engine`?
#'
#' Em geral, não: para o uso cotidiano, mantenha `sby_knn_engine = "auto"` e
#' `sby_knn_algorithm = "auto"`. Nesse modo, o pacote usa `FNN` com paralelismo exato por blocos quando `sby_knn_workers > 1L`; o algoritmo também
#' é escolhido automaticamente pela dimensionalidade dos preditores. Informe o
#' engine explicitamente quando quiser garantir uma biblioteca específica, usar
#' algoritmo que pertence a uma família específica. Na API atual, informar apenas
#' `sby_knn_engine = "RcppHNSW"` quando quiser a implementação HNSW do pacote
#' `RcppHNSW`.
#'
#' ## Engines disponíveis
#'
#' | Engine | Melhor para | Vantagens | Limitações e cuidados |
#' |---|---|---|---|
#' | `"FNN"` | Bases pequenas a médias, distância euclidiana, execução sequencial ou paralela por blocos e busca exata. | Simples, rápido para baixa/média dimensionalidade, determinístico e sem aproximação. | Aceita apenas `"euclidean"` neste pacote; só combina com `"auto"`, `"kd_tree"`, `"cover_tree"` ou `"brute"`. |
#' | `"RcppHNSW"` | Bases grandes, alta dimensionalidade, métricas `"cosine"`/`"ip"` e consultas em que velocidade é mais importante que exatidão perfeita. | Implementa HNSW de alto desempenho, usa `sby_knn_workers`, costuma escalar melhor que busca exata e suporta `"euclidean"`, `"cosine"` e `"ip"`. | A busca é aproximada; exige calibrar `sby_knn_hnsw_m` e `sby_knn_hnsw_ef`; consome memória para o grafo; resultados podem diferir de uma busca exata quando `ef` é baixo. |
#'
#' ## Algoritmos disponíveis
#'
#' | Algoritmo | Engine compatível | Tipo | Quando usar | Evite quando |
#' |---|---|---|---|---|
#' | `"kd_tree"` | `FNN` | Exato | Dados euclidianos com poucas ou médias dimensões; tende a ser eficiente quando as partições espaciais ainda discriminam bem os vizinhos. | Alta dimensionalidade, muitas variáveis ruidosas ou métrica não euclidiana. |
#' | `"cover_tree"` | `FNN` | Exato | Alternativa exata para dados euclidianos quando a estrutura intrínseca pode favorecer árvore de cobertura. | Quando testes rápidos mostram desempenho inferior a `"kd_tree"`/`"brute"`; não serve para cosseno ou produto interno. |
#' | `"brute"` | `FNN` | Exato | Alta dimensionalidade moderada, bases pequenas/médias, auditorias ou cenários em que simplicidade e previsibilidade importam mais que indexação. | Bases muito grandes, pois compara muitos pares e pode ficar lento. |
#'
#' ## Desempenho, memória e matrizes esparsas
#'
#' A velocidade deve ser interpretada como um atributo de qualidade da
#' configuração KNN, junto com fidelidade da vizinhança, consumo de memória e
#' compatibilidade de métrica. Em termos práticos:
#'
#' | Escolha | Velocidade esperada | Memória | Por que isso acontece |
#' |---|---|---|---|
#' | Paralelismo (`sby_knn_workers > 1L`) | Pode reduzir tempo de consulta em matrizes grandes. | Aumenta uso simultâneo de CPU e pode elevar pressão de memória. | O trabalho é dividido entre workers em `FNN` por blocos de consulta exatos e em `RcppHNSW` pelos threads nativos. |
#'
#' As rotinas atuais operam sobre `matrix` numérica densa após seleção de
#' preditores e padronização. Matrizes esparsas do pacote `Matrix` são rejeitadas
#' antes de densificação implícita para evitar estouro de memória. Portanto, em
#' dados muito esparsos, materialize conscientemente uma matriz densa somente se
#' houver memória suficiente, reduza dimensionalidade/seleção de variáveis antes
#' do balanceamento, ou use outro pré-processamento que produza preditores densos.
#' Em matrizes densas muito largas, prefira testar uma busca aproximada ou uma
#' busca exaustiva paralelizável, pois árvores exatas tendem a perder vantagem.
#'
#' ## Métricas de distância
#'
#' | Métrica | Interpretação | Compatibilidade | Recomendação |
#' |---|---|---|---|
#' | `"cosine"` | Distância angular; compara a orientação dos vetores e reduz a influência da norma após normalização L2. | `RcppHNSW`; não é aceita por `FNN`. | Use quando o padrão relativo entre variáveis importa mais que o tamanho absoluto, como composições, assinaturas de perfil e vetores de alta dimensionalidade já materializados como matriz densa. |
#' | `"ip"` | Produto interno convertido em distância; após normalização L2, fica próximo de uma comparação por similaridade angular. | Somente `RcppHNSW` neste pacote. | Use quando o modelo conceitual é similaridade por produto interno ou quando você precisa alinhar a busca a embeddings/vetores normalizados; requer busca aproximada. |
#'
#' ## Parâmetros exclusivos da rota HNSW
#'
#' Os argumentos `sby_knn_hnsw_m` e `sby_knn_hnsw_ef` só afetam a rota
#' `sby_knn_engine = "RcppHNSW"`; eles não configuram o `"Hnsw"` de
#' representa a conectividade máxima do grafo: valores maiores criam mais arestas,
#' melhoram o recall e tornam a busca mais robusta, mas aumentam memória e tempo
#' de construção. O padrão `16L` é conservador; valores como 24 ou 32 podem ser
#' úteis em bases grandes, ruidosas ou de alta dimensionalidade. `sby_knn_hnsw_ef`
#' representa quantos candidatos são explorados dinamicamente na busca: deve ser
#' pelo menos tão grande quanto o número de vizinhos solicitado e, internamente, é
#' limitado ao número de linhas da base de referência. O padrão `200L` favorece
#' qualidade; reduza para acelerar quando pequenas perdas de recall forem
#' aceitáveis, ou aumente quando NearMiss/ADASYN ficarem sensíveis a vizinhos
#' aproximados subótimos.
#'
#' ## Recomendações práticas
#'
#' - Comece com `sby_knn_engine = "auto"`, `sby_knn_algorithm = "auto"` e
#'   `sby_knn_distance_metric = "euclidean"`.
#' - Para auditoria, bases pequenas ou necessidade de vizinhos exatos, prefira
#' - Para bases grandes, embeddings, `"ip"` ou alta dimensionalidade, use
#'   `sby_knn_engine = "RcppHNSW"` e ajuste `sby_knn_hnsw_m`/`sby_knn_hnsw_ef`.
#' - Em ADASYN, vizinhos aproximados podem mudar quais regiões recebem amostras
#'   sintéticas; em NearMiss, podem mudar quais exemplos majoritários são retidos.
#'   Aumente `sby_knn_hnsw_ef` quando essa estabilidade for importante.
#'
#' @references
#' He, H., Bai, Y., Garcia, E. A., & Li, S. (2008). ADASYN: Adaptive synthetic
#' sampling approach for imbalanced learning. In *2008 IEEE International Joint
#' Conference on Neural Networks* (pp. 1322-1328). IEEE.
#'
#' Malkov, Y. A., & Yashunin, D. A. (2018). Efficient and robust approximate
#' nearest neighbor search using Hierarchical Navigable Small World graphs.
#' *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 42(4),
#' 824-836.
#'
#' @return Tibble balanceado quando `sby_audit = FALSE`; lista com resultados de sobreamostragem, subamostragem, dados finais e diagnósticos quando `sby_audit = TRUE`.
#' @export
sby_adanear <- function(
  sby_formula,
  sby_data,
  sby_over_ratio = 0.2,
  sby_under_ratio = 0.5,
  sby_knn_over_k = 5L,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_restore_types = TRUE,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "FNN", "RcppHNSW"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L,
  sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L
){
  sby_adanear_check_user_interrupt()

  sby_formula_data <- sby_extract_formula_data(sby_formula = sby_formula, sby_data = sby_data)
  sby_predictor_data <- sby_formula_data$sby_predictor_data
  sby_target_vector <- sby_formula_data$sby_target_vector

  sby_audit <- sby_validate_logical_scalar(sby_audit, "sby_audit")
  sby_restore_types <- sby_validate_logical_scalar(sby_restore_types, "sby_restore_types")

  sby_validate_sampling_inputs(sby_predictor_data, sby_target_vector, sby_seed)
  sby_x_matrix <- sby_adanear_as_numeric_matrix(sby_predictor_data)
  colnames(sby_x_matrix) <- sby_adanear_get_column_names(sby_predictor_data)
  sby_target_factor <- as.factor(sby_target_vector)
  sby_type_info <- sby_infer_numeric_column_types(sby_predictor_data)

  sby_matrix_result <- sby_adanear_matrix(
    sby_x_matrix = sby_x_matrix,
    sby_y_vector = sby_target_factor,
    sby_over_ratio = sby_over_ratio,
    sby_under_ratio = sby_under_ratio,
    sby_knn_over_k = sby_knn_over_k,
    sby_knn_under_k = sby_knn_under_k,
    sby_seed = sby_seed,
    sby_audit = sby_audit,
    sby_return_scaled = sby_audit,
    sby_return_original_scale = TRUE,
    sby_knn_algorithm = sby_knn_algorithm,
    sby_knn_engine = sby_knn_engine,
    sby_knn_distance_metric = sby_knn_distance_metric,
    sby_knn_workers = sby_knn_workers,
    sby_knn_hnsw_m = sby_knn_hnsw_m,
    sby_knn_hnsw_ef = sby_knn_hnsw_ef
  )

  sby_final_predictors <- if(isTRUE(sby_restore_types)){
    sby_restore_numeric_column_types(sby_matrix_result$sby_x_matrix, sby_type_info, TRUE)
  }else{
    sby_out <- as.data.frame(sby_matrix_result$sby_x_matrix, stringsAsFactors = FALSE)
    names(sby_out) <- sby_type_info$sby_column_name
    sby_out
  }
  sby_balanced_data <- sby_build_balanced_tibble(sby_final_predictors, sby_matrix_result$sby_y_vector)

  if(isTRUE(sby_audit)){
    return(list(
      sby_oversampling_result = sby_matrix_result$sby_oversampling_result,
      sby_undersampling_result = sby_matrix_result$sby_undersampling_result,
      sby_balanced_data = sby_balanced_data,
      sby_type_info = sby_type_info,
      sby_scaling_info = sby_matrix_result$sby_scaling_info,
      sby_diagnostics = sby_matrix_result$sby_diagnostics,
      sby_balanced_scaled = sby_matrix_result$sby_balanced_scaled
    ))
  }
  return(sby_balanced_data)
}
####
## Fim
#
