#' Adicionar uma etapa recipes de balanceamento NearMiss-1
#'
#' @description
#' `sby_step_nearmiss()` adiciona a uma `recipe` uma etapa supervisionada de
#' subamostragem NearMiss-1 para problemas de classificação binária. A etapa é
#' ajustada durante `prep()` a partir de uma única coluna de desfecho selecionada
#' e, por padrão, é ignorada em novos dados durante `bake()` porque modifica o
#' número de linhas do conjunto processado.
#'
#' @details
#' A função utiliza uma arquitetura KNN configurável sobre preditores numéricos
#' padronizados por Z-score. Três escolhas trabalham em conjunto:
#' `sby_knn_algorithm` define a estrutura/estratégia de busca,
#' `sby_knn_engine` define a implementação computacional e
#' `sby_knn_distance_metric` define a noção de proximidade.
#'
#' Em geral, não: para o uso cotidiano, mantenha `sby_knn_engine = "auto"` e
#' `sby_knn_algorithm = "auto"`. Nesse modo, com métrica euclidiana, o pacote
#' prefere a engine `native` exata quando as rotinas nativas estão carregadas e
#' usa `FNN` como fallback exato. Busca aproximada por `RcppHNSW` só é escolhida
#' automaticamente quando a opção `sbyadanear.sby_knn_allow_approx = TRUE` está
#' ativa ou quando o engine é escolhido explicitamente.
#'
#' ## Engines disponíveis
#'
#' | Engine | Melhor para | Vantagens | Limitações e cuidados |
#' |---|---|---|---|
#' | `"native"` | Busca exata euclidiana densa por kernel C/C++ interno. | Retorna índices 1-based e distâncias euclidianas reais; honra controle explícito de self-neighbor nas rotas internas. | Somente `"euclidean"`; só combina com `"auto"` ou `"brute"`. |
#' | `"FNN"` | Bases pequenas a médias, distância euclidiana e busca exata. | Usa `FNN::get.knnx()`; a rota `"brute"` pode acionar o kernel nativo de compatibilidade quando disponível. | Aceita apenas `"euclidean"` neste pacote; só combina com `"auto"`, `"kd_tree"`, `"cover_tree"` ou `"brute"`. |
#' | `"RcppHNSW"` | Bases grandes, alta dimensionalidade, métricas `"cosine"`/`"ip"` e consultas em que velocidade é mais importante que exatidão perfeita. | Implementa HNSW de alto desempenho, usa `sby_knn_workers`, costuma escalar melhor que busca exata e suporta `"euclidean"`, `"cosine"` e `"ip"`. | A busca é aproximada; exige calibrar `sby_knn_hnsw_m` e `sby_knn_hnsw_ef`; consome memória para o grafo; resultados podem diferir de uma busca exata quando `ef` é baixo. |
#' | `"KernelKnn"` | Comparação exata euclidiana via OpenMP externo. | Permite benchmark de `KernelKnn::knn.index.dist()` dentro do contrato comum `nn.index`/`nn.dist`. | Somente `"euclidean"`; usa OpenMP e pode competir com outros backends multithread. |
#' | `"bigKNN"` | Busca exata euclidiana baseada em `bigmemory::big.matrix`. | Usa `bigKNN::knn_bigmatrix()` e permite avaliar blocos em bases grandes. | Somente `"euclidean"`; converte a referência densa para `big.matrix`. |
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
#' Na rota exata `FNN` com `sby_knn_algorithm = "brute"`, quando os
#' kernels nativos estão disponíveis, o pacote usa produtos matriciais BLAS para
#' calcular top-k sem materializar uma matriz completa de distâncias. As chamadas
#' internas retornam apenas índices ou apenas distâncias quando a etapa precisa
#' de um único componente, reduzindo alocações em ADASYN e NearMiss.
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
#' Os argumentos `sby_knn_hnsw_m` e `sby_knn_hnsw_ef` só afetam a rota
#' `sby_knn_engine = "RcppHNSW"`. O parâmetro `sby_knn_hnsw_m`
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
#' - Para auditoria, bases pequenas ou necessidade de vizinhos exatos, prefira `sby_knn_engine = "FNN"` com `sby_knn_algorithm = "kd_tree"` ou `"brute"`.
#' - Para bases grandes, embeddings, `"ip"` ou alta dimensionalidade, use
#'   `sby_knn_engine = "RcppHNSW"` e ajuste `sby_knn_hnsw_m`/`sby_knn_hnsw_ef`.
#' - Em ADASYN, vizinhos aproximados podem mudar quais regiões recebem amostras
#'   sintéticas; em NearMiss, podem mudar quais exemplos majoritários são retidos.
#'   Aumente `sby_knn_hnsw_ef` quando essa estabilidade for importante.
#'
#' @param recipe Objeto `recipe` que receberá a etapa de balanceamento.
#'
#' @param ... Seletores `recipes` ou `tidyselect` usados para identificar exatamente uma coluna de desfecho binário.
#'
#' @param role Papel armazenado na etapa para compatibilidade com `recipes`.
#'
#' @param trained Indicador lógico escalar que informa se a etapa já passou por `prep()`.
#'
#' @param columns Vetor de caracteres ou `NULL` com o nome da coluna de desfecho resolvida durante `prep()`.
#'
#' @param sby_nearmiss_ratio Valor numérico escalar maior que zero que controla a quantidade de registros majoritários retidos pelo NearMiss-1 em relação ao tamanho da classe rara. O alvo é `floor(n_minoria * sby_nearmiss_ratio)`, limitado à maioria disponível.
#'
#' @param sby_knn_under_k Número inteiro positivo de vizinhos usados pela etapa NearMiss-1.
#'
#' @param sby_seed Valor numérico inteiro usado para inicializar o gerador pseudoaleatório.
#'
#' @param sby_audit Indicador lógico escalar que controla se metadados de auditoria devem ser preservados.
#'
#' @param sby_restore_types Indicador lógico escalar que define se tipos numéricos originais devem ser restaurados.
#'
#' @param sby_knn_algorithm String escalar que escolhe a estratégia de busca KNN: `"auto"`, `"kd_tree"`, `"cover_tree"` ou `"brute"`. Use `"auto"` para deixar o pacote escolher uma opção compatível com o engine e a dimensionalidade; informe uma alternativa explícita quando quiser controlar o compromisso entre exatidão, velocidade de execução, consumo de memória e suporte a métricas. Consulte os detalhes para recomendações por algoritmo.
#'
#' @param sby_knn_engine String escalar que escolhe a biblioteca usada para executar a busca KNN: `"auto"`, `"native"`, `"FNN"`, `"RcppHNSW"`, `"KernelKnn"` ou `"bigKNN"`. Na maioria dos casos, mantenha `"auto"`; informe explicitamente apenas quando precisar de uma implementação específica, engine nativa exata, compatibilidade `FNN` ou busca aproximada HNSW por `RcppHNSW`. Consulte os detalhes para saber quando o engine precisa ser declarado.
#'
#' @param sby_knn_distance_metric String escalar que define a geometria da vizinhança: `"euclidean"`, `"cosine"` ou `"ip"`. A escolha muda o significado de proximidade e também restringe engines e algoritmos disponíveis; `"euclidean"` é a opção mais geral, `"cosine"` privilegia direção angular e `"ip"` usa produto interno via `RcppHNSW`. Consulte os detalhes para recomendações.
#'
#' @param sby_knn_parallel_backend Backend de paralelismo KNN. Use `"parallel"` para o particionamento por blocos com o pacote base `parallel` ou `"RcppParallel"` para acionar threads nativos no kernel bruto exato (`sby_knn_engine = "native"` ou compatibilidade `"FNN"` + `"brute"`).
#'
#' @param sby_knn_workers Número de workers KNN configurado.
#'
#' @param sby_knn_hnsw_m Número inteiro positivo usado apenas quando o engine efetivo é `"RcppHNSW"`. Controla a conectividade máxima do grafo (`M`): valores maiores aumentam a chance de recuperar vizinhos melhores e tornam o índice mais robusto, mas consomem mais memória e tempo de construção. O padrão `16L` costuma ser um bom equilíbrio; aumente em bases grandes, ruidosas ou de alta dimensionalidade quando recall for mais importante que memória.
#'
#' @param sby_knn_query_chunk_size Número inteiro positivo que define quantas linhas de consulta KNN são processadas por bloco. O padrão é `1000L`. Valores maiores reduzem overhead de chamadas e podem favorecer kernels BLAS/MKL em matrizes densas, enquanto valores menores reduzem pico de memória em bases muito grandes.
#'
#' @param sby_knn_hnsw_ef Número inteiro positivo usado apenas quando o engine efetivo é `"RcppHNSW"`. Controla a largura da lista dinâmica de candidatos (`ef`) durante construção/consulta: valores maiores aproximam a busca do resultado exato e estabilizam ADASYN/NearMiss, mas deixam as consultas mais lentas. O padrão `200L` prioriza qualidade; reduza para velocidade ou aumente quando a vizinhança aproximada precisar de mais fidelidade.
#'
#' @param skip Indicador lógico escalar que define se a etapa deve ser ignorada em novos dados.
#'
#' @param id Identificador recipes da etapa.
#'
#' @return Objeto `recipe` com uma etapa `sby_step_nearmiss` adicionada ao pipeline.
#'
#' @export
sby_step_nearmiss <- function(
  recipe,
  ...,
  role = NA,
  trained = FALSE,
  columns = NULL,
  sby_nearmiss_ratio = 1,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE,
  sby_restore_types = TRUE,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "native", "FNN", "RcppHNSW", "KernelKnn", "bigKNN"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L,
  sby_knn_parallel_backend = c("parallel", "RcppParallel"),
  sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L,
  sby_knn_query_chunk_size = 1000L,
  skip = TRUE,
  id = recipes::rand_id("nearmiss")
){
  
  # Verifica se ha solicitacao de interrupcao antes de configurar a etapa
  sby_adanear_check_user_interrupt()

  # Valida dependencias declaradas para a etapa recipes
  recipes::recipes_pkg_check(
    required_pkgs.step_sby_step_nearmiss()
  )

  # Captura seletores de desfecho informados pelo chamador
  sby_terms <- rlang::enquos(...)

  # Valida indicadores logicos de controle da etapa
  sby_audit <- sby_validate_logical_scalar(sby_value = sby_audit, sby_name = "sby_audit")
  sby_restore_types <- sby_validate_logical_scalar(sby_value = sby_restore_types, sby_name = "sby_restore_types")
  skip <- sby_validate_logical_scalar(sby_value = skip, sby_name = "skip")

  # Valida a semente ainda na construcao da etapa para evitar falhas tardias
  sby_seed <- sby_validate_seed(
    sby_seed = sby_seed
  )

  # Valida hiperparametros KNN discretos ainda na construcao da etapa
  sby_knn_under_k <- sby_validate_positive_integer_scalar(
    sby_value = sby_knn_under_k,
    sby_name  = "sby_knn_under_k"
  )

  # Resolve opcoes declaradas de algoritmo, engine e metrica KNN
  sby_knn_algorithm <- match.arg(arg = sby_knn_algorithm)
  sby_knn_engine <- match.arg(arg = sby_knn_engine)
  sby_knn_distance_metric <- match.arg(arg = sby_knn_distance_metric)

  # Valida recursos paralelos e parametros HNSW
  sby_knn_workers <- sby_validate_knn_workers(sby_knn_workers = sby_knn_workers)
  sby_knn_parallel_backend <- sby_validate_knn_parallel_backend(
    sby_knn_parallel_backend = sby_knn_parallel_backend
  )
  sby_hnsw_params <- sby_validate_hnsw_params(
    sby_knn_hnsw_m = sby_knn_hnsw_m,
    sby_knn_hnsw_ef = sby_knn_hnsw_ef
  )
  sby_knn_query_chunk_size <- sby_validate_knn_query_chunk_size(
    sby_knn_query_chunk_size = sby_knn_query_chunk_size
  )

  # Adiciona etapa configurada ao objeto recipe
  return(recipes::add_step(
    rec = recipe,
    object = sby_step_sampling_new(
      sby_subclass                = "sby_step_nearmiss",
      sby_sampling_method         = "nearmiss",
      sby_terms                   = sby_terms,
      sby_role                    = role,
      sby_trained                 = trained,
      sby_columns                 = columns,
      sby_adasyn_ratio              = NA_real_,
      sby_nearmiss_ratio             = sby_nearmiss_ratio,
      sby_knn_over_k              = NA_integer_,
      sby_knn_under_k             = sby_knn_under_k,
      sby_seed                    = sby_seed,
      sby_audit                   = sby_audit,
      sby_restore_types           = sby_restore_types,
      sby_knn_algorithm           = sby_knn_algorithm,
      sby_knn_engine              = sby_knn_engine,
      sby_knn_distance_metric     = sby_knn_distance_metric,
      sby_knn_workers             = sby_knn_workers,
      sby_knn_parallel_backend     = sby_knn_parallel_backend,
      sby_knn_hnsw_m              = sby_hnsw_params$sby_knn_hnsw_m,
      sby_knn_hnsw_ef             = sby_hnsw_params$sby_knn_hnsw_ef,
      sby_knn_query_chunk_size = sby_knn_query_chunk_size,
      sby_skip                    = skip,
      sby_id                      = id
    )
  ))
}
####
## Fim
#
