# sbyadanear 0.4.0 (em desenvolvimento)

## Atalho HPC (oneAPI, AVX-512)

* Novas funcoes exportadas `sby_adanear_hpc()`, `sby_adasyn_hpc()` e
  `sby_nearmiss_hpc()`. Cada uma e um atalho de alto desempenho que executa o
  fluxo estritamente no espaco padronizado, eliminando a dupla normalizacao, e
  monta o tibble final por zero-copy diretamente em C++ via `Rcpp::List`.
* O motor HPC consolidado usa a Vector Statistics Library para as estatisticas
  iniciais, `cblas_dgemm` para a matriz de distancias
  (`D^2 = ||A||^2 + ||B||^2 - 2 A B^T`), `vdrnguniform` para a interpolacao do
  ADASYN e laco SIMD com FMA (`vfmadd213pd`) para a reversao do z-score.
* As tres funcoes controlam temporariamente apenas `MKL_NUM_THREADS` e
  `OMP_NUM_THREADS`, restaurando os valores originais por um bloco `on.exit()`
  inflexivel. As demais variaveis de ambiente ficam sob controle do servidor.
* O atalho HPC substitui internamente a rota `sby_knn_engine = "native"` como
  caminho rapido quando o motor consolidado esta compilado e carregado. As
  funcoes originais `sby_adanear()`, `sby_adasyn()` e `sby_nearmiss()` continuam
  acessiveis e o fluxo classico permanece inalterado quando o motor HPC nao esta
  disponivel ou quando se pede auditoria, restauro de tipos ou escala intermediaria.

## Correções de contrato KNN

* A rota de compatibilidade `sby_knn_engine = "FNN"` com
  `sby_knn_algorithm = "brute"` agora reutiliza a engine `native` parametrizada,
  preservando `sby_query_is_data`, `sby_exclude_self`, `sby_knn_return` e offsets
  de chunks da mesma forma que `sby_knn_engine = "native"`.
* Chamadas FNN com `sby_knn_algorithm = "auto"` agora resolvem o marcador para
  um algoritmo aceito por `FNN::get.knnx()` antes da chamada externa.
* Validadores de workers, tamanho de chunk e parâmetros HNSW passaram a rejeitar
  valores fracionários em vez de truncá-los silenciosamente.
* As rotas RcppParallel nativas agora preservam o formato de retorno parcial,
  retornando apenas `nn.index` ou apenas `nn.dist` quando solicitado.


## Mudancas de comportamento (breaking)

* `sby_knn_engine = "native"` foi adicionado como engine explicito para KNN
  euclidiano exato em matriz double densa, retornando `nn.index` e/ou
  `nn.dist` no mesmo contrato usado por `FNN::get.knnx()`.
* `sby_knn_engine = "auto"` agora e conservador para ADASYN e NearMiss:
  - para `sby_knn_distance_metric = "euclidean"`, prefere `native` quando a
    biblioteca nativa esta carregada e usa `FNN` como fallback exato;
  - metricas nao euclidianas (`"cosine"`, `"ip"`) nao selecionam busca
    aproximada automaticamente, exceto quando a opcao
    `options(sbyadanear.sby_knn_allow_approx = TRUE)` e ativada ou quando o
    usuario escolhe `sby_knn_engine = "RcppHNSW"` explicitamente.
  Para preservar o comportamento anterior, passe `sby_knn_engine = "FNN"` ou
  `sby_knn_engine = "RcppHNSW"` explicitamente.

## Correcoes

* A engine `native` agora preserva `sby_knn_query_chunk_size` mesmo quando
  `sby_exclude_self = TRUE`, passando o offset global da query para o kernel C++
  para remover self-neighbors corretamente sem forcar uma consulta unica gigante.
* Em ambientes Intel oneAPI/MKL, a configuracao temporaria de threads agora
  atua somente sobre `MKL_NUM_THREADS` e `OMP_NUM_THREADS`.

* `sby_adasyn_matrix()`, `sby_nearmiss_matrix()` e `sby_nearmiss_index()`
  agora rejeitam classes minoritarias com menos de duas observacoes,
  alinhando o contrato da API de matriz com o da API tabular. Antes, a API
  de matriz aceitava `n_minority = 1` e gerava amostras sinteticas
  invalidas (interpolacao do ponto consigo mesmo).
* Teste `native NearMiss selector matches R fallback without ties` deixou
  de usar `.Call("OU_SelectNearMissMajorityC", ..., PACKAGE = "sbyadanear")`
  (forma proibida quando `R_useDynamicSymbols(FALSE)`); agora chama o
  simbolo nativo registrado pelo namespace.
* Mensagem mais clara quando uma coluna preditora se chama literalmente
  `TARGET` (a coluna `TARGET` e reservada para a saida).

## Desempenho

* `sby_knn_parallel_backend = "RcppParallel"` agora registra nos diagnósticos
  `sby_knn_parallel_runtime`, indicando se o runtime efetivo é TBB/oneTBB ou
  TinyThread. O pacote não expõe um parâmetro separado para oneTBB porque essa
  decisão é feita pela instalação do `RcppParallel`.

* API tabular `sby_adasyn()` ficou ~24x mais rapida em `n = 10.000, p = 20`
  ao substituir o loop linha-a-linha (`[[.data.frame` + `do.call(rbind, ...)`)
  por indexacao vetorizada em `sby_build_preserved_predictors()`.
* `sby_drop_self_neighbor_index()` agora usa o kernel C
  `OU_DropSelfNeighborC` por padrao quando a lib nativa esta carregada, e
  um caminho R vetorizado em fallback (era um loop linha-a-linha).
* `sby_generate_adasyn_samples()` substituiu o ciclo `factor -> character
  -> factor` por manipulacao direta dos codigos inteiros do alvo.
* `sby_get_binary_class_roles()` substituiu `table()` por `tabulate()`,
  preservando o shape de retorno.
* Nova rota brute force exata via BLAS: `OU_BruteForceKnnC` (`dgemm` +
  max-heap top-k) e usada por padrao quando `sby_knn_algorithm = "brute"`
  e a lib nativa esta carregada. Pode ser desativada com
  `options(sbyadanear.sby_use_native_brute = FALSE)`.
* A rota brute force nativa agora possui variantes index-only e dist-only,
  evitando alocar componentes KNN descartados pelas etapas ADASYN e NearMiss.
* NearMiss-1 recebeu uma rota exata fundida para `FNN` + `brute`, calculando
  as medias dos k vizinhos minoritarios via BLAS e selecionando a maioria sem
  materializar `nn.dist` em R.
* `sby_knn_query_chunk_size` passou a ser argumento publico nas APIs tabulares,
  matriciais e de `recipes`, permitindo ajustar o tamanho dos blocos KNN sem
  depender de opcoes globais.
* Variante experimental do kernel ADASYN com escrita column-friendly:
  `OU_GenerateSyntheticAdasynColC`, disponivel via
  `options(sbyadanear.sby_adasyn_kernel = "col")`. Em testes empiricos
  vence em `p` muito alto (> 200) e `n_synthetic` muito alto (> 10^5);
  o default (`"row"`) continua sendo a variante anterior.

## Empacotamento

* `DESCRIPTION`: `Language: pt_BR` corrigido para `Language: pt-BR`
  (BCP47 valido); `glue` removido de `Imports` (nao era usado); `Rfast`
  movido para `Suggests` (so usado em fallback puro-R).
* `R CMD check`: removidos avisos de non-ASCII em
  `R/sby_resolve_knn_engine.R` e `R/sby_resolve_knn_algorithm.R`, e as
  duas subsecoes vazias em `man/*.Rd`.
* README: versao documentada atualizada para 0.3.0.

# sbyadanear 0.3.0

## Novidades

* API matricial publica para ADASYN, NearMiss-1 e ADANEAR:
  `sby_adasyn_matrix()`, `sby_nearmiss_matrix()` e `sby_adanear_matrix()`.
* Novo roteador industrial `sby_balance_matrix()` para estrategias
  `none`, `weight`, `adasyn`, `nearmiss`, `adanear` e `adanearWeight`.
* `sby_nearmiss_index()` expoe selecao por indices para fluxos que
  precisam preservar matrizes externas, inclusive esparsas fora do pacote.
* Selecao NearMiss parcial em C com criterio deterministico por menor
  distancia media e menor indice em empates.
* Auditoria leve para diagnosticos resumidos sem matrizes intermediarias
  pesadas.
* Melhor controle de memoria: `none` e `weight` nao densificam entradas,
  KNN em blocos aloca apenas componentes solicitados e NearMiss evita
  recalculo de z-score.
* Wrappers tabulares preservam linhas originais e desnormalizam/restauram
  tipos apenas nas linhas sinteticas do ADASYN.
