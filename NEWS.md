# sbyadanear 0.4.0 (em desenvolvimento)

## Correcoes do motor HPC

* Compilacao: `src/Makevars` e `src/Makevars.win` passam a propagar
  `$(SHLIB_OPENMP_CXXFLAGS)` e `$(SHLIB_OPENMP_FCFLAGS)`/`$(SHLIB_OPENMP_FFLAGS)`
  para C++, Fortran e link. Sem isso, todo `#pragma omp` e todo `!$omp` eram
  descartados em silencio e `sby_config_max_threads` nao tinha efeito algum.
* Compilacao: o oneMKL agora e de fato ligado (`-lmkl_rt`) quando `MKLROOT` ou
  `ONEAPI_ROOT` apontam para uma instalacao valida. A macro
  `SBYADANEAR_ONEAPI_MKL` so e definida nesse caso, e `sby_hpc_cpu_report()`
  expoe `compile_report$mkl_linked` para auditoria. Sem MKL, o pacote liga
  `$(BLAS_LIBS)`.
* Compilacao: as flags de arquitetura Cascade Lake sao detectadas por sondagem
  do compilador (`-march=cascadelake` no GCC/Clang, `-xCORE-AVX512` no Intel),
  com escape por `SBYADANEAR_ARCH_FLAGS` e `SBYADANEAR_NO_ARCH_FLAGS`.
* ADASYN: a alocacao das linhas sinteticas no motor HPC era round-robin sobre a
  minoria. Agora segue o ADASYN classico, ponderada pela fracao de vizinhos
  majoritarios de cada ponto raro (`r_i`), com o mesmo desempate por maior
  residuo da rota em R. Pontos raros sem vizinhanca majoritaria deixam de
  receber sinteticas.
* ADASYN e NearMiss: `Rcpp::runif()` passa a ser chamado sob `Rcpp::RNGScope`,
  o que restabelece a reprodutibilidade por `sby_seed` nas rotas HPC.
* Threads: `omp_set_num_threads()` era global e permanente. O motor agora usa um
  guarda RAII que restaura `omp_get_max_threads()` ao sair da chamada.
* Contagens: o motor HPC nao arredonda mais para cima uma linha sintetica nem um
  registro majoritario quando a razao pedida resulta em zero. O comportamento
  passa a ser identico ao da rota classica em R; `sby_nearmiss_ratio` que zera a
  maioria agora aborta com mensagem explicita.
* kNN: a busca exata deixa de preencher com o indice 1 quando encontra menos de
  `k` candidatos validos. As posicoes nao preenchidas usam sentinela `-1` e a
  chamada aborta em vez de enviesar as sinteticas para a primeira linha.
* Distancias: as normas `||A||^2` e `||B||^2` sao acumuladas em precisao dupla.
  Em float32, o cancelamento da identidade euclidiana produzia distancias
  negativas que o clamp em zero mascarava, invertendo a ordem dos vizinhos.
* Formulas: `sby_adanear_hpc()`, `sby_adasyn_hpc()` e `sby_nearmiss_hpc()`
  reordenam apenas as colunas presentes na saida. Formulas diferentes de
  `y ~ .` abortavam em `collapse::fselect()`.
* Validacao: `sby_nearmiss_ratio` em `sby_adanear_hpc()` agora rejeita `Inf`,
  como ja fazia `sby_adasyn_ratio`.
* Validacao: `sby_binary_class_counts_fast()` aborta com classes perfeitamente
  balanceadas, em vez de deixar `which.min()` e `which.max()` colapsarem no
  mesmo nivel. A mensagem e a mesma de `sby_get_binary_class_roles()`.
* Documentacao: removidas as mencoes a MKL VSL, `cblas_sgemm` e ao restauro de
  variaveis de ambiente por `on.exit()`, que nao correspondiam ao codigo.
  `cascade_lake_native` no relatorio de compilacao passa a refletir as macros
  AVX-512 reais do compilador.

## Desempenho do motor HPC

* O laco de manutencao do top-k do ADASYN era serial; agora e paralelo em
  OpenMP, como o equivalente do NearMiss.
* A blocagem do `sgemm` passa a ter piso de tile e a particionar as duas
  dimensoes, trocando milhares de chamadas finas por poucas chamadas densas.
* Os buffers deixam de ser zerados em paralelo logo antes de serem sobrescritos
  por `beta = 0`. O *first touch* NUMA e feito por pagina, sem inicializar todos
  os elementos.
* Os blocos de distancia sao orientados como (referencia x consulta), o que
  torna contigua a leitura do laco interno tanto no kNN quanto no NearMiss.
* A despadronizacao das sinteticas escreve direto na matriz de saida, sem buffer
  intermediario nem copia.
* `sby_hpc_resolve_threads()` respeita cotas de cgroup v1 e v2, evitando
  oversubscription em containers. Removidos o resolvedor morto de
  `MKL_NUM_STRIPES` e o injetor de ambiente que nenhuma rota chamava.
* README com orientacao NUMA para servidores de dois sockets
  (`OMP_PROC_BIND`, `OMP_PLACES`, `numactl --interleave=all`).

## Atalho HPC (oneAPI, AVX-512)

* Novas funcoes exportadas `sby_adanear_hpc()`, `sby_adasyn_hpc()` e
  `sby_nearmiss_hpc()`. Cada uma e um atalho de alto desempenho que executa o
  fluxo estritamente no espaco padronizado, eliminando a dupla normalizacao, e
  monta o tibble final por zero-copy diretamente em C++ via `Rcpp::List`.
* O motor HPC consolidado usa lacos SIMD paralelos para as estatisticas
  iniciais, `sgemm` para a matriz de distancias
  (`D^2 = ||A||^2 + ||B||^2 - 2 A B^T`), `Rcpp::runif()` sob `RNGScope` para a
  interpolacao do ADASYN e laco SIMD com FMA para a reversao do z-score. A
  interface Fortran padrao do BLAS e usada de proposito, para que o oneMKL seja
  aproveitado quando disponivel sem tornar a BLAS do R um caminho invalido.
* As consultas KNN internas e o NearMiss usam blocagem/streaming de SGEMM
  com top-k incremental para evitar materializar matrizes de distancia completas
  quando os blocos excedem o orcamento interno de memoria.
* As tres funcoes nao alteram variaveis de ambiente do runtime MKL/OpenMP. O
  numero de threads pedido em `sby_config_max_threads` vale apenas para a
  chamada corrente: o motor nativo salva e restaura `omp_get_max_threads()` em
  torno do kernel. A politica de afinidade e de memoria fica sob controle do
  servidor.
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

* `sby_adanear_hpc()` e `sby_nearmiss_hpc()` voltaram a preservar integralmente
  a classe rara original. As duas rotas selecionavam as linhas raras por
  `sby_class_counts$sby_minority_level`, um campo que `sby_binary_class_counts_fast()`
  nunca devolveu; `as.integer(NULL)` produzia `integer(0)`, o `which()` devolvia
  zero indices e a saida ficava apenas com a maioria retida mais as sinteticas.
  Em `sby_adanear_hpc()` isso disparava a pos-condicao
  `sby_assert_minority_not_reduced()` (por exemplo, entrada = 2874 e
  saida = 1149 com `sby_adasyn_ratio = 0.4`); em `sby_nearmiss_hpc()`, que nao
  tinha essa verificacao, a perda era silenciosa. A funcao de contagem agora
  expoe `sby_minority_level` e `sby_majority_level` e `sby_nearmiss_hpc()`
  tambem valida a pos-condicao.
* A engine `native` agora preserva `sby_knn_query_chunk_size` mesmo quando
  `sby_exclude_self = TRUE`, passando o offset global da query para o kernel C++
  para remover self-neighbors corretamente sem forcar uma consulta unica gigante.
* Em ambientes Intel oneAPI/MKL, a configuracao temporaria de threads agora
  atua somente sobre `MKL_NUM_THREADS`, `OMP_NUM_THREADS` e `MKL_NUM_STRIPES`.

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
