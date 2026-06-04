# Auditoria tecnica senior de desempenho e corretude

Data: 2026-06-04.

## Sintese executiva

O pacote ja concentra os caminhos quentes em matriz densa, KNN nativo opcional, RcppParallel, FNN e RcppHNSW. A auditoria encontrou uma oportunidade segura aplicada imediatamente no KNN bruto nativo: a rotina de resgate numerico para distancias quase zero alocava tres buffers temporarios via `R_alloc()` para cada par resgatado, e os vetores de vizinhos eram realocados a cada linha de consulta. O patch remove essas alocacoes por par e reutiliza o buffer de vizinhos por bloco ou tarefa paralela.

A API publica foi preservada. Nenhuma funcao exportada foi renomeada. As demais oportunidades envolvendo troca de backend KNN, matriz completa de similaridade ou novas dependencias foram classificadas como investigar ou nao aplicar por risco de memoria, semantica ou CRAN.

## Bugs e riscos encontrados

| Item | Severidade | Arquivo | Funcao | Achado | Patch |
| --- | --- | --- | --- | --- | --- |
| Alocacao no resgate numerico de distancia quase zero | Media | `src/sbyadanear.cpp` | `euclidean_d2_with_rescue()` | Pares identicos ou quase identicos acionavam `R_alloc()` repetidamente dentro do loop KNN. Em self-KNN, isso ocorre em muitos pontos e aumenta pressao de memoria e tempo. | Aplicado |
| Realocacao de vetor de vizinhos por query | Media | `src/sbyadanear.cpp` | `brute_force_knn_impl()`, `nearmiss_brute_select_c()` | `std::vector<knn_neighbor>` era criado para cada linha de consulta. | Aplicado |
| Realocacao por query no caminho RcppParallel | Media | `src/sbyadanear_rcppparallel.cpp` | `sby_brute_force_knn_worker::operator()` | Cada query alocava um vetor de vizinhos dentro da tarefa paralela. | Aplicado |
| Oversubscription potencial | Media | `R/sby_get_knnx.R`, `R/sby_oneapi_mkl.R` | `sby_get_knnx()`, `sby_configure_blas_threads()` | Ha controle de BLAS e workers, mas FNN por blocos, RcppHNSW e RcppParallel podem competir com BLAS/MKL se o usuario configurar workers altos. | Investigar com benchmark por maquina |
| Materializacao densa inevitavel | Alta | `R/sby_adanear_as_numeric_matrix.R`, `R/sby_validate_dense_double_matrix.R` | Conversao e validacao de matriz | O pacote bloqueia Matrix esparsa, o que evita densificacao acidental. O risco remanescente e o usuario converter antes da chamada. | Manter bloqueio |
| Troca semantica em KNN aproximado | Media | `R/sby_resolve_knn_engine.R`, `R/sby_get_knnx.R` | Resolucao de engine | RcppHNSW e metricas angulares mudam vizinhos e podem alterar ADASYN/NearMiss. | Apenas benchmark e testes de equivalencia quando aplicavel |

## Duplicacoes e gargalos

| Tipo | Arquivo | Funcao | Observacao | Acao |
| --- | --- | --- | --- | --- |
| Validacao repetida | `R/sby_adasyn_matrix.R`, `R/sby_nearmiss_index.R`, `R/sby_adanear_matrix.R` | APIs matriciais | Validacoes de fator, workers, engine e escala sao repetidas, mas legiveis e baratas frente ao KNN. | Manter por clareza |
| KNN bruto completo | `src/sbyadanear.cpp` | `brute_force_knn_impl()` | Calcula bloco de produto cruzado e top-k por query. O caminho e quente e sensivel a alocacoes. | Patch aplicado |
| NearMiss com distancias | `R/sby_nearmiss_index.R` | `sby_nearmiss_index()` | Caminho nao nativo materializa `nn.dist` para depois tirar media. | Manter nativo quando brute euclidiano; investigar streaming para outros engines |
| Conversao data frame para matriz | `R/sby_adanear_as_numeric_matrix.R` | `sby_adanear_as_numeric_matrix()` | `data.matrix()` e seguro, mas pode copiar grande volume. | Benchmark opcional com `collapse::qM()` |

## Tabela pacote a pacote

| Pacote | Funcoes candidatas | Arquivo do sbyadanear afetado | Funcao afetada | Classificacao | Ganho esperado | Risco | Acao recomendada | Precisa benchmark | Precisa teste |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RcppParallel | `parallelFor`, TBB | `src/sbyadanear_rcppparallel.cpp` | `brute_force_knn_rcpp_parallel_c()` | Alta | Melhor paralelismo, menor overhead apos reuse de vetor | Oversubscription com MKL | Manter Imports/LinkingTo e benchmarkar workers | Sim | Sim |
| kit | `iif`, `nif`, `topn`, `funique` | `R/sby_restore_numeric_column_types.R`, `R/sby_resolve_knn_engine.R` | Restauracao e roteamento | Baixa | Menor overhead em operacoes vetorizadas pequenas | Dependencia ja existe | Manter uso atual; nao migrar KNN | Nao | Sim |
| Rfast | `rowmeans`, `rowsums`, `eachrow` | `R/sby_generate_adasyn_samples.R`, `R/sby_nearmiss_index.R`, `R/sby_apply_z_score_scaling_matrix.R` | Razoes, medias, z-score | Media | Velocidade em matriz densa | Threads internos podem competir | Manter, benchmarkar com MKL threads | Sim | Sim |
| coop | `covar`, `cosine` | `R/sby_fast_stat_helpers.R`, KNN potencial | Diagnostico, nao KNN | Risco | BLAS rapido para matriz completa | Estouro de memoria `n x n` | Nao usar para top-k em producao | Sim se diagnostico | Sim |
| collapse | `fnrow`, `fncol`, `fdim`, `qM` | APIs R | Validacao e conversao | Media | Menor overhead de dimensoes e conversao | `qM` pode mudar casting | Manter atual; testar `qM()` em benchmark | Sim | Sim |
| data.table | `rbindlist` | `R/sby_build_preserved_predictors.R` | Montagem tabular | Media | Menor copia na montagem | Ordem e atributos | Manter Suggests/Imports atual | Sim | Sim |
| fastmap | Cache | Nenhum atual | Nenhuma | Nao aplicavel | Nenhum ganho claro | Cache invalido por seed/parametros | Nao adicionar | Nao | Nao |
| cheapr | Conversoes rapidas | Conversao tabular | `sby_adanear_as_numeric_matrix()` | Baixa | Possivel menor copia | Dependencia pouco essencial | Nao adicionar agora | Sim | Sim |
| r2c | Compilacao R para C | Loops R fallback | Fallback ADASYN | Risco | Velocidade teorica | Manutencao e CRAN | Nao adicionar | Sim | Sim |
| nCompiler | Compilacao | Loops R fallback | Fallback ADASYN | Risco | Teorico | Dependencia/portabilidade | Nao adicionar | Sim | Sim |
| tidycpp | C++ tidy | Nenhum | Nenhuma | Nao aplicavel | Nenhum | API externa | Nao adicionar | Nao | Nao |
| dtplyr | Lazy data.table | Nenhum pipeline dplyr | Nenhuma | Nao aplicavel | Nenhum | Dependencia extra | Nao adicionar | Nao | Nao |
| fastplyr | Pipelines | Nenhum | Nenhuma | Nao aplicavel | Nenhum | Dependencia extra | Nao adicionar | Nao | Nao |
| tidyfast | Pipelines | Nenhum | Nenhuma | Nao aplicavel | Nenhum | Dependencia extra | Nao adicionar | Nao | Nao |
| tidytable | Pipelines | Nenhum | Nenhuma | Nao aplicavel | Nenhum | Dependencia extra | Nao adicionar | Nao | Nao |
| vctrs | Casting/restauracao | `R/sby_restore_numeric_column_types.R` | Restauracao | Baixa | Estabilidade de casting | Dependencia extra, overhead | Investigar apenas se API tibble quebrar | Sim | Sim |
| rrapply | Auditoria aninhada | Diagnosticos | Listas de auditoria | Baixa | Codigo menor | Sem ganho hot path | Nao adicionar | Nao | Nao |
| vcrpart | Particionamento | Nenhum | Nenhuma | Nao aplicavel | Nenhum | Sem encaixe em KNN | Nao adicionar | Nao | Nao |
| bbknnR | KNN batch | `R/sby_get_knnx.R` | KNN | Risco | Possivel velocidade | Semantica de batch correction nao e ADASYN/NearMiss | Nao substituir | Sim | Sim |
| tsfknn | KNN series temporais | Nenhum | Nenhuma | Nao aplicavel | Nenhum | Dominio errado | Nao adicionar | Nao | Nao |
| knn.covertree | Cover tree | `R/sby_get_knnx.R` | KNN | Baixa | Alternativa exata | API e manutencao | Investigar somente benchmark | Sim | Sim |
| KernelKnn | `knn.index.dist` | `R/sby_get_knnx.R` | KNN | Media | Alternativa euclidiana OpenMP | Mudanca de ordenacao/threads | Implementado como engine opcional em Suggests | Sim | Sim |
| bigKNN | `knn_bigmatrix` | `R/sby_get_knnx.R` | KNN | Media | Busca exata em bigmemory por blocos | Dependencia opcional/API externa | Implementado como engine opcional em Suggests | Sim | Sim |
| fastverse | Colecao de pacotes | Nenhum | Nenhuma | Nao aplicavel | Nenhum | Dependencia ampla | Nao adicionar | Nao | Nao |
| ClassifyR | Pipelines classificacao | Nenhum | Nenhuma | Nao aplicavel | Validacao externa | Escopo diferente | Nao adicionar | Nao | Nao |
| CRAN Task View HPC | OpenMP, BLAS, bigmemory | `src`, `R/sby_oneapi_mkl.R` | Nativo e threads | Media | Guia de escolhas | Nao e dependencia | Usar como referencia | Sim | Nao |
| mlr3pipelines | PipeOps | Recipes comparacao | Steps | Baixa | Ideias de API | Dependencia pesada | Nao substituir recipes | Sim | Sim |
| mlr3 PipeOp NearMiss | NearMiss referencia | `R/sby_nearmiss_index.R` | NearMiss | Baixa | Validacao semantica | API externa | Usar somente em testes comparativos opcionais | Sim | Sim |
| smotefamily | ADASYN referencia | `R/sby_generate_adasyn_samples.R` | ADASYN | Media | Teste de comportamento | Resultados aleatorios diferentes | Suggests opcional em testes, nao producao | Sim | Sim |
| bugsparallel | Paralelismo | Nenhum | Nenhuma | Nao aplicavel | Nenhum | Escopo incerto | Nao adicionar | Nao | Nao |
| quickr | Aceleracao R | Nenhum | Nenhuma | Nao aplicavel | Incerto | Portabilidade | Nao adicionar | Nao | Nao |
| tinythread++ | Threads C++ | `src` | KNN nativo | Risco | Alternativa leve | Ja ha RcppParallel/TBB | Nao adicionar | Sim | Sim |

## Tabela funcao por funcao

| Funcao | Arquivo | Problema encontrado | Tipo | Severidade | Otimizacao proposta | Impacto esperado | Risco semantico | Teste necessario | Benchmark necessario | Aplicar agora ou investigar |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `sby_adasyn` | `R/sby_adasyn.R` | Caminho tabular depende de conversao para matriz | Performance | Baixa | Benchmark `collapse::qM()` | Menor copia | Casting | Formula e tibble | Sim | Investigar |
| `sby_nearmiss` | `R/sby_nearmiss.R` | Mesmo custo tabular | Performance | Baixa | Mesmo acima | Menor copia | Casting | Formula e tibble | Sim | Investigar |
| `sby_adanear` | `R/sby_adanear_sampling.R` | Combina duas passagens KNN | Performance | Media | Manter reuso de escala ja existente | Menor recomputacao | Baixo | Seed e auditoria | Sim | Manter |
| `sby_adasyn_matrix` | `R/sby_adasyn_matrix.R` | KNN e geracao sintetica sao hot path | Performance | Media | Patch no KNN nativo e benchmark | Velocidade e memoria | Baixo | KNN brute/RcppParallel | Sim | Aplicado parcialmente |
| `sby_nearmiss_matrix` | `R/sby_nearmiss_matrix.R` | Depende de `sby_nearmiss_index()` | Performance | Media | Patch no selector nativo | Menor memoria | Baixo | Distribuicoes | Sim | Aplicado parcialmente |
| `sby_adanear_matrix` | `R/sby_adanear_matrix.R` | Potencial oversubscription em duas etapas | Execucao | Media | Benchmark workers/MKL | Estabilidade | Medio | workers 1/2 | Sim | Investigar |
| `sby_balance_matrix` | `R/sby_balance_matrix.R` | Wrapper preserva estrategias | API | Baixa | Nao alterar | Nenhum | Alto se alterar | Estrategias | Sim | Manter |
| `sby_nearmiss_index` | `R/sby_nearmiss_index.R` | Caminho nao nativo materializa distancias | Memoria | Media | Streaming para engines futuros | Menor memoria | Medio | Ordem e ties | Sim | Investigar |
| `sby_step_adasyn` | `R/sby_step_adasyn.R` | Overhead recipes fora hot path | API | Baixa | Nao alterar | Estabilidade | Alto se alterar | Recipes | Nao | Manter |
| `sby_step_nearmiss` | `R/sby_step_nearmiss.R` | Mesmo acima | API | Baixa | Nao alterar | Estabilidade | Alto | Recipes | Nao | Manter |
| `sby_step_adanear` | `R/sby_step_adanear.R` | Mesmo acima | API | Baixa | Nao alterar | Estabilidade | Alto | Recipes | Nao | Manter |
| `euclidean_d2_with_rescue` | `src/sbyadanear.cpp` | `R_alloc()` por par resgatado | Memoria/performance | Media | Somar diretamente em `long double` | Menor alocacao | Baixo | KNN nativo | Sim | Aplicado |
| `brute_force_knn_impl` | `src/sbyadanear.cpp` | Vetor por query | Memoria/performance | Media | Reusar vetor por bloco | Menor alocacao | Baixo | Simbolos C | Sim | Aplicado |
| `nearmiss_brute_select_c` | `src/sbyadanear.cpp` | Vetor por query | Memoria/performance | Media | Reusar vetor por bloco | Menor alocacao | Baixo | NearMiss | Sim | Aplicado |
| `sby_brute_force_knn_worker::operator()` | `src/sbyadanear_rcppparallel.cpp` | Vetor por query | Memoria/performance | Media | Reusar vetor por tarefa | Menor alocacao | Baixo | RcppParallel | Sim | Aplicado |

## Patches propostos e aplicados

1. Remover buffers temporarios do resgate numerico de distancia quase zero em `src/sbyadanear.cpp`.
2. Reusar `std::vector<knn_neighbor>` por bloco no KNN bruto nativo em `src/sbyadanear.cpp`.
3. Reusar `std::vector<sby_neighbor>` por tarefa no worker RcppParallel em `src/sbyadanear_rcppparallel.cpp`.
4. Adicionar teste de regressao para caminho nativo com distancias zero e caminho RcppParallel.
5. Adicionar benchmark focado em pressao de alocacao do KNN nativo.

## Testes e benchmarks novos

| Arquivo | Objetivo |
| --- | --- |
| `tests/testthat/test-native-allocation-regression.R` | Verifica que NearMiss nativo lida com distancias zero e que ADASYN no caminho RcppParallel retorna formatos coerentes. |
| `tools/benchmarks/benchmark-native-knn-allocation.R` | Mede tempo e contagens finais em cenarios pequeno, medio e alta dimensionalidade para KNN bruto NearMiss. |

## Dependencias

| Acao | Dependencia | Justificativa |
| --- | --- | --- |
| Manter Imports/LinkingTo | RcppParallel | Necessario para kernel nativo paralelo e deteccao TBB. |
| Manter Imports | FNN, RcppHNSW | Engines KNN publicos. |
| Manter Imports | collapse, kit, Rfast, data.table, coop | Ja usados em caminhos reais ou helpers. |
| Nao adicionar | fastmap, dtplyr, fastplyr, tidyfast, tidytable, rrapply, fastverse | Sem ganho no hot path atual. |
| Adicionado como engine opcional | KernelKnn, bigKNN | Solicitados para benchmark sem alterar default; exigem teste semantico antes de uso operacional. |
| Nao adicionar agora | knn.covertree, smotefamily, mlr3pipelines | Exigem benchmark e teste semantico antes. |

## Mudancas que nao devem ser feitas agora

- Nao substituir top-k por `coop::cosine()` ou correlacao completa, pois isso materializa matriz `n_query x n_ref`.
- Nao adicionar dependencia pesada para conversao tabular antes de confirmar ganho em benchmark.
- Nao ativar paralelismo aninhado entre `parallel`, RcppParallel, HNSW e MKL sem controle explicito de threads.
- Nao trocar FNN brute exato por HNSW aproximado como padrao, pois muda vizinhos e resultados sinteticos.

## Plano de implementacao

1. Concluir patch de alocacao nativa ja aplicado.
2. Rodar `R CMD check`, `testthat` e benchmark nativo em ambiente com R.
3. Medir workers 1 e 2 com BLAS/MKL em 1 thread e multiplas threads.
4. Investigar `collapse::qM()` contra `data.matrix()` em data frames grandes.
5. Criar benchmark comparativo opcional para KernelKnn/bigKNN/knn.covertree sem alterar defaults.
6. Revisar se `coop` deve permanecer em Imports ou migrar para Suggests se o helper nao for API necessaria.
