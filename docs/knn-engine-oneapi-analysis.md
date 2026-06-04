# Analise de engines KNN, RcppParallel, bigmemory e oneAPI

Data: 2026-06-04.

## Fontes consultadas

- Rcpp Gallery, `parallelFor` com `RcppParallel::Worker` e `RMatrix` thread safe para ler e escrever matrizes em threads de fundo.
- Rcpp Gallery, `bigmemory` com `XPtr<BigMatrix>` e `MatrixAccessor`, destacando acesso column-major por `mat[col][row]`.
- Manual `KernelKnn`, funcao `knn.index.dist()` com `data`, `TEST_data`, `k`, `method` e `threads` via OpenMP.
- Manual `bigKNN`, funcao `knn_bigmatrix()` para busca exata em `bigmemory::big.matrix`, consulta densa opcional, `block_size` e saida `index`/`distance`.
- Repositorio e pagina `fastmatrix`, com foco em rotinas matriciais estatisticas, decomposicoes e interfaces C, mas sem backend top-k KNN.

## Aplicacao ao sbyadanear

O padrao do artigo de RcppParallel confirma que o caminho nativo deve encapsular trabalho independente por faixa de linhas em um `Worker`, ler dados por `RMatrix` e evitar chamadas R/Rcpp dentro das threads. O pacote ja segue essa direcao em `sby_brute_force_knn_worker`; o ajuste aplicado reutiliza o buffer de vizinhos por tarefa, o que combina melhor com TBB, oneTBB e tbbmalloc.

O artigo de bigmemory reforca que acesso eficiente a `big.matrix` em C/C++ e column-major e que a ponte correta e um ponteiro externo. A integracao feita aqui usa `bigKNN` de forma opcional, convertendo a referencia densa para `bigmemory::big.matrix` e delegando a busca exata a `bigKNN::knn_bigmatrix()` quando o usuario seleciona `sby_knn_engine = "bigKNN"`.

`KernelKnn` foi integrado como engine opcional via `KernelKnn::knn.index.dist()` para distancia euclidiana. Ele usa OpenMP internamente, portanto deve ser benchmarkado com `sby_knn_workers` controlado e com MKL/oneTBB/TCM sem paralelismo externo simultaneo.

`fastmatrix` nao foi adicionado como dependencia. As rotinas publicas revisadas sao voltadas a matrizes estruturadas, decomposicoes, regressao, covariancia, Mahalanobis e interfaces C. Elas nao substituem diretamente o gargalo atual, que e top-k KNN streaming e geracao ADASYN.

## Ranking estimado de gargalos

| Ranking | Gargalo | Motivo | Acao atual |
| --- | --- | --- | --- |
| 1 | KNN global do ADASYN | Consulta cada linha minoritaria contra todo o conjunto para estimar dificuldade local. Custo cresce com `n_minority * n * p`. | Manter kernel nativo brute, FNN, RcppHNSW, KernelKnn e bigKNN como rotas benchmarkaveis. |
| 2 | KNN NearMiss maioria contra minoria | Consulta cada linha majoritaria contra a minoria e ordena medias para selecionar retencao. | Usar selector nativo brute quando euclidiano e investigar bigKNN em bases grandes. |
| 3 | KNN minoritario do ADASYN | Self-KNN da minoria para interpolacao sintetica. Pode acionar muitos empates e distancias zero. | Patch removeu alocacoes no resgate numerico e reusa buffers. |
| 4 | Padronizacao e reversao de escala | Passa por toda matriz densa e cria matriz escalada. | Manter Rfast e avaliar BLAS/MKL threading. |
| 5 | Montagem do retorno tabular | `as.data.frame`, restauração de tipos, `rbindlist` e tibble podem copiar grande volume apos KNN. | Usar `collapse::qM()` na conversao tabular numerica e manter montagem fora do hot path quando API matricial e usada. |
| 6 | Geração sintetica ADASYN | Interpolacao depende do numero de sinteticos e `p`; menos pesada que KNN em bases grandes. | Manter kernel nativo e reavaliar caminho column-friendly apenas por benchmark. |
| 7 | Auditoria full | Pode reter artefatos grandes, indices e escalas. | Usar `sby_audit_level = "light"` para diagnostico operacional. |

## Regras de uso recomendadas

- Para exatidao euclidiana e alta dimensionalidade moderada, iniciar com `sby_knn_engine = "FNN"` e `sby_knn_algorithm = "brute"`, usando o kernel nativo quando disponivel.
- Para comparar OpenMP externo, usar `sby_knn_engine = "KernelKnn"` somente com `sby_knn_distance_metric = "euclidean"`.
- Para bases que ja justificam `bigmemory`, testar `sby_knn_engine = "bigKNN"` com `sby_knn_query_chunk_size` alinhado a cache e memoria disponivel.
- Para cosine ou produto interno, manter `RcppHNSW`, pois as novas rotas foram limitadas a euclidean para evitar divergencia semantica.
- Em oneAPI, evitar combinar simultaneamente MKL multithread, OpenMP de KernelKnn, TBB de RcppParallel e workers de `parallel`. Escolher um nivel principal de paralelismo por benchmark.
