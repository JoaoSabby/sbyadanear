# sbyadanear — balanceamento binário ADASYN + NearMiss-1

`sbyadanear` é um pacote R para **engenharia de instâncias em problemas de
classificação binária desbalanceados**. O pacote oferece rotinas de
sobreamostragem **ADASYN**, subamostragem **NearMiss-1** e um pipeline híbrido
**ADASYN + NearMiss-1** chamado `sby_adanear()`.

Internamente, as rotinas trabalham com matrizes numéricas, padronização Z-score,
consultas KNN configuráveis, fallback com `Rfast` e kernels nativos em C para
partes críticas quando disponíveis.

## Estado da API

A versão atual do pacote é **0.3.0**. O pacote está em desenvolvimento e a API
pública foi organizada para usar:

- Funções e parâmetros com prefixo `sby_` e padrão `snake_case`.
- Interface principal por **fórmula + dados**: `sby_formula` e `sby_data`.
- Preditores numéricos e alvo binário com exatamente duas classes.
- Retorno padrão como `tibble` balanceado.
- Coluna de desfecho padronizada como `TARGET` no objeto retornado.
- `sby_audit = FALSE` para retornar apenas os dados balanceados.
- `sby_audit = TRUE` para retornar lista com dados balanceados, diagnósticos,
  informações de escala e metadados de tipos.

> **Importante:** em chamadas como `sby_adanear(sby_y ~ ., sby_data)`, o lado
> esquerdo da fórmula identifica o desfecho e o lado direito identifica os
> preditores. Use `~ .` para usar todas as demais colunas como preditores.

### Contrato de fórmulas e dados

As funções públicas aceitam fórmulas para **selecionar colunas já existentes** em
`sby_data`. Transformações, interações e offsets da sintaxe de fórmulas do R
(por exemplo, `log(x)`, `x1:x2` ou `offset(z)`) devem ser calculados antes da
chamada e armazenados como colunas explícitas. Essa restrição evita divergências
entre a seleção por fórmula, a restauração de tipos e o pipeline matricial usado
pelas rotinas KNN e nativas.

Preditores devem ser numéricos, finitos e densos. Matrizes esparsas do pacote
`Matrix` são rejeitadas com erro explícito para evitar densificação acidental de
bases grandes. Colunas constantes também são rejeitadas porque a padronização
Z-score exige desvio padrão positivo.

Para bases pequenas, `sby_over_ratio` positivo sempre gera ao menos uma linha
sintética; use uma razão positiva explícita e uma `sby_seed` inteira para
reprodutibilidade.

## KNN, métricas e engines

As rotinas usam KNN para estimar vizinhanças locais. Os principais controles são:

- `sby_knn_engine`: engine de busca (`"auto"`, `"native"`, `"FNN"`, `"RcppHNSW"`, `"KernelKnn"`, `"bigKNN"`).
- `sby_knn_algorithm`: algoritmo exato do FNN (`"auto"`, `"kd_tree"`,
  `"cover_tree"`, `"brute"`); a engine `native` aceita `"auto"` ou `"brute"` e resolve internamente a rota exata.
- `sby_knn_distance_metric`: métrica (`"euclidean"`, `"cosine"`, `"ip"`).
- `sby_knn_workers`: número de workers. Em `native`, os workers podem acionar o kernel RcppParallel exato; em `FNN`, consultas exatas são
  paralelizadas por blocos; em `RcppHNSW`, os workers são repassados aos
  threads nativos do índice aproximado.
- `sby_knn_parallel_backend`: backend do paralelismo exato. Use `"parallel"`
  para manter o particionamento por blocos do R ou `"RcppParallel"` para
  acionar threads nativos no kernel exato bruto (`native` ou compatibilidade `FNN` + `brute`).
- `sby_knn_hnsw_m` e `sby_knn_hnsw_ef`: parâmetros do HNSW quando
  `sby_knn_engine = "RcppHNSW"`.
- `sby_knn_query_chunk_size`: quantidade de linhas de consulta processadas por
  bloco nas rotas KNN. O padrão `1000L` equilibra overhead de chamadas e pico
  de memória; valores maiores podem favorecer BLAS/MKL em matrizes densas,
  enquanto valores menores reduzem pressão de memória.

Resumo de compatibilidade:

| Engine | Tipo de busca | Métricas suportadas | Algoritmos aceitos no pacote |
|---|---|---|---|
| `native` | Exata densa via kernel C/C++ interno | `euclidean` | `auto`, `brute` |
| `FNN` | Exata via `FNN::get.knnx()` | `euclidean` | `auto`, `kd_tree`, `cover_tree`, `brute` |
| `RcppHNSW` | Aproximada por HNSW | `euclidean`, `cosine`, `ip` | `auto` |
| `KernelKnn` | Exata via `KernelKnn::knn.index.dist()` | `euclidean` | `auto`, `brute` |
| `bigKNN` | Exata via `bigKNN::knn_bigmatrix()` | `euclidean` | `auto`, `brute` |

A seleção automática evita busca aproximada por padrão. Com
`sby_knn_engine = "auto"` e `sby_knn_distance_metric = "euclidean"`, o pacote
prefere `native` quando as rotinas nativas estão carregadas; se elas não estiverem
disponíveis, usa `FNN` como fallback exato. Para permitir que `auto` escolha
`RcppHNSW` em métricas não euclidianas quando a aproximação é aceitável, use
`options(sbyadanear.sby_knn_allow_approx = TRUE)`. Essa configuração é uma opção
global, não um argumento público das funções.

O contrato interno comum de KNN é uma lista com `nn.index` e/ou `nn.dist`:
`nn.index` usa índices 1-based compatíveis com R, e `nn.dist` representa a
distância retornada pela engine efetiva. Nas rotas exatas euclidianas nativas, as
distâncias são a raiz quadrada da soma de quadrados, não a distância quadrática.
Para `cosine` e `ip`, o pacote normaliza as linhas por norma L2 antes da busca e
usa a escala de distância retornada por `RcppHNSW`.


Consultas KNN longas são executadas em blocos para permitir interrupção por
`Ctrl + C` entre blocos e para controlar o pico de memória. Ajuste esse
comportamento diretamente na chamada:

```r
sby_adanear(
  sby_formula = alvo ~ .,
  sby_data = dados,
  sby_knn_query_chunk_size = 2000L
)
```

Para cálculo exato em matrizes densas de alta dimensionalidade, prefira a engine
`native` explícita. A rota de compatibilidade `FNN` com `sby_knn_algorithm =
"brute"` usa a mesma implementação nativa quando os kernels estão disponíveis e a
opção `sbyadanear.sby_use_native_brute` permanece ativa:

```r
sby_adanear(
  sby_formula = alvo ~ .,
  sby_data = dados,
  sby_knn_engine = "native",
  sby_knn_algorithm = "brute",
  sby_knn_distance_metric = "euclidean",
  sby_knn_workers = 1L,
  sby_knn_parallel_backend = "parallel"
)
```

Quando o R estiver ligado ao Intel oneAPI/MKL, `sby_knn_workers = 1L` permite que
o BLAS use seus próprios threads. Quando `sby_knn_workers > 1L`, o pacote reduz
threads BLAS por processo para evitar competição excessiva de CPU. Para trocar o
paralelismo por blocos do R por threads nativos no caminho exato bruto, combine
`sby_knn_engine = "native"`, `sby_knn_algorithm = "brute"` e
`sby_knn_parallel_backend = "RcppParallel"`.

`RcppParallel` decide o runtime concreto: em plataformas suportadas ele usa TBB
/ oneTBB e, nas demais, cai para TinyThread. Por isso o pacote não expõe um
parâmetro separado `oneTBB`; quando `sby_audit = TRUE`, os diagnósticos incluem
`sby_knn_parallel_runtime` para indicar se a execução efetiva foi
`"parallel"`, `"RcppParallel::TBB"` ou `"RcppParallel::TinyThread"`.

O kernel `RcppParallel` do `sbyadanear` não contém regiões OpenMP internas. Ainda
assim, duplo paralelismo pode ocorrer se o usuário envolver a chamada em outro
backend paralelo ou se uma biblioteca numérica externa abrir threads ao mesmo
tempo. Em servidores com Intel oneAPI/MKL, o pacote mitiga esse cenário
reduzindo temporariamente `OMP_NUM_THREADS` e `MKL_NUM_THREADS` para `1` quando
`sby_knn_workers > 1L`; em pipelines já paralelos, prefira
`sby_knn_workers = 1L` por tarefa externa.

Para HNSW com maior proximidade em relação ao resultado exato, aumente `M` e
`ef`. A configuração abaixo prioriza fidelidade sobre velocidade e memória:

```r
sby_adanear(
  sby_formula = alvo ~ .,
  sby_data = dados,
  sby_knn_engine = "RcppHNSW",
  sby_knn_algorithm = "auto",
  sby_knn_distance_metric = "euclidean",
  sby_knn_hnsw_m = 32L,
  sby_knn_hnsw_ef = 1000L,
  sby_knn_workers = parallel::detectCores(logical = FALSE)
)
```

Em bases muito sensíveis à vizinhança local, valores como `sby_knn_hnsw_m = 48L`
e `sby_knn_hnsw_ef = 2000L` podem aproximar mais a seleção do resultado exato,
com maior consumo de memória e tempo de construção do índice.

## Funções principais

```r
# Pipeline híbrido: primeiro gera amostras sintéticas com ADASYN e depois
# reduz a classe majoritária com NearMiss-1
sby_adanear(
  sby_formula,
  sby_data,
  sby_over_ratio = 0.2,
  sby_under_ratio = 1,
  sby_knn_over_k = 5L,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE
)

# Somente sobreamostragem ADASYN da classe minoritária
sby_adasyn(
  sby_formula,
  sby_data,
  sby_over_ratio = 0.2,
  sby_knn_over_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE
)

# Somente subamostragem NearMiss-1 da classe majoritária
sby_nearmiss(
  sby_formula,
  sby_data,
  sby_under_ratio = 1,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE
)
```

## Etapas para `recipes`

O pacote também oferece etapas supervisionadas para pipelines `recipes`:

```r
# Etapa ADASYN para recipes.
sby_step_adasyn(
  recipe,
  ...,
  sby_over_ratio = 0.2,
  sby_knn_over_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE
)

# Etapa NearMiss-1 para recipes.
sby_step_nearmiss(
  recipe,
  ...,
  sby_under_ratio = 1,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE
)

# Etapa combinada ADASYN + NearMiss-1 para recipes.
sby_step_adanear(
  recipe,
  ...,
  sby_over_ratio = 0.2,
  sby_under_ratio = 1,
  sby_knn_over_k = 5L,
  sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE
)
```

Por padrão, as etapas usam `skip = TRUE`, pois alteram o número de linhas do
conjunto processado e normalmente devem ser aplicadas apenas no treinamento

## Exemplo rápido com `sby_adanear()`

```r
library(sbyadanear)

# A semente aqui controla apenas a criação do exemplo reproduzível
set.seed(42)

# Cria dois preditores numéricos. As rotinas de sampling esperam preditores
# numéricos; variáveis categóricas devem ser tratadas antes do balanceamento
sby_x <- tibble::tibble(
  sby_a = rnorm(40),
  sby_b = rnorm(40)
)

# Cria um alvo binário desbalanceado: 10 observações minoritárias e 30
# majoritárias. A coluna do alvo pode ter qualquer nome na entrada
sby_y <- factor(c(rep("minority", 10), rep("majority", 30)))

# Junta preditores e alvo em um único data frame, pois a API pública usa fórmula + dados
sby_data <- tibble::add_column(sby_x, sby_y = sby_y)

# Aplica o pipeline híbrido:
# - sby_formula = sby_y ~ . informa que sby_y é o alvo e as demais colunas são
#   preditores;
# - sby_over_ratio controla a geração sintética ADASYN;
# - sby_under_ratio controla a razão minoria/maioria final no NearMiss-1;
#   use 1 para reduzir a majoritária até igualar a minoritária disponível;
# - sby_seed fixo torna a geração e desempates reproduzíveis.
sby_balanced <- sby_adanear(
  sby_formula = sby_y ~ .,
  sby_data = sby_data,
  sby_over_ratio = 0.5,
  sby_under_ratio = 0.8,
  sby_seed = 123
)

# O retorno padrão é um tibble. A coluna alvo é padronizada como TARGET.
sby_balanced
```

## Auditoria

```r
# Com sby_audit = TRUE, a função retorna uma lista com dados finais,
# resultados intermediários e diagnósticos de contagem/configuração
sby_audit <- sby_adanear(
  sby_formula = sby_y ~ .,
  sby_data = sby_data,
  sby_over_ratio = 0.5,
  sby_under_ratio = 0.8,
  sby_seed = 123,
  sby_audit = TRUE
)

# Diagnósticos incluem contagens de linhas, distribuição de classes e
# parâmetros KNN resolvidos
sby_audit$sby_diagnostics

# Dados balanceados finais.
sby_audit$sby_balanced_data
```

## Exemplos individuais

```r
# Apenas ADASYN: aumenta adaptativamente a classe minoritária e mantém todos os
# exemplos originais
sby_only_over <- sby_adasyn(
  sby_formula = sby_y ~ .,
  sby_data = sby_data,
  sby_over_ratio = 0.5,
  sby_seed = 123
)

# Apenas NearMiss-1: reduz a classe majoritária priorizando exemplos próximos à
# classe minoritária
sby_only_under <- sby_nearmiss(
  sby_formula = sby_y ~ .,
  sby_data = sby_data,
  sby_under_ratio = 0.8,
  sby_seed = 123
)
```

## Exemplo com `recipes`

```r
library(recipes)

# Define uma recipe simples. O desfecho é sby_y e os preditores são sby_a/sby_b.
sby_rec <- recipe(sby_y ~ ., data = sby_data)

# Seleciona explicitamente o desfecho para a etapa supervisionada
sby_rec <- sby_step_adanear(
  recipe = sby_rec,
  all_outcomes(),
  sby_over_ratio = 0.5,
  sby_under_ratio = 0.8,
  sby_seed = 123
)

# prep() treina a etapa e resolve a coluna de desfecho selecionada.
sby_rec_prepped <- prep(sby_rec, training = sby_data)

# bake() aplica a etapa ao conjunto informado. Como a etapa altera linhas, use
# com cuidado fora do treinamento.
sby_rec_balanced <- bake(sby_rec_prepped, new_data = sby_data)
```

## Instalação local

```sh
R CMD INSTALL .
```

## Dependências

Dependências importadas pelo pacote:

```r
install.packages(c(
  "Rcpp", "RcppHNSW", "RcppParallel", "FNN", "cli", "generics",
  "recipes", "rlang", "tibble", "collapse", "data.table", "kit",
  "Rfast", "coop"
))
```

Dependências opcionais usadas em testes, benchmarks ou engines opcionais:

```r
install.packages(c("modeldata", "Matrix", "KernelKnn", "bigKNN", "bigmemory", "testthat"))
```


## Ambiente de desenvolvimento

O repositório inclui um `Dockerfile` com R, toolchain de compilação e as
dependências de sistema necessárias para desenvolvimento e validação do pacote

```sh
docker build -t sbyadanear-r .
docker run --rm -it -v "$PWD":/workspace/sbyadanear sbyadanear-r
```

Dentro do container, valide o pacote com:

```sh
R CMD build .
R CMD check sbyadanear_0.3.0.tar.gz
```

## Arquivos principais

- `DESCRIPTION`: metadados, versão e dependências do pacote R.
- `NAMESPACE`: funções exportadas, métodos S3 e carregamento da biblioteca
  nativa.
- `R/`: funções R, helpers internos e métodos S3 das etapas `recipes`.
- `src/sbyadanear.c`: kernels nativos em C compilados na instalação do pacote.
- `man/`: documentação gerada a partir dos blocos roxygen2.

## Validação recomendada

```sh
R CMD build .
R CMD check sbyadanear_0.3.0.tar.gz
R CMD INSTALL .
```

Quando o binário do R não estiver disponível no ambiente, valide ao menos a
estrutura textual com `git diff --check` e buscas com `rg`.


## Ambiente R para validacao operacional completa

Para executar validacao completa localmente e no GitHub, o pacote requer um
ambiente com R, toolchain de compilacao C e dependencias opcionais para testes
especificos. O caminho recomendado e usar o Docker do proprio repositorio.

### Opcao 1: validacao local com Docker

```bash
docker build -f docker/oraclelinux97-r453/Dockerfile -t sbyadanear:oraclelinux97-r453 .
docker run --rm -it sbyadanear:oraclelinux97-r453 R --version
docker run --rm -it -v "$PWD":/workspace/r-package-validation sbyadanear:oraclelinux97-r453 Rscript tools/docker/run_many_tests.R
```

### Opcao 2: validacao automatizada no GitHub Actions

O workflow `.github/workflows/main.yml` ja executa:

- build da imagem Oracle Linux com R 4.5.3
- execucao repetida de `tools/docker/run_many_tests.R`
- upload dos artefatos CSV em `test-results`

### Sobre MKL no Docker

Sim, o Docker pode executar testes com MKL quando a imagem estiver configurada
com oneAPI/MKL e variaveis de ambiente adequadas. No pacote `sbyadanear`, MKL
nao e dependencia obrigatoria do pacote em si. O beneficio vem do ambiente R e
do backend BLAS/LAPACK configurado na imagem.

Para diagnostico dentro do container:

```r
Sys.getenv("OMP_NUM_THREADS")
Sys.getenv("MKL_NUM_THREADS")
```

Se `RhpcBLASctl` estiver instalado no ambiente:

```r
RhpcBLASctl::blas_get_num_procs()
RhpcBLASctl::omp_get_num_procs()
```

### Recomendacao de threads

Quando houver paralelismo externo no pipeline, limitar threads de BLAS/OpenMP
normalmente reduz oversubscription:

```r
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)
```

Quando o processo for unico e computacionalmente intenso, aumentar threads pode
ser util, dependendo do hardware e do backend numerico.
