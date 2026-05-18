# instenginer 0.4.0 (em desenvolvimento)

## Mudancas de comportamento (breaking)

* `sby_knn_engine = "auto"` agora seleciona o engine de forma sensivel ao
  contexto, em vez de retornar sempre `"FNN"`:
  - Metricas nao euclidianas (`"cosine"`, `"ip"`): seleciona `RcppHNSW`
    automaticamente (o `FNN` nao suporta essas metricas neste pacote).
  - Bases grandes em alta dimensionalidade
    (`n * p >= 5e6` AND `p >= 50` por padrao): seleciona `RcppHNSW` por
    custo de busca exata. Limite configuravel via
    `options(instenginer.sby_auto_engine_hnsw_min_cells = ...)`.
  - Demais casos: continua selecionando `FNN` exato.
  Para preservar exatamente o comportamento anterior em qualquer caso,
  passe `sby_knn_engine = "FNN"` explicitamente.

## Correcoes

* `sby_adasyn_matrix()`, `sby_nearmiss_matrix()` e `sby_nearmiss_index()`
  agora rejeitam classes minoritarias com menos de duas observacoes,
  alinhando o contrato da API de matriz com o da API tabular. Antes, a API
  de matriz aceitava `n_minority = 1` e gerava amostras sinteticas
  invalidas (interpolacao do ponto consigo mesmo).
* Teste `native NearMiss selector matches R fallback without ties` deixou
  de usar `.Call("OU_SelectNearMissMajorityC", ..., PACKAGE = "instenginer")`
  (forma proibida quando `R_useDynamicSymbols(FALSE)`); agora chama o
  simbolo nativo registrado pelo namespace.
* Mensagem mais clara quando uma coluna preditora se chama literalmente
  `TARGET` (a coluna `TARGET` e reservada para a saida).

## Desempenho

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
  `options(instenginer.sby_use_native_brute = FALSE)`.
* Variante experimental do kernel ADASYN com escrita column-friendly:
  `OU_GenerateSyntheticAdasynColC`, disponivel via
  `options(instenginer.sby_adasyn_kernel = "col")`. Em testes empiricos
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

# instenginer 0.3.0

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
