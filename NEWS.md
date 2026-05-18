# instenginer 0.3.0

## Novidades

* API matricial pública para ADASYN, NearMiss-1 e ADANEAR: `sby_adasyn_matrix()`, `sby_nearmiss_matrix()` e `sby_adanear_matrix()`.
* Novo roteador industrial `sby_balance_matrix()` para estratégias `none`, `weight`, `adasyn`, `nearmiss`, `adanear` e `adanearWeight`.
* `sby_nearmiss_index()` expõe seleção por índices para fluxos que precisam preservar matrizes externas, inclusive esparsas fora do pacote.
* Seleção NearMiss parcial em C com critério determinístico por menor distância média e menor índice em empates.
* Auditoria leve para diagnósticos resumidos sem matrizes intermediárias pesadas.
* Melhor controle de memória: `none` e `weight` não densificam entradas, KNN em blocos aloca apenas componentes solicitados e NearMiss evita recálculo de z-score.
* Wrappers tabulares preservam linhas originais e desnormalizam/restauram tipos apenas nas linhas sintéticas do ADASYN.
