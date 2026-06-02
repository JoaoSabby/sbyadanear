# Revisão de oportunidades de desempenho: kit, coop e collapse

> Nota de escopo: não localizei pacote CRAN chamado `kid` nas fontes consultadas. Como `kit`, `coop` e `collapse` compõem o ecossistema fastverse e `kit` é um pacote de rotinas C para manipulação de dados, esta revisão interpreta `kid` como `kit`. Se o alvo era outro pacote chamado `kid` fora do CRAN, esta análise deve ser refeita com a origem correta.

## Fontes consultadas

- `kit` 0.0.20: índice rdrr com 65 funções e 12 páginas de manual; descrição oficial: funções básicas em C para grandes manipulações de dados, `iif`/`nif`/`vswitch`, funções paralelas tipo `psum`/`pprod`, e símbolos chamáveis em C.
- `coop` 0.6-3: manual r-universe; descrição oficial: implementações rápidas e eficientes em memória de covariância, correlação e similaridade cosseno, com despacho S3; rotas densas baseadas em BLAS `dsyrk`/`dgemm`; suporte a matrizes densas e alguns formatos esparsos.
- `collapse` 2.1.7: manual r-universe e índice rdrr; descrição oficial: pacote C/C++ amplo, com estatísticas agrupadas/ponderadas, OpenMP, agrupamento/ordenação rápidos, transformação, reshape, séries temporais/painéis, conversões rápidas e preservação de classes.

## Critério de avaliação para o `sbyadanear`

- **Alto**: pode substituir caminho quente atual com expectativa clara de ganho em tempo/memória sem mudar semântica central.
- **Médio**: útil em cenários específicos ou como fallback opcional; exige benchmark antes de virar dependência.
- **Baixo**: benefício marginal porque o pacote já possui C nativo, BLAS, ou porque a função resolve problema periférico.
- **Não**: não se aplica ao fluxo ADASYN/NearMiss/KNN atual, ou aumenta dependência/risco sem ganho esperado.
- **Risco**: pode até acelerar caso isolado, mas cria oversubscription, alocação grande ou diferença semântica relevante.

## Síntese executiva

1. **Não recomendo adicionar `collapse` ou `kit` como `Imports` só por desempenho** neste momento. O `sbyadanear` já usa kernels C próprios para z-score, ADASYN, seleção NearMiss, top-k KNN bruto e `rbind` denso. A maior parte dos ganhos desses pacotes ficaria em bordas tabulares, não no miolo numérico.
2. **Recomendo avaliar `collapse::qM()` como `Suggests` opcional** para conversão data frame -> matrix em bases tabulares grandes, comparando contra `data.matrix()` dentro de `sby_adanear_as_numeric_matrix()`.
3. **Recomendo avaliar `kit::funique()` / `kit::uniqLen()` apenas se surgirem gargalos de validação/diagnóstico** com data frames grandes. A contagem binária de classes já usa `tabulate()` e provavelmente não ganha com `kit::count()`.
4. **Não recomendo `coop::cosine()` para o KNN atual**: ele é rápido para matriz completa de similaridade, mas materializa `n x n` ou `n_query x n_ref`; o `sbyadanear` precisa de top-k, então a alocação completa é o gargalo.
5. **Cuidado com duplo paralelismo**: `collapse` usa OpenMP em várias rotas; `coop` depende de BLAS; `sbyadanear` já tem BLAS/RcppParallel/RcppHNSW. Misturar sem controle pode piorar tempo por oversubscription.

## Oportunidades recomendadas para benchmark

| Prioridade | Ideia | Funções candidatas | Local provável no `sbyadanear` | Ganho esperado | Observação |
|---|---|---|---|---|---|
| P1 | Conversão rápida tabular -> matrix | `collapse::qM`, `collapse::qDF`, `collapse::qDT` | `R/sby_adanear_as_numeric_matrix.R` | Médio | Só como `Suggests`; validar nomes, fatores e tipos numéricos. |
| P1 | Auditoria/diagnóstico de duplicados e cardinalidade | `kit::funique`, `kit::fduplicated`, `kit::uniqLen`; `collapse::funique`, `collapse::fduplicated`, `collapse::fnunique` | validações e testes de dados | Baixo/Médio | Útil apenas se validações passarem a consumir tempo relevante. |
| P2 | Resumo tabular rápido em auditoria | `collapse::qsu`, `collapse::descr`, `collapse::fsummarise` | futuras rotinas de auditoria/relatórios | Médio | Não usar no caminho padrão para não aumentar dependência. |
| P2 | Ordenação/particionamento | `collapse::radixorder`, `roworder`, `kit::psort`, `kit::topn` | seleção/diagnóstico, não KNN bruto atual | Baixo | O top-k KNN já é heap C; `topn` não evita cálculo de distâncias. |
| P3 | Similaridade/correlação completa | `coop::cosine`, `coop::pcor`, `coop::covar` | exploração/diagnóstico, não produção KNN | Risco | Matriz completa pode estourar memória; melhor manter top-k streaming. |

## Checklist: `kit`

| Função | Benefício para `sbyadanear`? | Uso potencial / motivo |
|---|---:|---|
| `charToFact` | Baixo | Conversão de caracteres para fator; fluxos atuais validam fator alvo explicitamente. |
| `clearData` | Não | Gerência de dados compartilhados entre sessões; fora do fluxo. |
| `count` | Baixo | Contagens gerais; classe binária já usa `tabulate()` com custo mínimo. |
| `countNA` | Baixo | Validação de NA poderia usar, mas `anyNA()` já é simples e base. |
| `countOccur` | Baixo | Útil para auditoria, não caminho quente. |
| `fduplicated` | Médio | Pode acelerar validações de duplicatas em data frames grandes, se necessárias. |
| `fpos` | Não | Procura posição de matriz dentro de matriz maior; não há uso KNN/ADASYN. |
| `funique` | Médio | Potencial para auditoria/validação de linhas únicas; não melhora KNN. |
| `getData` | Não | Ligado a `shareData`; fora do pacote. |
| `iif` | Baixo | Pode substituir `ifelse()` vetorizado, mas há pouco uso crítico. |
| `nif` | Baixo | Nested if vetorizado; sem gargalo atual identificado. |
| `nswitch` | Baixo | Switch vetorizado; sem gargalo atual identificado. |
| `pall` | Baixo | Reduções paralelas por linha entre vetores; validações atuais são simples. |
| `pallNA` | Baixo | Detecção paralela de NA entre vetores; `anyNA()`/matriz C já resolvem. |
| `pallv` | Baixo | Comparações por valor entre vetores; uso periférico. |
| `pany` | Baixo | Similar a `pall`; não substitui kernel numérico. |
| `panyNA` | Baixo | Útil apenas para auditoria de NA multi-coluna. |
| `panyv` | Baixo | Uso periférico. |
| `pcount` | Baixo | Contagem por valor; `tabulate()` é suficiente para alvo binário. |
| `pcountNA` | Baixo | Validação de NA periférica. |
| `pfirst` | Não | Primeiro não ausente entre vetores; sem uso. |
| `plast` | Não | Último não ausente entre vetores; sem uso. |
| `pmean` | Baixo | Média paralela entre vetores, não por coluna; z-score já é C nativo. |
| `pprod` | Não | Produto paralelo entre vetores; sem uso. |
| `psort` | Baixo | Ordenação paralela; top-k KNN usa heap C e evita sort completo. |
| `psum` | Baixo | Soma entre vetores; z-score/ADASYN já usam C. |
| `setlevels` | Baixo | Alteração de níveis por referência; risco de efeitos colaterais no alvo. |
| `shareData` | Não | Compartilhamento entre sessões; fora de escopo. |
| `topn` | Médio | Pode ajudar em seleções top-n genéricas; não substitui top-k por distância sem calcular matriz completa. |
| `uniqLen` | Médio | Contagem rápida de únicos para auditoria; benefício se houver validação pesada. |
| `vswitch` | Baixo | Switch vetorizado; uso periférico. |

## Checklist: `coop`

| Função | Benefício para `sbyadanear`? | Uso potencial / motivo |
|---|---:|---|
| `cosine` | Risco | Acelera similaridade completa, mas materializa matriz completa; KNN precisa top-k. |
| `tcosine` | Risco | Mesma limitação de `cosine`, orientada a linhas. |
| `covar` | Não | Covariância completa não é necessária ao ADASYN/NearMiss. |
| `tcovar` | Não | Idem, por linhas. |
| `pcor` | Não | Correlação completa não é usada no fluxo. |
| `tpcor` | Não | Idem, por linhas. |
| `scaler` | Baixo | Poderia substituir scaling, mas `sby_compute_z_score_params()` e `sby_apply_z_score_scaling_matrix()` já são C nativo e auditáveis. |
| `sparsity` | Médio | Útil para diagnosticar sparse/densificação antes de abortar, mas não melhora desempenho de rota densa. |
| `cosine_wt` | Não | Similaridade ponderada completa; fora do algoritmo atual. |
| `covar_wt` | Não | Covariância ponderada completa; fora do algoritmo atual. |
| `pcor_wt` | Não | Correlação ponderada completa; fora do algoritmo atual. |

## Checklist: `collapse` por famílias públicas/documentadas

| Família / funções | Benefício para `sbyadanear`? | Uso potencial / motivo |
|---|---:|---|
| `across` | Baixo | Conveniência para múltiplas colunas; não substitui kernels nativos. |
| Operadores aritméticos `%r+%`, `%r-%`, `%r*%`, `%r/%`, `%rr%`, `%c+%`, `%c-%`, `%c*%`, `%c/%`, `%cr%` | Baixo | Poderiam simplificar operações row/col, mas KNN/z-score já são C/BLAS. |
| Operadores in-place `%+=%`, `%-=%`, `%*=%`, `%/=%`, `%=%`, `%==%`, `%!=%` | Risco | Mutação por referência pode dificultar garantias de cópia/atributos. |
| Matching `%in%`, `%!in%`, `%iin%`, `%!iin%`, `fmatch`, `ckmatch` | Médio | Potencial para validações e índices; ganho pequeno frente a KNN. |
| Split-apply `BY` | Baixo | Útil em relatórios agrupados; não no caminho quente atual. |
| Agregação `collap`, `collapv`, `collapg` | Médio | Boa opção para auditorias agrupadas futuras; não necessário na geração sintética. |
| `colorder`, `colorderv` | Baixo | Conveniência de ordenação de colunas; sem impacto. |
| `dapply` | Baixo | Aplicação em data frames; não substitui C nativo. |
| `descr` | Médio | Auditoria descritiva rápida opcional. |
| Helpers `alloc`, `copyv`, `setv`, `setop`, `massign`, `allNA`, `allv`, `anyv`, `whichNA`, `whichv`, `missing_cases`, `na_rm`, `na_omit`, `na_insert`, `na_locf`, `na_focb`, `cinv`, `vec`, `vgcd`, `vlengths`, `vtypes`, `seq_row`, `seq_col` | Baixo | Úteis para programação/memória; ganho só em validações periféricas. |
| Conversões `qM`, `qDF`, `qDT`, `qTBL`, `qF`, `qG`, `as_integer_factor`, `as_numeric_factor`, `as_character_factor`, `as_factor_GRP`, `as_factor_qG` | Médio | `qM` é o principal candidato para benchmark em entrada tabular. |
| Seleção/variáveis `num_vars`, `cat_vars`, `char_vars`, `date_vars`, `fact_vars`, `logi_vars`, `get_vars`, `add_vars`, `av`, `nv`, `slt`, `gvr` | Baixo | Pode acelerar helpers de seleção; recipes/base já cuidam disso. |
| Renomeação/labels `frename`, `rnm`, `relabel`, `setrename`, `setrelabel`, `setLabels`, `vlabels`, `namlab`, `rm_stub`, `add_stub` | Baixo | Sem gargalo; cuidado com atributos. |
| Agrupamento `GRP`, `GRPN`, `GRPid`, `GRPnames`, `group`, `groupv`, `gby`, `fgroup_by`, `fgroup_vars`, `groupid`, `seqid`, `timeid`, `finteraction`, `qG`, `itn`, `is_GRP`, `is_qG` | Médio | Poderia acelerar agrupamentos/auditoria, mas classes são binárias e já baratas. |
| Ordenação `radixorder`, `radixorderv`, `roworder`, `roworderv`, `greorder` | Baixo/Médio | Útil para ordenações grandes; top-k KNN não deve ordenar tudo. |
| Estatísticas rápidas `fsum`, `fmean`, `fmedian`, `fnth`, `fmode`, `fmin`, `fmax`, `fprod`, `fvar`, `fsd`, `fquantile`, `frange`, `fnobs`, `fndistinct`, `fnunique`, `fNobs`, `fNdistinct` | Baixo/Médio | Fortes genericamente, mas z-score/contagens já têm C nativo; úteis para auditoria. |
| Transformações estatísticas `fbetween`, `fwithin`, `fscale`, `standardize`, `B`, `W`, `STD`, `TRA`, `setTRA` | Baixo | Sobrepõe scaling/centering existentes; benchmark só se remover C próprio. |
| Cumulativas/diferenças/lag `fcumsum`, `fdiff`, `D`, `Dlog`, `fgrowth`, `G`, `flag`, `L`, `F` | Não | Séries temporais/painel; fora do algoritmo. |
| Modelagem `flm`, `fFtest` | Não | Fora do balanceamento. |
| Distância `fdist` | Risco | Distâncias completas podem alocar demais; top-k streaming é preferível. |
| Uniques/duplicatas `funique`, `fduplicated`, `any_duplicated` | Médio | Útil para validações de duplicatas/cardinalidade. |
| Contagem agrupada `fcount`, `fcountv`, `qtab`, `qtable` | Baixo/Médio | Auditoria/distribuição de classes; `tabulate()` já basta para binário. |
| Data manipulation `fselect`, `fslice`, `fslicev`, `fsubset`, `ss`, `sbt`, `ftransform`, `ftransformv`, `settransform`, `settransformv`, `tfm`, `tfmv`, `fcompute`, `fcomputev`, `fmutate`, `fsummarise`, `fsummarize`, `smr`, `fungroup` | Baixo/Médio | Pode acelerar wrappers tabulares; não deve entrar no núcleo numérico sem benchmark. |
| Joins/reshape/listas `join`, `pivot`, `rowbind`, `unlist2d`, `rbindlist`, `rsplit`, `gsplit`, `rapply2d`, `t_list`, `ldepth`, `list_elem`, `atomic_elem`, `reg_elem`, `irreg_elem`, `has_elem`, `get_elem`, `extract_list` | Baixo | Útil em infraestrutura de dados, não na computação KNN. |
| Painel/séries `indexing`, `findex`, `findex_by`, `ix`, `iby`, `to_plm`, `unindex`, `reindex`, `psmat`, `psacf`, `pspacf`, `psccf`, `tscomp`, `varying` | Não | Fora de escopo. |
| Correlação/covariância `pwcor`, `pwcov`, `pwnobs`, `pwNobs` | Não/Risco | Mesma lógica de `coop`: matriz completa não substitui top-k. |
| Recodificação/outliers `recode_char`, `recode_num`, `replace_NA`, `replace_Inf`, `replace_na`, `replace_inf`, `replace_outliers` | Baixo | Pré-processamento de usuário; pacote deve validar, não imputar/recodificar silenciosamente. |
| Fatores/níveis `fdroplevels`, `fnlevels`, `is_categorical`, `is_date`, `is_irregular`, `is_unlistable`, `is.Date`, `is.GRP`, `is.qG`, `is.categorical`, `is.unlistable` | Baixo | Útil em validações periféricas. |
| Opções/documentação `set_collapse`, `get_collapse`, `collapse`, `collapse-package`, `collapse-documentation`, `collapse-options`, `collapse-renamed`, `A*` tópicos | Não | Metadados/documentação. |
| Dados `GGDC10S`, `wlddev` | Não | Datasets. |

## Checklist detalhado de nomes revisados

### `kit`

- [x] `charToFact` — Baixo
- [x] `clearData` — Não
- [x] `count` — Baixo
- [x] `countNA` — Baixo
- [x] `countOccur` — Baixo
- [x] `fduplicated` — Médio
- [x] `fpos` — Não
- [x] `funique` — Médio
- [x] `getData` — Não
- [x] `iif` — Baixo
- [x] `nif` — Baixo
- [x] `nswitch` — Baixo
- [x] `pall` — Baixo
- [x] `pallNA` — Baixo
- [x] `pallv` — Baixo
- [x] `pany` — Baixo
- [x] `panyNA` — Baixo
- [x] `panyv` — Baixo
- [x] `pcount` — Baixo
- [x] `pcountNA` — Baixo
- [x] `pfirst` — Não
- [x] `plast` — Não
- [x] `pmean` — Baixo
- [x] `pprod` — Não
- [x] `psort` — Baixo
- [x] `psum` — Baixo
- [x] `setlevels` — Baixo/Risco
- [x] `shareData` — Não
- [x] `topn` — Médio
- [x] `uniqLen` — Médio
- [x] `vswitch` — Baixo

### `coop`

- [x] `cosine` — Risco
- [x] `tcosine` — Risco
- [x] `covar` — Não
- [x] `tcovar` — Não
- [x] `pcor` — Não
- [x] `tpcor` — Não
- [x] `scaler` — Baixo
- [x] `sparsity` — Médio
- [x] `cosine_wt` — Não
- [x] `covar_wt` — Não
- [x] `pcor_wt` — Não

### `collapse`

- [x] `%!=%`, `%!iin%`, `%!in%`, `%*=%`, `%+=%`, `%-=%`, `%/=%`, `%=%`, `%==%`, `%c*%`, `%c+%`, `%c-%`, `%c/%`, `%cr%`, `%iin%`, `%r*%`, `%r+%`, `%r-%`, `%r/%`, `%rr%` — Baixo/Risco conforme operador
- [x] `B`, `W`, `STD`, `TRA`, `setTRA`, `fbetween`, `fwithin`, `fscale`, `standardize` — Baixo
- [x] `BY` — Baixo
- [x] `D`, `Dlog`, `F`, `G`, `L`, `flag`, `fdiff`, `fgrowth`, `fcumsum` — Não/Baixo
- [x] `GRP`, `GRPN`, `GRPid`, `GRPnames`, `group`, `groupv`, `gby`, `fgroup_by`, `fgroup_vars`, `groupid`, `seqid`, `timeid`, `qG`, `qF`, `itn`, `is_GRP`, `is_qG` — Médio/Baixo
- [x] `HDB`, `HDW`, `fHDbetween`, `fHDwithin`, `fhdbetween`, `fhdwithin` — Não/Baixo
- [x] `across` — Baixo
- [x] `add_stub`, `rm_stub`, `add_vars`, `av`, `cat_vars`, `char_vars`, `date_vars`, `fact_vars`, `logi_vars`, `num_vars`, `nv`, `slt`, `gvr`, `get_vars` — Baixo
- [x] `allNA`, `all_funs`, `all_identical`, `all_obj_equal`, `alloc`, `allv`, `anyv`, `copyAttrib`, `copyMostAttrib`, `copyv`, `massign`, `setop`, `setv`, `seq_col`, `seq_row`, `vec`, `vgcd`, `vlengths`, `vtypes`, `whichNA`, `whichv` — Baixo
- [x] `any_duplicated`, `fduplicated`, `funique`, `fnunique`, `fnlevels`, `uniq`-equivalentes — Médio
- [x] `as_character_factor`, `as_integer_factor`, `as_numeric_factor`, `as_factor_GRP`, `as_factor_qG`, `qDF`, `qDT`, `qM`, `qTBL` — Médio para `qM`; demais Baixo
- [x] `atomic_elem`, `atomic_elem<-`, `list_elem`, `list_elem<-`, `reg_elem`, `irreg_elem`, `has_elem`, `get_elem`, `extract_list`, `ldepth`, `t_list`, `rapply2d` — Baixo
- [x] `cinv`, `ckmatch`, `fmatch`, `%in%` rápidos — Médio/Baixo
- [x] `collap`, `collapv`, `collapg` — Médio para auditoria; Baixo no núcleo
- [x] `colorder`, `colorderv`, `roworder`, `roworderv`, `radixorder`, `radixorderv`, `greorder` — Baixo/Médio
- [x] `dapply`, `fcompute`, `fcomputev`, `ftransform`, `ftransformv`, `settransform`, `settransformv`, `tfm`, `tfmv`, `fselect`, `fslice`, `fslicev`, `fsubset`, `ss`, `sbt`, `fmutate`, `fsummarise`, `fsummarize`, `smr`, `fungroup` — Baixo/Médio
- [x] `descr`, `qsu` — Médio para auditoria
- [x] `fFtest`, `flm` — Não
- [x] `fNdistinct`, `fNobs`, `fcount`, `fcountv`, `qtab`, `qtable`, `fnobs`, `fndistinct` — Baixo/Médio
- [x] `fdim`, `fnrow`, `fncol`, `fNCOL` — Baixo
- [x] `fdist` — Risco
- [x] `fdroplevels`, `is_categorical`, `is_date`, `is_irregular`, `is_unlistable`, `is.Date`, `is.GRP`, `is.qG`, `is.categorical`, `is.unlistable` — Baixo
- [x] `ffirst`, `flast`, `fmin`, `fmax`, `fsum`, `fmean`, `fmedian`, `fnth`, `fmode`, `fprod`, `fquantile`, `frange`, `fvar`, `fsd` — Baixo/Médio; fortes em geral, mas redundantes com C nativo atual
- [x] `frename`, `rnm`, `relabel`, `setrename`, `setrelabel`, `setLabels`, `vlabels`, `vlabels<-`, `namlab`, `setAttrib`, `setattrib`, `setColnames`, `setDimnames`, `setRownames`, `unattrib` — Baixo
- [x] `indexing`, `ix`, `iby`, `to_plm`, `unindex`, `reindex`, `findex`, `findex_by`, `psmat`, `psacf`, `pspacf`, `psccf`, `tscomp`, `varying` — Não
- [x] `join`, `pivot`, `rowbind`, `unlist2d`, `rsplit`, `gsplit` — Baixo
- [x] `mctl`, `mrtl` — Baixo
- [x] `missing_cases`, `na_focb`, `na_insert`, `na_locf`, `na_omit`, `na_rm` — Baixo
- [x] `pad` — Não/Baixo
- [x] `pwcor`, `pwcov`, `pwnobs`, `pwNobs` — Não/Risco
- [x] `recode_char`, `recode_num`, `replace_Inf`, `replace_NA`, `replace_inf`, `replace_na`, `replace_outliers` — Baixo/Risco semântica
- [x] `set_collapse`, `get_collapse` — Não
- [x] `wlddev`, `GGDC10S` — Não

## Conclusão prática

Para uma próxima PR de implementação, eu faria apenas benchmarks pequenos e isolados:

1. `collapse::qM()` versus `data.matrix()` em `sby_adanear_as_numeric_matrix()` com data frames numéricos grandes e mistos.
2. `kit::uniqLen()`/`kit::funique()` versus base em validações que realmente apareçam em perfis.
3. `coop::cosine()` apenas como experimento de diagnóstico em matrizes pequenas; não como backend KNN padrão.

Sem benchmark demonstrando ganho material, a recomendação é **não adicionar dependências novas**: o custo de instalação, o risco de conflito de threads e a superfície de manutenção superam os ganhos prováveis no caminho crítico atual.
