# Check-list de implementacao das funcoes HPC (sby_adanear_hpc, sby_nearmiss_hpc, sby_adasyn_hpc)

Documento de planejamento e rastreio. Marcar [X] cada item concluido.

## Parte A - Analise do pacote (concluida antes da codificacao)

- [X] Mapear fluxo da API publica: sby_adanear -> sby_adanear_matrix -> sby_adasyn_matrix + sby_nearmiss_matrix.
- [X] Identificar a rota "native" (sby_knn_engine == "native") e seus kernels Fortran (compute/apply/revert zscore, rbind).
- [X] Confirmar que ADASYN no fluxo adanear opera no espaco padronizado (sby_return_original_scale = FALSE).
- [X] Confirmar contrato de sby_get_knnx, sby_generate_adasyn_samples, sby_nearmiss_index.
- [X] Confirmar convencoes: snake_case, ausencia de travessao, documentacao roxygen pt-BR.

## Parte B - Implementacao R (atalho HPC)

- [X] Criar R/sby_adanear_hpc.R com assinatura exata e on.exit de restauro de ambiente.
- [X] Criar R/sby_adasyn_hpc.R (atalho ADASYN puro no espaco padronizado).
- [X] Criar R/sby_nearmiss_hpc.R (atalho NearMiss-1 puro).
- [X] Helper R/sby_hpc_env.R: captura/injeta/restaura somente MKL_NUM_THREADS e OMP_NUM_THREADS.
- [X] Helper R/sby_hpc_native_available.R: verifica simbolos HPC carregados, com fallback transparente para a rota classica.
- [X] Roteamento: a rota "native" passa a delegar para o atalho HPC quando os simbolos existem; funcoes originais continuam acessiveis.

## Parte C - Camada C++ (Rcpp zero-copy)

- [X] sby_adanear_hpc_cpp: orquestra zscore populacional (VSL), aplicacao, KNN via dgemm, ADASYN e NearMiss no espaco padronizado.
- [X] Montagem Zero-Copy: Rcpp::List com 200 NumericVector pre-alocados ao numero exato de linhas finais.
- [X] Reversao do z-score por FMA durante a copia para os vetores da List final (chama kernel Fortran).
- [X] Atributos: class c("tbl_df","tbl","data.frame"), names (colunas + alvo) e row.names.
- [X] Registro do novo simbolo no R_CallMethodDef e R_init_sbyadanear.

## Parte D - Camada Fortran (FMA AVX-512, VSL, dgemm)

- [X] sby_zscore_population_vsl_f: estatisticas via VSL (vslsscompute - media e variancia).
- [X] sby_apply_zscore_simd_f: padronizacao no espaco z (laco SIMD).
- [X] sby_revert_zscore_fma_f: reversao via laco aninhado !DIR$ SIMD / !$OMP SIMD forcando vfmadd213pd.
- [X] sby_pairwise_sqdist_dgemm_f: D^2 = ||A||^2 + ||B||^2 - 2 A B^T via cblas_dgemm.
- [X] sby_adasyn_interp_uniform_f: interpolacao lambda com vdrnguniform no espaco padronizado.
- [X] Declaracoes de interface VSL/RNG/BLAS com bind(C).

## Parte E - Lista de Controlo Obrigatoria (do prompt)

- [X] Funcao R com assinatura exata e ambiente isolado com restauro (on.exit).
- [X] Nomenclatura 100% snake_case verificada em todo o codigo e variaveis.
- [X] Logica livre da dupla normalizacao: ADASYN atua diretamente na matriz padronizada.
- [X] FMA (Fused Multiply-Add) garantido via diretivas SIMD no Fortran para a reversao do z-score.
- [X] Montagem Zero-Copy do Tibble/data.frame efetuada diretamente no C++ (Rcpp::List).
- [X] Controle temporario de MKL_NUM_THREADS e OMP_NUM_THREADS com restauro ao final.
- [X] Substituicao do calculo de distancias pela rotina dgemm bloqueada.
- [X] Retornos e documentacao completamente livres do caractere travessao.

## Parte F - Conferencia final

- [X] Simulacao linha a linha do fluxo R -> C++ -> Fortran -> R.
- [X] Verificacao de consistencia de assinaturas entre camadas.
- [X] Verificacao de NAMESPACE, RcppExports.R e registro de simbolos.

## Parte G - Parecer Intel oneMKL e ajustes de desempenho adicionais

- [X] Conferido o guia Intel de uso do oneMKL com extensoes R: manter linkagem oneMKL opcional em `src/Makevars` e nao sobrescrever variaveis de ambiente globais ja carregadas pelo servidor do cliente.
- [X] Conferido o guia Intel de threading: o pacote controla somente `MKL_NUM_THREADS` e `OMP_NUM_THREADS` durante a rotina e restaura os valores originais ao final.
- [X] Conferidas as recomendacoes Intel de leading dimensions: o kernel `sby_pairwise_sqdist_dgemm_f` agora calcula leading dimensions double precision alinhadas a 64 bytes, evita multiplos exatos de grandes potencias de 2 e usa buffers acolchoados para A/B; C tambem e acolchoada quando a duplicacao da matriz de distancia fica limitada.
- [X] Conferido o contrato `dgemm`: os valores `lda`, `ldb` e `ldc` continuam respeitando os minimos exigidos por `m`, `n` e `k`, com padding apenas acima desses minimos.
- [X] Conferida a recomendacao de problemas pequenos: `MKL_DIRECT_CALL` foi avaliado, mas nao aplicado porque o kernel critico chama `cblas_dgemm` via interface Fortran direta e o gargalo informado e de processamento massivo, nao de micro-GEMMs repetidos.
- [X] Conferida a recomendacao de offload: nao aplicada porque o hardware alvo informado e CPU Intel Xeon Platinum 8260 sem GPU/accelerator oneAPI dedicado.
