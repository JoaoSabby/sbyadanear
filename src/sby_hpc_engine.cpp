/*
 * sby_hpc_engine.cpp
 *
 * Motor HPC consolidado do pacote sbyadanear. Este wrapper orquestra o
 * pipeline ADASYN + NearMiss-1 estritamente no espaco padronizado por z-score
 * e monta o tibble final por zero-copy diretamente em C++, sem passar por
 * as.data.frame() ou rbind() da camada R.
 *
 * Decisoes de desempenho:
 *   - estatisticas iniciais da populacao via Vector Statistics Library (Fortran).
 *   - matriz de distancias por D^2 = ||A||^2 + ||B||^2 - 2 A B^T com cblas_dgemm.
 *   - interpolacao lambda do ADASYN com vdrnguniform no espaco padronizado.
 *   - reversao do z-score por FMA AVX-512 durante a copia para os vetores finais.
 *
 * Toda a nomenclatura segue snake_case. Nenhum caractere de travessao e usado.
 *
 * Autor: Joao Batista Goncalves de Brito
 */

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

// Interfaces dos kernels Fortran do motor HPC
extern "C" {
  void sby_zscore_population_vsl_f(const double *x, int n, int p,
                                   double *means, double *sds, int *status);
  void sby_apply_zscore_simd_f(const double *x, int n, int p,
                               const double *means, const double *sds,
                               double *x_out, int *status);
  void sby_revert_zscore_fma_f(const double *x, int n, int p,
                               const double *means, const double *sds,
                               double *x_out, int *status);
  void sby_pairwise_sqdist_dgemm_f(const double *a, int n_a,
                                   const double *b, int n_b, int p,
                                   double *d_out, int *status);
  void sby_adasyn_interp_uniform_f(const double *minority, int n_min, int p,
                                   const int *base_idx, const int *nbr_idx,
                                   const double *lambda, int n_syn,
                                   double *syn_out, int *status);
}

// -------------------------------------------------------------------
// sby_resolve_minority_role
// Determina o codigo do nivel minoritario de um fator binario com base nas
// contagens. Retorna o codigo 1-based do nivel menos frequente.
// -------------------------------------------------------------------
static int sby_resolve_minority_role(const Rcpp::IntegerVector& y_codes, int n_levels){
  std::vector<long> counts(n_levels + 1, 0L);
  for(R_xlen_t i = 0; i < y_codes.size(); ++i){
    int code = y_codes[i];
    if(code >= 1 && code <= n_levels){
      counts[code] += 1L;
    }
  }
  int minority_code = 1;
  long minority_count = -1L;
  for(int c = 1; c <= n_levels; ++c){
    if(minority_count < 0L || counts[c] < minority_count){
      minority_count = counts[c];
      minority_code = c;
    }
  }
  return minority_code;
}

// -------------------------------------------------------------------
// sby_zscore_population
// Computa centros e escalas populacionais chamando o kernel Fortran VSL.
// -------------------------------------------------------------------
static void sby_zscore_population(const double* x, int n, int p,
                                  std::vector<double>& means,
                                  std::vector<double>& sds){
  means.assign(p, 0.0);
  sds.assign(p, 0.0);
  int status = 0;
  sby_zscore_population_vsl_f(x, n, p, means.data(), sds.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha no calculo de z-score populacional (status=%d)", status);
  }
  // Protege contra escalas nulas para nao destruir a reversao posterior
  for(int j = 0; j < p; ++j){
    if(!(sds[j] > 0.0) || !std::isfinite(sds[j])){
      sds[j] = 1.0;
    }
  }
}

// -------------------------------------------------------------------
// sby_apply_zscore
// Aplica o z-score por kernel Fortran SIMD sobre a matriz n x p.
// -------------------------------------------------------------------
static void sby_apply_zscore(const double* x, int n, int p,
                             const std::vector<double>& means,
                             const std::vector<double>& sds,
                             std::vector<double>& x_scaled){
  x_scaled.assign((size_t) n * (size_t) p, 0.0);
  int status = 0;
  sby_apply_zscore_simd_f(x, n, p, means.data(), sds.data(), x_scaled.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha na aplicacao de z-score (status=%d)", status);
  }
}

// -------------------------------------------------------------------
// sby_knn_topk_against_reference
// Calcula a matriz de distancias ao quadrado de uma matriz de consulta contra
// uma matriz de referencia usando cblas_dgemm e devolve os k vizinhos mais
// proximos (indices 1-based) por linha de consulta. Ambas as matrizes estao em
// layout column major n x p no espaco padronizado.
// -------------------------------------------------------------------
static void sby_knn_topk_against_reference(
    const std::vector<double>& query, int n_query,
    const std::vector<double>& reference, int n_ref, int p,
    int k, bool drop_self,
    std::vector< std::vector<int> >& out_index){
  std::vector<double> dist2((size_t) n_query * (size_t) n_ref, 0.0);
  int status = 0;
  sby_pairwise_sqdist_dgemm_f(query.data(), n_query, reference.data(), n_ref, p,
                              dist2.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha no calculo de distancias por dgemm (status=%d)", status);
  }

  out_index.assign(n_query, std::vector<int>());
  int effective_k = k;
  if(effective_k > n_ref){
    effective_k = n_ref;
  }

  for(int i = 0; i < n_query; ++i){
    std::vector<int> order(n_ref);
    std::iota(order.begin(), order.end(), 0);
    const double* drow = dist2.data() + (size_t) i * (size_t) n_ref;
    int take = effective_k;
    if(take < (int) order.size()){
      std::partial_sort(order.begin(), order.begin() + take, order.end(),
        [&](int a, int b){
          if(drow[a] != drow[b]) return drow[a] < drow[b];
          return a < b;
        });
    } else {
      std::sort(order.begin(), order.end(),
        [&](int a, int b){
          if(drow[a] != drow[b]) return drow[a] < drow[b];
          return a < b;
        });
    }

    std::vector<int>& dst = out_index[i];
    dst.reserve(take);
    for(int t = 0; t < take && (int) dst.size() < take; ++t){
      int cand = order[t];
      if(drop_self && cand == i){
        continue;
      }
      dst.push_back(cand + 1); // 1-based
    }
    // Garante pelo menos um vizinho quando drop_self removeu o unico candidato
    if(dst.empty() && n_ref > 0){
      dst.push_back((order[0] == i && n_ref > 1) ? (order[1] + 1) : (order[0] + 1));
    }
  }
}

// -------------------------------------------------------------------
// sby_assemble_tibble_zero_copy
// Monta o tibble final diretamente em C++. Recebe um buffer consolidado no
// espaco padronizado (linhas finais x p, column major), os parametros de
// z-score, os codigos do alvo e os metadados de nomes. Aloca uma Rcpp::List
// com p + 1 colunas, cada uma um NumericVector pre-alocado ao numero exato de
// linhas finais, e reverte o z-score por FMA durante a copia.
// -------------------------------------------------------------------
static SEXP sby_assemble_tibble_zero_copy(
    const std::vector<double>& final_scaled, int n_final, int p,
    const std::vector<double>& means, const std::vector<double>& sds,
    const std::vector<int>& y_codes_final,
    const Rcpp::CharacterVector& column_names,
    const std::string& target_name,
    const Rcpp::CharacterVector& target_levels){

  // Reversao do z-score por FMA diretamente no buffer column major final.
  std::vector<double> final_original((size_t) n_final * (size_t) p, 0.0);
  int status = 0;
  sby_revert_zscore_fma_f(final_scaled.data(), n_final, p,
                          means.data(), sds.data(), final_original.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha na reversao de z-score por FMA (status=%d)", status);
  }

  // Aloca a lista com tamanho exato de colunas (p preditores + 1 alvo).
  Rcpp::List out(p + 1);

  // Cada coluna preditora e um NumericVector pre-alocado ao numero exato de
  // linhas finais. A copia e column major contigua, portanto zero-copy logico.
  for(int j = 0; j < p; ++j){
    Rcpp::NumericVector col(n_final);
    const double* src = final_original.data() + (size_t) j * (size_t) n_final;
    std::copy(src, src + n_final, col.begin());
    out[j] = col;
  }

  // Coluna do alvo reconstruida como fator preservando os niveis originais.
  Rcpp::IntegerVector target(n_final);
  for(int i = 0; i < n_final; ++i){
    target[i] = y_codes_final[i];
  }
  target.attr("levels") = target_levels;
  target.attr("class") = "factor";
  out[p] = target;

  // Nomes das colunas: preditores seguidos do nome do alvo.
  Rcpp::CharacterVector out_names(p + 1);
  for(int j = 0; j < p; ++j){
    out_names[j] = column_names[j];
  }
  out_names[p] = target_name;
  out.attr("names") = out_names;

  // Atributo de classe que entrega o objeto pronto como tibble a camada R.
  out.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");

  // row.names compacto no formato interno do R.
  Rcpp::IntegerVector row_names = Rcpp::IntegerVector::create(NA_INTEGER, -n_final);
  out.attr("row.names") = row_names;

  return out;
}

// -------------------------------------------------------------------
// sby_run_adasyn_stage
// Executa o ADASYN no espaco padronizado. Anexa as linhas sinteticas ao final
// do buffer escalado e atualiza os codigos do alvo. Devolve o numero de linhas
// apos a expansao.
// -------------------------------------------------------------------
static int sby_run_adasyn_stage(
    std::vector<double>& x_scaled, int n, int p,
    std::vector<int>& y_codes, int minority_code,
    int k_neighbor, double over_ratio){

  // Indices minoritarios 1-based no conjunto completo
  std::vector<int> minority_index;
  minority_index.reserve(n);
  for(int i = 0; i < n; ++i){
    if(y_codes[i] == minority_code){
      minority_index.push_back(i);
    }
  }
  int n_min = (int) minority_index.size();
  if(n_min < 2){
    Rcpp::stop("ADASYN exige ao menos 2 observacoes minoritarias");
  }

  // Matriz minoritaria em layout column major n_min x p
  std::vector<double> minority((size_t) n_min * (size_t) p, 0.0);
  for(int j = 0; j < p; ++j){
    const double* col = x_scaled.data() + (size_t) j * (size_t) n;
    double* dst = minority.data() + (size_t) j * (size_t) n_min;
    for(int r = 0; r < n_min; ++r){
      dst[r] = col[minority_index[r]];
    }
  }

  // Quantidade sintetica seguindo a regra floor(min_count * over_ratio), min 1
  int synthetic_count = (int) std::floor((double) n_min * over_ratio);
  if(synthetic_count < 1){
    synthetic_count = 1;
  }

  // Vizinhos minoritarios para interpolacao (dgemm + topk, removendo self)
  int effective_k = k_neighbor;
  if(effective_k > n_min - 1){
    effective_k = n_min - 1;
  }
  if(effective_k < 1){
    effective_k = 1;
  }
  std::vector< std::vector<int> > minority_neighbors;
  sby_knn_topk_against_reference(minority, n_min, minority, n_min, p,
                                 effective_k, true, minority_neighbors);

  // Distribuicao uniforme das linhas sinteticas entre as bases minoritarias
  std::vector<int> base_idx(synthetic_count);
  std::vector<int> nbr_idx(synthetic_count);
  std::vector<double> lambda(synthetic_count);

  Rcpp::NumericVector unif_lambda = Rcpp::runif(synthetic_count);
  Rcpp::NumericVector unif_pick = Rcpp::runif(synthetic_count);
  for(int s = 0; s < synthetic_count; ++s){
    int base = s % n_min;            // 0-based base minoritaria
    base_idx[s] = base + 1;          // 1-based para o Fortran
    const std::vector<int>& nbrs = minority_neighbors[base];
    int pick = (int) std::floor(unif_pick[s] * (double) nbrs.size());
    if(pick < 0) pick = 0;
    if(pick >= (int) nbrs.size()) pick = (int) nbrs.size() - 1;
    nbr_idx[s] = nbrs[pick];         // ja 1-based
    lambda[s] = unif_lambda[s];
  }

  // Geracao sintetica por interpolacao FMA no espaco padronizado
  std::vector<double> syn_out((size_t) synthetic_count * (size_t) p, 0.0);
  int status = 0;
  sby_adasyn_interp_uniform_f(minority.data(), n_min, p,
                              base_idx.data(), nbr_idx.data(), lambda.data(),
                              synthetic_count, syn_out.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha na interpolacao sintetica ADASYN (status=%d)", status);
  }

  // Anexa as linhas sinteticas ao buffer escalado, preservando column major.
  int n_new = n + synthetic_count;
  std::vector<double> expanded((size_t) n_new * (size_t) p, 0.0);
  for(int j = 0; j < p; ++j){
    const double* src_col = x_scaled.data() + (size_t) j * (size_t) n;
    const double* syn_col = syn_out.data() + (size_t) j * (size_t) synthetic_count;
    double* dst_col = expanded.data() + (size_t) j * (size_t) n_new;
    std::copy(src_col, src_col + n, dst_col);
    std::copy(syn_col, syn_col + synthetic_count, dst_col + n);
  }
  x_scaled.swap(expanded);

  y_codes.reserve(n_new);
  for(int s = 0; s < synthetic_count; ++s){
    y_codes.push_back(minority_code);
  }
  return n_new;
}

// -------------------------------------------------------------------
// sby_run_nearmiss_stage
// Executa o NearMiss-1 no espaco padronizado. Retem todas as linhas
// minoritarias e seleciona as linhas majoritarias com menor distancia media
// aos k vizinhos minoritarios mais proximos. Devolve os indices retidos
// 0-based ordenados.
// -------------------------------------------------------------------
static std::vector<int> sby_run_nearmiss_stage(
    const std::vector<double>& x_scaled, int n, int p,
    const std::vector<int>& y_codes, int minority_code,
    int k_neighbor, double under_ratio){

  std::vector<int> minority_index;
  std::vector<int> majority_index;
  minority_index.reserve(n);
  majority_index.reserve(n);
  for(int i = 0; i < n; ++i){
    if(y_codes[i] == minority_code){
      minority_index.push_back(i);
    } else {
      majority_index.push_back(i);
    }
  }
  int n_min = (int) minority_index.size();
  int n_maj = (int) majority_index.size();

  // Quantidade majoritaria retida: floor(n_min / under_ratio), limitada por n_maj
  int retained_majority = (int) std::floor((double) n_min / under_ratio);
  if(retained_majority > n_maj){
    retained_majority = n_maj;
  }
  if(retained_majority < 0){
    retained_majority = 0;
  }

  // Matrizes column major das classes
  std::vector<double> majority((size_t) n_maj * (size_t) p, 0.0);
  std::vector<double> minority((size_t) n_min * (size_t) p, 0.0);
  for(int j = 0; j < p; ++j){
    const double* col = x_scaled.data() + (size_t) j * (size_t) n;
    double* mdst = majority.data() + (size_t) j * (size_t) n_maj;
    double* ndst = minority.data() + (size_t) j * (size_t) n_min;
    for(int r = 0; r < n_maj; ++r){ mdst[r] = col[majority_index[r]]; }
    for(int r = 0; r < n_min; ++r){ ndst[r] = col[minority_index[r]]; }
  }

  // Distancias majoritaria contra minoritaria por dgemm
  std::vector<double> dist2((size_t) n_maj * (size_t) n_min, 0.0);
  int status = 0;
  sby_pairwise_sqdist_dgemm_f(majority.data(), n_maj, minority.data(), n_min, p,
                              dist2.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha no calculo de distancias NearMiss por dgemm (status=%d)", status);
  }

  int effective_k = k_neighbor;
  if(effective_k > n_min){
    effective_k = n_min;
  }
  if(effective_k < 1){
    effective_k = 1;
  }

  // Score NearMiss-1: media dos k menores quadrados de distancia minoritaria
  std::vector< std::pair<double,int> > scores(n_maj);
  for(int m = 0; m < n_maj; ++m){
    const double* drow = dist2.data() + (size_t) m * (size_t) n_min;
    std::vector<double> row(drow, drow + n_min);
    std::partial_sort(row.begin(), row.begin() + effective_k, row.end());
    double acc = 0.0;
    for(int t = 0; t < effective_k; ++t){
      acc += row[t];
    }
    scores[m] = std::make_pair(acc / (double) effective_k, majority_index[m]);
  }
  std::sort(scores.begin(), scores.end(),
    [](const std::pair<double,int>& a, const std::pair<double,int>& b){
      if(a.first != b.first) return a.first < b.first;
      return a.second < b.second;
    });

  std::vector<int> retained;
  retained.reserve(n_min + retained_majority);
  for(int i = 0; i < n_min; ++i){
    retained.push_back(minority_index[i]);
  }
  for(int t = 0; t < retained_majority; ++t){
    retained.push_back(scores[t].second);
  }
  std::sort(retained.begin(), retained.end());
  return retained;
}

// -------------------------------------------------------------------
// sby_extract_factor_codes
// Extrai os codigos inteiros 1-based de um fator R.
// -------------------------------------------------------------------
static Rcpp::IntegerVector sby_extract_factor_codes(SEXP y){
  Rcpp::IntegerVector codes(y);
  return codes;
}

//' @title Motor HPC consolidado do pipeline ADASYN mais NearMiss-1
//' @description Executa todo o pipeline no espaco padronizado e monta o tibble por zero-copy.
// [[Rcpp::export]]
extern "C" SEXP sby_adanear_hpc_cpp(SEXP x_matrix, SEXP y_factor,
                                    SEXP k_adanear, SEXP k_nearmiss,
                                    SEXP over_ratio, SEXP under_ratio,
                                    SEXP max_threads, SEXP column_names,
                                    SEXP target_name, SEXP target_levels){
  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow();
  int p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels = levels.size();
  int k_over = Rcpp::as<int>(k_adanear);
  int k_under = Rcpp::as<int>(k_nearmiss);
  double ratio_over = Rcpp::as<double>(over_ratio);
  double ratio_under = Rcpp::as<double>(under_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  // 1. Estatisticas populacionais e padronizacao (espaco z unico para tudo)
  std::vector<double> means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  std::vector<double> x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());

  // 2. ADASYN no espaco padronizado (expande o buffer e os codigos do alvo)
  int n_after_over = sby_run_adasyn_stage(x_scaled, n, p, y_codes,
                                          minority_code, k_over, ratio_over);

  // 3. NearMiss-1 no espaco padronizado sobre o conjunto expandido
  std::vector<int> retained = sby_run_nearmiss_stage(x_scaled, n_after_over, p,
                                                     y_codes, minority_code,
                                                     k_under, ratio_under);
  int n_final = (int) retained.size();

  // 4. Consolida o buffer final column major e os codigos finais do alvo
  std::vector<double> final_scaled((size_t) n_final * (size_t) p, 0.0);
  std::vector<int> y_final(n_final);
  for(int j = 0; j < p; ++j){
    const double* src = x_scaled.data() + (size_t) j * (size_t) n_after_over;
    double* dst = final_scaled.data() + (size_t) j * (size_t) n_final;
    for(int r = 0; r < n_final; ++r){
      dst[r] = src[retained[r]];
    }
  }
  for(int r = 0; r < n_final; ++r){
    y_final[r] = y_codes[retained[r]];
  }

  // 5. Montagem zero-copy do tibble com reversao FMA do z-score
  return sby_assemble_tibble_zero_copy(final_scaled, n_final, p, means, sds,
                                       y_final, Rcpp::CharacterVector(column_names),
                                       Rcpp::as<std::string>(target_name), levels);
}

//' @title Motor HPC consolidado do pipeline ADASYN mais NearMiss-1 com metadados
//' @description Executa o pipeline no espaco padronizado e retorna matriz final escalada, indices retidos, alvo e parametros de z-score.
// [[Rcpp::export]]
extern "C" SEXP sby_adanear_hpc_result_cpp(SEXP x_matrix, SEXP y_factor,
                                           SEXP k_adanear, SEXP k_nearmiss,
                                           SEXP over_ratio, SEXP under_ratio,
                                           SEXP max_threads, SEXP column_names,
                                           SEXP target_levels){
  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow();
  int p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels = levels.size();
  int k_over = Rcpp::as<int>(k_adanear);
  int k_under = Rcpp::as<int>(k_nearmiss);
  double ratio_over = Rcpp::as<double>(over_ratio);
  double ratio_under = Rcpp::as<double>(under_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  std::vector<double> means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  std::vector<double> x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());
  int n_after_over = sby_run_adasyn_stage(x_scaled, n, p, y_codes,
                                          minority_code, k_over, ratio_over);

  std::vector<int> retained = sby_run_nearmiss_stage(x_scaled, n_after_over, p,
                                                     y_codes, minority_code,
                                                     k_under, ratio_under);
  int n_final = (int) retained.size();

  Rcpp::NumericMatrix final_scaled(n_final, p);
  for(int j = 0; j < p; ++j){
    const double* src = x_scaled.data() + (size_t) j * (size_t) n_after_over;
    for(int r = 0; r < n_final; ++r){
      final_scaled(r, j) = src[retained[r]];
    }
  }
  final_scaled.attr("dimnames") = Rcpp::List::create(R_NilValue, column_names);

  Rcpp::IntegerVector y_final(n_final);
  Rcpp::IntegerVector retained_index(n_final);
  for(int r = 0; r < n_final; ++r){
    y_final[r] = y_codes[retained[r]];
    retained_index[r] = retained[r] + 1;
  }
  y_final.attr("levels") = levels;
  y_final.attr("class") = "factor";

  Rcpp::NumericVector centers(p);
  Rcpp::NumericVector scales(p);
  for(int j = 0; j < p; ++j){
    centers[j] = means[j];
    scales[j] = sds[j];
  }

  return Rcpp::List::create(
    Rcpp::Named("sby_final_scaled") = final_scaled,
    Rcpp::Named("sby_y_vector") = y_final,
    Rcpp::Named("sby_retained_index") = retained_index,
    Rcpp::Named("sby_scaling_info") = Rcpp::List::create(
      Rcpp::Named("centers") = centers,
      Rcpp::Named("scales") = scales
    )
  );
}

//' @title Motor HPC consolidado do oversampling ADASYN
//' @description Executa apenas o ADASYN no espaco padronizado e monta o tibble por zero-copy.
// [[Rcpp::export]]
extern "C" SEXP sby_adasyn_hpc_cpp(SEXP x_matrix, SEXP y_factor,
                                   SEXP k_adanear, SEXP over_ratio,
                                   SEXP max_threads, SEXP column_names,
                                   SEXP target_name, SEXP target_levels){
  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow();
  int p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels = levels.size();
  int k_over = Rcpp::as<int>(k_adanear);
  double ratio = Rcpp::as<double>(over_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  std::vector<double> means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  std::vector<double> x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());
  int n_after_over = sby_run_adasyn_stage(x_scaled, n, p, y_codes,
                                          minority_code, k_over, ratio);

  return sby_assemble_tibble_zero_copy(x_scaled, n_after_over, p, means, sds,
                                       y_codes, Rcpp::CharacterVector(column_names),
                                       Rcpp::as<std::string>(target_name), levels);
}

//' @title Motor HPC consolidado do undersampling NearMiss-1
//' @description Executa apenas o NearMiss-1 no espaco padronizado e monta o tibble por zero-copy.
// [[Rcpp::export]]
extern "C" SEXP sby_nearmiss_hpc_cpp(SEXP x_matrix, SEXP y_factor,
                                     SEXP k_nearmiss, SEXP under_ratio,
                                     SEXP max_threads, SEXP column_names,
                                     SEXP target_name, SEXP target_levels){
  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow();
  int p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels = levels.size();
  int k_under = Rcpp::as<int>(k_nearmiss);
  double ratio = Rcpp::as<double>(under_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  std::vector<double> means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  std::vector<double> x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());
  std::vector<int> retained = sby_run_nearmiss_stage(x_scaled, n, p, y_codes,
                                                     minority_code, k_under, ratio);
  int n_final = (int) retained.size();

  std::vector<double> final_scaled((size_t) n_final * (size_t) p, 0.0);
  std::vector<int> y_final(n_final);
  for(int j = 0; j < p; ++j){
    const double* src = x_scaled.data() + (size_t) j * (size_t) n;
    double* dst = final_scaled.data() + (size_t) j * (size_t) n_final;
    for(int r = 0; r < n_final; ++r){
      dst[r] = src[retained[r]];
    }
  }
  for(int r = 0; r < n_final; ++r){
    y_final[r] = y_codes[retained[r]];
  }

  return sby_assemble_tibble_zero_copy(final_scaled, n_final, p, means, sds,
                                       y_final, Rcpp::CharacterVector(column_names),
                                       Rcpp::as<std::string>(target_name), levels);
}
