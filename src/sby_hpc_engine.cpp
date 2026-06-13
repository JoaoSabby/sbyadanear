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
 *   - matriz de distancias por D^2 = ||A||^2 + ||B||^2 - 2 A B^T com cblas_sgemm.
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
#include <string>
#include <limits>
#include <cstdlib>
#include <new>
#include <type_traits>
#include <utility>
#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


// -------------------------------------------------------------------
// Alocador alinhado a 64 bytes para buffers numericos quentes.
// A documentacao Intel oneMKL recomenda alinhamento em fronteiras de 64 bytes
// para melhor desempenho em kernels vetorizados e chamadas BLAS.
// -------------------------------------------------------------------
template <typename T, std::size_t Alignment>
class sby_aligned_allocator {
public:
  using value_type = T;

  sby_aligned_allocator() noexcept = default;
  template <class U>
  sby_aligned_allocator(const sby_aligned_allocator<U, Alignment>&) noexcept {}

  T* allocate(std::size_t n){
    if(n > std::numeric_limits<std::size_t>::max() / sizeof(T)){
      throw std::bad_array_new_length();
    }
    if(n == 0){
      return nullptr;
    }
    void* ptr = nullptr;
#if defined(_MSC_VER)
    ptr = _aligned_malloc(n * sizeof(T), Alignment);
    if(ptr == nullptr){
      throw std::bad_alloc();
    }
#else
    if(posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0){
      throw std::bad_alloc();
    }
#endif
    return static_cast<T*>(ptr);
  }

  void deallocate(T* p, std::size_t) noexcept{
#if defined(_MSC_VER)
    _aligned_free(p);
#else
    free(p);
#endif
  }

  template <class U>
  struct rebind { using other = sby_aligned_allocator<U, Alignment>; };

  template <class U, class... Args>
  void construct(U* p, Args&&... args){
    ::new((void*) p) U(std::forward<Args>(args)...);
  }

  template <class U>
  void construct(U* p){
    if constexpr (std::is_trivially_default_constructible<U>::value){
      ::new((void*) p) U;
    } else {
      ::new((void*) p) U();
    }
  }
};

template <class T, class U, std::size_t Alignment>
bool operator==(const sby_aligned_allocator<T, Alignment>&,
                const sby_aligned_allocator<U, Alignment>&) noexcept { return true; }

template <class T, class U, std::size_t Alignment>
bool operator!=(const sby_aligned_allocator<T, Alignment>&,
                const sby_aligned_allocator<U, Alignment>&) noexcept { return false; }

using sby_double_buffer = std::vector<double, sby_aligned_allocator<double, 64> >;
using sby_float_buffer = std::vector<float, sby_aligned_allocator<float, 64> >;

template <typename Buffer>
static void sby_resize_first_touch(Buffer& buffer, size_t n, typename Buffer::value_type value){
  buffer.resize(n);
  typename Buffer::value_type* data = buffer.data();
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static)
#endif
  for(size_t i = 0; i < n; ++i){
    data[i] = value;
  }
}

template <typename T>
static void sby_parallel_copy(const T* src, T* dst, size_t n){
#ifdef _OPENMP
#pragma omp simd
#endif
  for(size_t i = 0; i < n; ++i){
    dst[i] = src[i];
  }
}


static int sby_resolve_gemm_ref_block(int n_query, int n_ref, int p){
  const size_t target_bytes = (size_t) 256 * (size_t) 1024 * (size_t) 1024;
  size_t row_count = (size_t) std::max(1, n_query);
  size_t by_dist = target_bytes / (sizeof(float) * row_count);
  if(by_dist < 1){
    by_dist = 1;
  }
  size_t by_copy = target_bytes / (sizeof(float) * (size_t) std::max(1, p));
  if(by_copy < 1){
    by_copy = 1;
  }
  size_t block = std::min(by_dist, by_copy);
  if(block > (size_t) n_ref){
    block = (size_t) n_ref;
  }
  if(block < 1){
    block = 1;
  }
  return (int) block;
}

static void sby_copy_column_block(const sby_float_buffer& source, int n_source,
                                  int p, int col_start, int n_block,
                                  sby_float_buffer& block){
  sby_resize_first_touch(block, (size_t) n_block * (size_t) p, 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* src = source.data() + (size_t) j * (size_t) n_source + (size_t) col_start;
    float* dst = block.data() + (size_t) j * (size_t) n_block;
    sby_parallel_copy(src, dst, (size_t) n_block);
  }
}

static int sby_resolve_native_threads(SEXP max_threads){
  int requested = Rf_asInteger(max_threads);
  if(requested == NA_INTEGER || requested < 1){
    return -1;
  }
#ifdef _OPENMP
  int max_available = omp_get_num_procs();
  if(max_available < 1){
    max_available = omp_get_max_threads();
  }
  if(max_available > 0 && requested > max_available){
    requested = max_available;
  }
#endif
  return requested;
}

static void sby_apply_native_threads(SEXP max_threads){
#ifdef _OPENMP
  int requested = sby_resolve_native_threads(max_threads);
  if(requested > 0){
    omp_set_num_threads(requested);
  }
#else
  (void) max_threads;
#endif
}

// Interfaces dos kernels Fortran do motor HPC
extern "C" {
  void sby_zscore_population_vsl_f(const double *x, int n, int p,
                                   double *means, double *sds, int *status);
  void sby_apply_zscore_simd_f(const double *x, int n, int p,
                               const double *means, const double *sds,
                               float *x_out, int *status);
  void sby_revert_zscore_fma_f(const float *x, int n, int p,
                               const double *means, const double *sds,
                               double *x_out, int *status);
  void sby_pairwise_sqdist_sgemm_f(const float *a, int n_a,
                                   const float *b, int n_b, int p,
                                   float *d_out, int *status);
  void sby_adasyn_interp_uniform_f(const float *minority, int n_min, int p,
                                   const int *base_idx, const int *nbr_idx,
                                   const float *lambda, int n_syn,
                                   float *syn_out, int *status);
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
                                  sby_double_buffer& means,
                                  sby_double_buffer& sds){
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
                             const sby_double_buffer& means,
                             const sby_double_buffer& sds,
                             sby_float_buffer& x_scaled){
  sby_resize_first_touch(x_scaled, (size_t) n * (size_t) p, 0.0f);
  int status = 0;
  sby_apply_zscore_simd_f(x, n, p, means.data(), sds.data(), x_scaled.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha na aplicacao de z-score (status=%d)", status);
  }
}

// -------------------------------------------------------------------
// sby_knn_topk_against_reference
// Calcula a matriz de distancias ao quadrado de uma matriz de consulta contra
// uma matriz de referencia usando cblas_sgemm e devolve os k vizinhos mais
// proximos (indices 1-based) por linha de consulta. Ambas as matrizes estao em
// layout column major n x p no espaco padronizado.
// -------------------------------------------------------------------
static void sby_knn_topk_against_reference(
    const sby_float_buffer& query, int n_query,
    const sby_float_buffer& reference, int n_ref, int p,
    int k, bool drop_self,
    std::vector< std::vector<int> >& out_index){
  out_index.assign(n_query, std::vector<int>());
  int effective_k = k;
  if(effective_k > n_ref){
    effective_k = n_ref;
  }
  if(effective_k < 1 || n_query < 1 || n_ref < 1){
    return;
  }

  int keep_k = effective_k + (drop_self ? 1 : 0);
  if(keep_k > n_ref){
    keep_k = n_ref;
  }

  sby_double_buffer top_dist;
  sby_resize_first_touch(top_dist, (size_t) n_query * (size_t) keep_k,
                         std::numeric_limits<double>::infinity());
  std::vector<int> top_index((size_t) n_query * (size_t) keep_k, -1);
  std::vector<int> worst_pos(n_query, 0);
  sby_double_buffer worst_val;
  sby_resize_first_touch(worst_val, (size_t) n_query,
                         std::numeric_limits<double>::infinity());

  int ref_block_size = sby_resolve_gemm_ref_block(n_query, n_ref, p);
  sby_float_buffer reference_block;
  sby_float_buffer dist_block;
  int status = 0;

  for(int ref_start = 0; ref_start < n_ref; ref_start += ref_block_size){
    int n_block = std::min(ref_block_size, n_ref - ref_start);
    sby_copy_column_block(reference, n_ref, p, ref_start, n_block, reference_block);
    sby_resize_first_touch(dist_block, (size_t) n_query * (size_t) n_block, 0.0f);
    sby_pairwise_sqdist_sgemm_f(query.data(), n_query, reference_block.data(), n_block, p,
                                dist_block.data(), &status);
    if(status != 0){
      Rcpp::stop("Falha no calculo blocado de distancias por sgemm (status=%d)", status);
    }

    for(int b = 0; b < n_block; ++b){
      int ref_idx = ref_start + b;
      const float* col = dist_block.data() + (size_t) b * (size_t) n_query;
      for(int i = 0; i < n_query; ++i){
        float candidate = col[i];
        if(candidate >= worst_val[i]){
          continue;
        }
        size_t row_offset = (size_t) i * (size_t) keep_k;
        int pos = worst_pos[i];
        top_dist[row_offset + (size_t) pos] = (double) candidate;
        top_index[row_offset + (size_t) pos] = ref_idx;

        int new_worst_pos = 0;
        double new_worst_val = top_dist[row_offset];
        for(int t = 1; t < keep_k; ++t){
          double val = top_dist[row_offset + (size_t) t];
          if(val > new_worst_val){
            new_worst_val = val;
            new_worst_pos = t;
          }
        }
        worst_pos[i] = new_worst_pos;
        worst_val[i] = new_worst_val;
      }
    }
  }

  for(int i = 0; i < n_query; ++i){
    std::vector<int> order(keep_k);
    std::iota(order.begin(), order.end(), 0);
    size_t row_offset = (size_t) i * (size_t) keep_k;
    std::sort(order.begin(), order.end(), [&](int a, int b){
      double da = top_dist[row_offset + (size_t) a];
      double db = top_dist[row_offset + (size_t) b];
      int ia = top_index[row_offset + (size_t) a];
      int ib = top_index[row_offset + (size_t) b];
      if(da != db) return da < db;
      return ia < ib;
    });

    std::vector<int>& dst = out_index[i];
    dst.reserve(effective_k);
    for(int pos : order){
      int cand = top_index[row_offset + (size_t) pos];
      if(cand < 0){
        continue;
      }
      if(drop_self && cand == i){
        continue;
      }
      dst.push_back(cand + 1);
      if((int) dst.size() >= effective_k){
        break;
      }
    }
    if(dst.empty() && n_ref > 0){
      dst.push_back(1);
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
    const sby_float_buffer& final_scaled, int n_final, int p,
    const sby_double_buffer& means, const sby_double_buffer& sds,
    const std::vector<int>& y_codes_final,
    const Rcpp::CharacterVector& column_names,
    const std::string& target_name,
    const Rcpp::CharacterVector& target_levels){

  // Reversao do z-score por FMA diretamente no buffer column major final.
  sby_double_buffer final_original;
  sby_resize_first_touch(final_original, (size_t) n_final * (size_t) p, 0.0);
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
  std::vector<double*> column_ptrs(p);
  for(int j = 0; j < p; ++j){
    Rcpp::NumericVector col(n_final);
    column_ptrs[j] = col.begin();
    out[j] = col;
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const double* src = final_original.data() + (size_t) j * (size_t) n_final;
    sby_parallel_copy(src, column_ptrs[j], (size_t) n_final);
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
    sby_float_buffer& x_scaled, int n, int p,
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
  sby_float_buffer minority;
  sby_resize_first_touch(minority, (size_t) n_min * (size_t) p, 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* col = x_scaled.data() + (size_t) j * (size_t) n;
    float* dst = minority.data() + (size_t) j * (size_t) n_min;
    for(int r = 0; r < n_min; ++r){
      dst[r] = col[minority_index[r]];
    }
  }

  // Quantidade sintetica seguindo a regra floor(min_count * over_ratio), min 1
  int synthetic_count = (int) std::floor((double) n_min * over_ratio);
  if(synthetic_count < 1){
    synthetic_count = 1;
  }

  // Vizinhos minoritarios para interpolacao (sgemm + topk, removendo self)
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
  sby_float_buffer lambda;
  lambda.resize((size_t) synthetic_count);

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
    lambda[s] = (float) unif_lambda[s];
  }

  // Geracao sintetica por interpolacao FMA no espaco padronizado
  sby_float_buffer syn_out;
  sby_resize_first_touch(syn_out, (size_t) synthetic_count * (size_t) p, 0.0f);
  int status = 0;
  sby_adasyn_interp_uniform_f(minority.data(), n_min, p,
                              base_idx.data(), nbr_idx.data(), lambda.data(),
                              synthetic_count, syn_out.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha na interpolacao sintetica ADASYN (status=%d)", status);
  }

  // Anexa as linhas sinteticas ao buffer escalado, preservando column major.
  int n_new = n + synthetic_count;
  sby_float_buffer expanded;
  sby_resize_first_touch(expanded, (size_t) n_new * (size_t) p, 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* src_col = x_scaled.data() + (size_t) j * (size_t) n;
    const float* syn_col = syn_out.data() + (size_t) j * (size_t) synthetic_count;
    float* dst_col = expanded.data() + (size_t) j * (size_t) n_new;
    sby_parallel_copy(src_col, dst_col, (size_t) n);
    sby_parallel_copy(syn_col, dst_col + n, (size_t) synthetic_count);
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
    const sby_float_buffer& x_scaled, int n, int p,
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
  sby_float_buffer majority;
  sby_float_buffer minority;
  sby_resize_first_touch(majority, (size_t) n_maj * (size_t) p, 0.0f);
  sby_resize_first_touch(minority, (size_t) n_min * (size_t) p, 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* col = x_scaled.data() + (size_t) j * (size_t) n;
    float* mdst = majority.data() + (size_t) j * (size_t) n_maj;
    float* ndst = minority.data() + (size_t) j * (size_t) n_min;
    for(int r = 0; r < n_maj; ++r){ mdst[r] = col[majority_index[r]]; }
    for(int r = 0; r < n_min; ++r){ ndst[r] = col[minority_index[r]]; }
  }

  int effective_k = k_neighbor;
  if(effective_k > n_min){
    effective_k = n_min;
  }
  if(effective_k < 1){
    effective_k = 1;
  }

  // Score NearMiss-1 em streaming: processa blocos da minoria por SGEMM e
  // atualiza os k menores valores por linha majoritaria sem materializar a
  // matriz completa n_maj x n_min.
  std::vector< std::pair<double,int> > scores(n_maj);
  sby_double_buffer score_acc;
  sby_resize_first_touch(score_acc, (size_t) n_maj, 0.0);

  sby_double_buffer topk;
  sby_double_buffer worst_val;
  std::vector<int> worst_pos;
  if(effective_k < n_min){
    sby_resize_first_touch(topk, (size_t) n_maj * (size_t) effective_k,
                           std::numeric_limits<double>::infinity());
    worst_pos.assign(n_maj, 0);
    sby_resize_first_touch(worst_val, (size_t) n_maj,
                           std::numeric_limits<double>::infinity());
  }

  int minority_block_size = sby_resolve_gemm_ref_block(n_maj, n_min, p);
  sby_float_buffer minority_block;
  sby_float_buffer dist_block;
  int status = 0;

  for(int min_start = 0; min_start < n_min; min_start += minority_block_size){
    int n_block = std::min(minority_block_size, n_min - min_start);
    sby_copy_column_block(minority, n_min, p, min_start, n_block, minority_block);
    sby_resize_first_touch(dist_block, (size_t) n_maj * (size_t) n_block, 0.0f);
    sby_pairwise_sqdist_sgemm_f(majority.data(), n_maj, minority_block.data(), n_block, p,
                                dist_block.data(), &status);
    if(status != 0){
      Rcpp::stop("Falha no calculo blocado de distancias NearMiss por sgemm (status=%d)", status);
    }

    if(effective_k == n_min){
      for(int c = 0; c < n_block; ++c){
        const float* col = dist_block.data() + (size_t) c * (size_t) n_maj;
#ifdef _OPENMP
#pragma omp simd
#endif
        for(int m = 0; m < n_maj; ++m){
          score_acc[m] += col[m];
        }
      }
    } else {
      for(int c = 0; c < n_block; ++c){
        const float* col = dist_block.data() + (size_t) c * (size_t) n_maj;
        for(int m = 0; m < n_maj; ++m){
          float candidate = col[m];
          if(candidate >= worst_val[m]){
            continue;
          }
          double* row_topk = topk.data() + (size_t) m * (size_t) effective_k;
          row_topk[worst_pos[m]] = candidate;

          int new_worst_pos = 0;
          double new_worst_val = row_topk[0];
          for(int t = 1; t < effective_k; ++t){
            if(row_topk[t] > new_worst_val){
              new_worst_val = row_topk[t];
              new_worst_pos = t;
            }
          }
          worst_pos[m] = new_worst_pos;
          worst_val[m] = new_worst_val;
        }
      }
    }
  }

  if(effective_k < n_min){
    for(int m = 0; m < n_maj; ++m){
      const double* row_topk = topk.data() + (size_t) m * (size_t) effective_k;
      double acc = 0.0;
#ifdef _OPENMP
#pragma omp simd reduction(+:acc)
#endif
      for(int t = 0; t < effective_k; ++t){
        acc += row_topk[t];
      }
      score_acc[m] = acc;
    }
  }

  for(int m = 0; m < n_maj; ++m){
    scores[m] = std::make_pair(score_acc[m] / (double) effective_k, majority_index[m]);
  }

  auto score_less = [](const std::pair<double,int>& a, const std::pair<double,int>& b){
    if(a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  };
  if(retained_majority < n_maj){
    std::partial_sort(scores.begin(), scores.begin() + retained_majority, scores.end(), score_less);
  } else {
    std::sort(scores.begin(), scores.end(), score_less);
  }

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

//' @title Relatorio de compilacao do motor HPC
//' @description Retorna macros de compilacao usadas para validar AVX-512, AVX2, FMA e OpenMP.
// [[Rcpp::export]]
extern "C" SEXP sby_hpc_compile_report_cpp(){
  Rcpp::List out;

#if defined(SBYADANEAR_CASCADELAKE_NATIVE)
  out["cascade_lake_native"] = true;
#else
  out["cascade_lake_native"] = false;
#endif

#if defined(__AVX512F__)
  out["avx512f"] = true;
#else
  out["avx512f"] = false;
#endif

#if defined(__AVX512CD__)
  out["avx512cd"] = true;
#else
  out["avx512cd"] = false;
#endif

#if defined(__AVX512BW__)
  out["avx512bw"] = true;
#else
  out["avx512bw"] = false;
#endif

#if defined(__AVX512DQ__)
  out["avx512dq"] = true;
#else
  out["avx512dq"] = false;
#endif

#if defined(__AVX512VL__)
  out["avx512vl"] = true;
#else
  out["avx512vl"] = false;
#endif

#if defined(__AVX2__)
  out["avx2"] = true;
#else
  out["avx2"] = false;
#endif

#if defined(__FMA__)
  out["fma"] = true;
#else
  out["fma"] = false;
#endif

#ifdef _OPENMP
  out["openmp"] = true;
  out["openmp_version"] = _OPENMP;
  out["openmp_max_threads"] = omp_get_max_threads();
#else
  out["openmp"] = false;
  out["openmp_version"] = NA_INTEGER;
  out["openmp_max_threads"] = NA_INTEGER;
#endif

  return out;
}

//' @title Motor HPC consolidado do pipeline ADASYN mais NearMiss-1
//' @description Executa todo o pipeline no espaco padronizado e monta o tibble por zero-copy.
// [[Rcpp::export]]
extern "C" SEXP sby_adanear_hpc_cpp(SEXP x_matrix, SEXP y_factor,
                                    SEXP k_adanear, SEXP k_nearmiss,
                                    SEXP over_ratio, SEXP under_ratio,
                                    SEXP max_threads, SEXP column_names,
                                    SEXP target_name, SEXP target_levels){
  sby_apply_native_threads(max_threads);

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
  sby_double_buffer means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  sby_float_buffer x_scaled;
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
  sby_float_buffer final_scaled;
  sby_resize_first_touch(final_scaled, (size_t) n_final * (size_t) p, 0.0f);
  std::vector<int> y_final(n_final);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* src = x_scaled.data() + (size_t) j * (size_t) n_after_over;
    float* dst = final_scaled.data() + (size_t) j * (size_t) n_final;
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
  sby_apply_native_threads(max_threads);

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

  sby_double_buffer means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  sby_float_buffer x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());
  int n_after_over = sby_run_adasyn_stage(x_scaled, n, p, y_codes,
                                          minority_code, k_over, ratio_over);

  std::vector<int> retained = sby_run_nearmiss_stage(x_scaled, n_after_over, p,
                                                     y_codes, minority_code,
                                                     k_under, ratio_under);
  int n_final = (int) retained.size();

  Rcpp::NumericMatrix final_scaled(n_final, p);
  double* final_scaled_ptr = final_scaled.begin();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* src = x_scaled.data() + (size_t) j * (size_t) n_after_over;
    double* dst = final_scaled_ptr + (size_t) j * (size_t) n_final;
    for(int r = 0; r < n_final; ++r){
      dst[r] = (double) src[retained[r]];
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
  sby_apply_native_threads(max_threads);

  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow();
  int p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels = levels.size();
  int k_over = Rcpp::as<int>(k_adanear);
  double ratio = Rcpp::as<double>(over_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  sby_double_buffer means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  sby_float_buffer x_scaled;
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
  sby_apply_native_threads(max_threads);

  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow();
  int p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels = levels.size();
  int k_under = Rcpp::as<int>(k_nearmiss);
  double ratio = Rcpp::as<double>(under_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  sby_double_buffer means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  sby_float_buffer x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());
  std::vector<int> retained = sby_run_nearmiss_stage(x_scaled, n, p, y_codes,
                                                     minority_code, k_under, ratio);
  int n_final = (int) retained.size();

  sby_float_buffer final_scaled;
  sby_resize_first_touch(final_scaled, (size_t) n_final * (size_t) p, 0.0f);
  std::vector<int> y_final(n_final);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* src = x_scaled.data() + (size_t) j * (size_t) n;
    float* dst = final_scaled.data() + (size_t) j * (size_t) n_final;
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
