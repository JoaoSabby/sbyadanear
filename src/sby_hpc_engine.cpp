/*
 * sby_hpc_engine.cpp
 *
 * Motor HPC consolidado do pacote sbyadanear.
 *
 * Decisoes de desempenho:
 *   - estatisticas populacionais via Vector Statistics Library (Fortran).
 *   - matriz de distancias por D^2 = ||A||^2 + ||B||^2 - 2 A B^T com cblas_sgemm.
 *   - interpolacao lambda do ADASYN com uniformes gerados por Rcpp::runif no espaco padronizado.
 *   - despadronizacao das sinteticas por FMA AVX-512 inteiramente no C++.
 *   - NearMiss processa a minoria em blocos (GEMM particionado) sem materializar
 *     a matriz ampliada (minoria + sinteticas).
 *   - reconstrucao final do tibble ocorre na camada R.
 *
 * Contrato de retorno:
 *   sby_adanear_hpc_result_cpp -> List(
 *     sby_synthetic_rows        = NumericMatrix  (double, despadronizado),
 *     sby_retained_majority_idx = IntegerVector  (1-based),
 *     sby_target_synthetic      = IntegerVector  (codigos de nivel),
 *     sby_scaling_info          = List(centers, scales)
 *   )
 *   sby_adasyn_hpc_cpp -> List(
 *     sby_synthetic_rows   = NumericMatrix  (double, despadronizado),
 *     sby_target_synthetic = IntegerVector  (codigos de nivel),
 *     sby_scaling_info     = List(centers, scales)
 *   )
 *   sby_nearmiss_hpc_cpp -> List(
 *     sby_retained_majority_idx = IntegerVector (1-based),
 *     sby_scaling_info          = List(centers, scales)
 *   )
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
using sby_float_buffer  = std::vector<float,  sby_aligned_allocator<float,  64> >;

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
  if(by_dist < 1) by_dist = 1;
  size_t by_copy = target_bytes / (sizeof(float) * (size_t) std::max(1, p));
  if(by_copy < 1) by_copy = 1;
  size_t block = std::min(by_dist, by_copy);
  if(block > (size_t) n_ref) block = (size_t) n_ref;
  if(block < 1) block = 1;
  return (int) block;
}

static int sby_resolve_nearmiss_majority_block(int n_maj, int n_min, int p){
  const size_t target_dist_bytes = (size_t) 128 * (size_t) 1024 * (size_t) 1024;
  const size_t target_copy_bytes = (size_t) 64 * (size_t) 1024 * (size_t) 1024;
  int min_block = sby_resolve_gemm_ref_block(std::max(1, std::min(n_maj, 65536)), n_min, p);
  size_t by_dist = target_dist_bytes /
    (sizeof(float) * (size_t) std::max(1, min_block));
  size_t by_copy = target_copy_bytes /
    (sizeof(float) * (size_t) std::max(1, p));
  size_t block = std::min(by_dist, by_copy);
  if(block < 1024) block = 1024;
  if(block > (size_t) n_maj) block = (size_t) n_maj;
  if(block < 1) block = 1;
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
  for(int j = 0; j < p; ++j){
    if(!(sds[j] > 0.0) || !std::isfinite(sds[j])){
      sds[j] = 1.0;
    }
  }
}

// -------------------------------------------------------------------
// sby_apply_zscore
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
// sby_destandarize_synthetic
// Reverte o z-score das linhas sinteticas (float32) e escreve o
// resultado como double no NumericMatrix de saida. Executado inteiramente
// no C++; a camada R recebe os dados ja nas unidades originais.
// -------------------------------------------------------------------
static void sby_destandarize_synthetic(
    const sby_float_buffer& syn_scaled,   // column major n_syn x p
    int n_syn, int p,
    const sby_double_buffer& means,
    const sby_double_buffer& sds,
    Rcpp::NumericMatrix& out){

  sby_double_buffer tmp;
  sby_resize_first_touch(tmp, (size_t) n_syn * (size_t) p, 0.0);
  int status = 0;
  sby_revert_zscore_fma_f(syn_scaled.data(), n_syn, p,
                          means.data(), sds.data(), tmp.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha na despadronizacao das sinteticas (status=%d)", status);
  }
  // Copia column major -> NumericMatrix (column major identico)
  double* dst = out.begin();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const double* src_col = tmp.data() + (size_t) j * (size_t) n_syn;
    double*       dst_col = dst         + (size_t) j * (size_t) n_syn;
    sby_parallel_copy(src_col, dst_col, (size_t) n_syn);
  }
}

// -------------------------------------------------------------------
// sby_knn_topk_against_reference
// -------------------------------------------------------------------
static void sby_knn_topk_against_reference(
    const sby_float_buffer& query, int n_query,
    const sby_float_buffer& reference, int n_ref, int p,
    int k, bool drop_self,
    std::vector<int>& out_index,
    int& out_k){
  out_k = k;
  if(out_k > n_ref) out_k = n_ref;
  if(out_k < 1 || n_query < 1 || n_ref < 1){
    out_k = 0;
    out_index.clear();
    return;
  }
  out_index.assign((size_t) n_query * (size_t) out_k, 1);

  int keep_k = out_k + (drop_self ? 1 : 0);
  if(keep_k > n_ref) keep_k = n_ref;

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
        if(candidate >= worst_val[i]) continue;
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

  std::vector<int> order(keep_k);
  for(int i = 0; i < n_query; ++i){
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

    int written = 0;
    for(int pos : order){
      int cand = top_index[row_offset + (size_t) pos];
      if(cand < 0) continue;
      if(drop_self && cand == i) continue;
      out_index[(size_t) i * (size_t) out_k + (size_t) written] = cand + 1;
      ++written;
      if(written >= out_k) break;
    }
  }
}

// -------------------------------------------------------------------
// sby_run_adasyn_stage
// Retorna o buffer de sinteticas em float32 (column major) e os codigos
// de nivel correspondentes. NAO anexa as sinteticas ao buffer original
// (o NearMiss opera apenas sobre os originais).
// -------------------------------------------------------------------
static int sby_run_adasyn_stage(
    const sby_float_buffer& x_scaled_orig, int n, int p,
    const std::vector<int>& y_codes_orig, int minority_code,
    int k_neighbor, double over_ratio,
    sby_float_buffer& syn_scaled_out,
    std::vector<int>& syn_codes_out){

  std::vector<int> minority_index;
  minority_index.reserve(n);
  for(int i = 0; i < n; ++i){
    if(y_codes_orig[i] == minority_code){
      minority_index.push_back(i);
    }
  }
  int n_min = (int) minority_index.size();
  if(n_min < 2){
    Rcpp::stop("ADASYN exige ao menos 2 observacoes minoritarias");
  }

  sby_float_buffer minority;
  sby_resize_first_touch(minority, (size_t) n_min * (size_t) p, 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* col = x_scaled_orig.data() + (size_t) j * (size_t) n;
    float* dst = minority.data() + (size_t) j * (size_t) n_min;
    for(int r = 0; r < n_min; ++r){
      dst[r] = col[minority_index[r]];
    }
  }

  // Numero de sinteticas: floor(n_min * over_ratio), minimo 1.
  // Mantem a semantica historica de expandir a minoria sem reduzir
  // observacoes raras originais.
  int synthetic_count = (int) std::floor((double) n_min * over_ratio);
  if(synthetic_count < 1) synthetic_count = 1;

  int effective_k = k_neighbor;
  if(effective_k > n_min - 1) effective_k = n_min - 1;
  if(effective_k < 1) effective_k = 1;

  std::vector<int> minority_neighbors;
  int minority_neighbor_k = 0;
  sby_knn_topk_against_reference(minority, n_min, minority, n_min, p,
                                 effective_k, true, minority_neighbors,
                                 minority_neighbor_k);

  std::vector<int> base_idx(synthetic_count);
  std::vector<int> nbr_idx(synthetic_count);
  sby_float_buffer lambda;
  lambda.resize((size_t) synthetic_count);

  Rcpp::NumericVector unif_lambda = Rcpp::runif(synthetic_count);
  Rcpp::NumericVector unif_pick   = Rcpp::runif(synthetic_count);
  for(int s = 0; s < synthetic_count; ++s){
    int base = s % n_min;
    base_idx[s] = base + 1;
    int pick = (int) std::floor(unif_pick[s] * (double) minority_neighbor_k);
    if(pick < 0) pick = 0;
    if(pick >= minority_neighbor_k) pick = minority_neighbor_k - 1;
    nbr_idx[s] = minority_neighbors[(size_t) base * (size_t) minority_neighbor_k + (size_t) pick];
    lambda[s]  = (float) unif_lambda[s];
  }

  sby_resize_first_touch(syn_scaled_out, (size_t) synthetic_count * (size_t) p, 0.0f);
  int status = 0;
  sby_adasyn_interp_uniform_f(minority.data(), n_min, p,
                              base_idx.data(), nbr_idx.data(), lambda.data(),
                              synthetic_count, syn_scaled_out.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha na interpolacao sintetica ADASYN (status=%d)", status);
  }

  syn_codes_out.assign(synthetic_count, minority_code);
  return synthetic_count;
}

// -------------------------------------------------------------------
// sby_run_nearmiss_stage
// Recebe apenas os originais (sem sinteticas). Retorna os indices 0-based
// das linhas majoritarias retidas (dentro de x_scaled_orig).
// GEMM particionado sobre a minoria original: sem materializar concatenacao.
// -------------------------------------------------------------------
static std::vector<int> sby_run_nearmiss_stage(
    const sby_float_buffer& x_scaled_orig, int n, int p,
    const std::vector<int>& y_codes_orig, int minority_code,
    int k_neighbor, double under_ratio,
    const sby_float_buffer* extra_minority = nullptr, int n_extra_minority = 0){

  std::vector<int> minority_index;
  std::vector<int> majority_index;
  minority_index.reserve(n);
  majority_index.reserve(n);
  for(int i = 0; i < n; ++i){
    if(y_codes_orig[i] == minority_code){
      minority_index.push_back(i);
    } else {
      majority_index.push_back(i);
    }
  }
  int n_min_orig = (int) minority_index.size();
  int n_maj = (int) majority_index.size();
  if(extra_minority == nullptr || n_extra_minority < 1){
    n_extra_minority = 0;
  }
  int n_min = n_min_orig + n_extra_minority;

  // sby_under_ratio opera sobre a quantidade minoritaria apos oversampling
  // quando sinteticas sao fornecidas. Valores menores que 1 podem reter uma
  // maioria menor que a minoria; as linhas raras nunca sao descartadas aqui.
  int retained_majority = (int) std::floor((double) n_min * under_ratio);
  if(retained_majority > n_maj) retained_majority = n_maj;
  if(retained_majority < 1 && n_maj > 0) retained_majority = 1;

  sby_float_buffer minority;
  sby_resize_first_touch(minority, (size_t) n_min * (size_t) p, 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* col = x_scaled_orig.data() + (size_t) j * (size_t) n;
    float* ndst = minority.data() + (size_t) j * (size_t) n_min;
    for(int r = 0; r < n_min_orig; ++r){ ndst[r] = col[minority_index[r]]; }
    if(n_extra_minority > 0){
      const float* extra_col = extra_minority->data() + (size_t) j * (size_t) n_extra_minority;
      for(int r = 0; r < n_extra_minority; ++r){
        ndst[n_min_orig + r] = extra_col[r];
      }
    }
  }

  int effective_k = k_neighbor;
  if(effective_k > n_min) effective_k = n_min;
  if(effective_k < 1)     effective_k = 1;

  std::vector< std::pair<double,int> > scores(n_maj);
  int majority_block_size = sby_resolve_nearmiss_majority_block(n_maj, n_min, p);
  sby_float_buffer majority_block;
  sby_float_buffer minority_block;
  sby_float_buffer dist_block;
  int status = 0;

  for(int maj_start = 0; maj_start < n_maj; maj_start += majority_block_size){
    int n_maj_block = std::min(majority_block_size, n_maj - maj_start);
    sby_resize_first_touch(majority_block, (size_t) n_maj_block * (size_t) p, 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int j = 0; j < p; ++j){
      const float* col = x_scaled_orig.data() + (size_t) j * (size_t) n;
      float* mdst = majority_block.data() + (size_t) j * (size_t) n_maj_block;
      for(int r = 0; r < n_maj_block; ++r){
        mdst[r] = col[majority_index[maj_start + r]];
      }
    }

    int minority_block_size = sby_resolve_gemm_ref_block(n_maj_block, n_min, p);
    sby_double_buffer score_acc;
    sby_resize_first_touch(score_acc, (size_t) n_maj_block, 0.0);

    sby_double_buffer topk;
    sby_double_buffer worst_val;
    std::vector<int> worst_pos;
    if(effective_k < n_min){
      sby_resize_first_touch(topk, (size_t) n_maj_block * (size_t) effective_k,
                             std::numeric_limits<double>::infinity());
      worst_pos.assign(n_maj_block, 0);
      sby_resize_first_touch(worst_val, (size_t) n_maj_block,
                             std::numeric_limits<double>::infinity());
    }

    for(int min_start = 0; min_start < n_min; min_start += minority_block_size){
      int n_block = std::min(minority_block_size, n_min - min_start);
      sby_copy_column_block(minority, n_min, p, min_start, n_block, minority_block);
      sby_resize_first_touch(dist_block, (size_t) n_maj_block * (size_t) n_block, 0.0f);
      sby_pairwise_sqdist_sgemm_f(majority_block.data(), n_maj_block,
                                  minority_block.data(), n_block, p,
                                  dist_block.data(), &status);
      if(status != 0){
        Rcpp::stop("Falha no calculo blocado de distancias NearMiss por sgemm (status=%d)", status);
      }

      if(effective_k == n_min){
        for(int c = 0; c < n_block; ++c){
          const float* col = dist_block.data() + (size_t) c * (size_t) n_maj_block;
#ifdef _OPENMP
#pragma omp simd
#endif
          for(int m = 0; m < n_maj_block; ++m){
            score_acc[m] += col[m];
          }
        }
      } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int m = 0; m < n_maj_block; ++m){
          double* row_topk = topk.data() + (size_t) m * (size_t) effective_k;
          int local_worst_pos = worst_pos[m];
          double local_worst_val = worst_val[m];
          for(int c = 0; c < n_block; ++c){
            const float candidate = dist_block[(size_t) c * (size_t) n_maj_block + (size_t) m];
            if(candidate >= local_worst_val) continue;
            row_topk[local_worst_pos] = candidate;

            int new_worst_pos = 0;
            double new_worst_val = row_topk[0];
            for(int t = 1; t < effective_k; ++t){
              if(row_topk[t] > new_worst_val){
                new_worst_val = row_topk[t];
                new_worst_pos = t;
              }
            }
            local_worst_pos = new_worst_pos;
            local_worst_val = new_worst_val;
          }
          worst_pos[m] = local_worst_pos;
          worst_val[m] = local_worst_val;
        }
      }
    }

    if(effective_k < n_min){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for(int m = 0; m < n_maj_block; ++m){
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

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int m = 0; m < n_maj_block; ++m){
      scores[maj_start + m] = std::make_pair(
        score_acc[m] / (double) effective_k,
        majority_index[maj_start + m]
      );
    }
  }

  auto score_less = [](const std::pair<double,int>& a, const std::pair<double,int>& b){
    if(a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  };
  if(retained_majority > 0 && retained_majority < n_maj){
    auto middle = scores.begin() + retained_majority;
    std::nth_element(scores.begin(), middle, scores.end(), score_less);
    std::sort(scores.begin(), middle, score_less);
  } else {
    std::sort(scores.begin(), scores.end(), score_less);
  }

  std::vector<int> retained_maj_0based;
  retained_maj_0based.reserve(retained_majority);
  for(int t = 0; t < retained_majority; ++t){
    retained_maj_0based.push_back(scores[t].second);
  }
  return retained_maj_0based;
}

static Rcpp::IntegerVector sby_extract_factor_codes(SEXP y){
  return Rcpp::IntegerVector(y);
}

static Rcpp::List sby_build_scaling_info(const sby_double_buffer& means,
                                         const sby_double_buffer& sds,
                                         int p){
  Rcpp::NumericVector centers(p), scales(p);
  for(int j = 0; j < p; ++j){
    centers[j] = means[j];
    scales[j]  = sds[j];
  }
  return Rcpp::List::create(
    Rcpp::Named("centers") = centers,
    Rcpp::Named("scales")  = scales
  );
}


//' @title Relatorio de compilacao do motor HPC
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
  out["openmp"]          = true;
  out["openmp_version"]  = _OPENMP;
  out["openmp_max_threads"] = omp_get_max_threads();
#else
  out["openmp"]             = false;
  out["openmp_version"]     = NA_INTEGER;
  out["openmp_max_threads"] = NA_INTEGER;
#endif
  return out;
}


//' @title Pipeline ADASYN + NearMiss-1 HPC
//' @description Retorna lista com sinteticas despadronizadas, indices da maioria
//'   retida e metadados de escala. Tibble montado na camada R.
// [[Rcpp::export]]
extern "C" SEXP sby_adanear_hpc_result_cpp(
    SEXP x_matrix, SEXP y_factor,
    SEXP k_adanear, SEXP k_nearmiss,
    SEXP over_ratio, SEXP under_ratio,
    SEXP max_threads, SEXP column_names,
    SEXP target_levels){

  sby_apply_native_threads(max_threads);

  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow(), p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels   = levels.size();
  int k_over     = Rcpp::as<int>(k_adanear);
  int k_under    = Rcpp::as<int>(k_nearmiss);
  double r_over  = Rcpp::as<double>(over_ratio);
  double r_under = Rcpp::as<double>(under_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  sby_double_buffer means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  sby_float_buffer x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());

  // ADASYN sobre os originais — sinteticas ficam em buffer separado
  sby_float_buffer syn_scaled;
  std::vector<int> syn_codes;
  int n_syn = sby_run_adasyn_stage(x_scaled, n, p, y_codes, minority_code,
                                   k_over, r_over, syn_scaled, syn_codes);

  // NearMiss usa a minoria apos oversampling como referencia de distancia e
  // como base para a contagem de maioria retida, sem descartar raras.
  std::vector<int> retained_maj_0based = sby_run_nearmiss_stage(
    x_scaled, n, p, y_codes, minority_code, k_under, r_under,
    &syn_scaled, n_syn
  );
  int n_ret_maj = (int) retained_maj_0based.size();

  // Despadronizacao das sinteticas no C++
  Rcpp::NumericMatrix syn_out(n_syn, p);
  syn_out.attr("dimnames") = Rcpp::List::create(R_NilValue, column_names);
  if(n_syn > 0){
    sby_destandarize_synthetic(syn_scaled, n_syn, p, means, sds, syn_out);
  }

  // Codigos de nivel das sinteticas
  Rcpp::IntegerVector target_synthetic(n_syn);
  for(int s = 0; s < n_syn; ++s){
    target_synthetic[s] = syn_codes[s];
  }

  // Indices 1-based da maioria retida
  Rcpp::IntegerVector retained_majority_idx(n_ret_maj);
  for(int r = 0; r < n_ret_maj; ++r){
    retained_majority_idx[r] = retained_maj_0based[r] + 1;
  }

  return Rcpp::List::create(
    Rcpp::Named("sby_synthetic_rows")        = syn_out,
    Rcpp::Named("sby_retained_majority_idx") = retained_majority_idx,
    Rcpp::Named("sby_target_synthetic")      = target_synthetic,
    Rcpp::Named("sby_scaling_info")          = sby_build_scaling_info(means, sds, p)
  );
}


//' @title ADASYN HPC
//' @description Retorna sinteticas despadronizadas e metadados de escala.
//   Tibble montado na camada R.
// [[Rcpp::export]]
extern "C" SEXP sby_adasyn_hpc_cpp(
    SEXP x_matrix, SEXP y_factor,
    SEXP k_adanear, SEXP over_ratio,
    SEXP max_threads, SEXP column_names,
    SEXP target_levels){

  sby_apply_native_threads(max_threads);

  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow(), p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels = levels.size();
  int k_over   = Rcpp::as<int>(k_adanear);
  double r_over = Rcpp::as<double>(over_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  sby_double_buffer means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  sby_float_buffer x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());

  sby_float_buffer syn_scaled;
  std::vector<int> syn_codes;
  int n_syn = sby_run_adasyn_stage(x_scaled, n, p, y_codes, minority_code,
                                   k_over, r_over, syn_scaled, syn_codes);

  Rcpp::NumericMatrix syn_out(n_syn, p);
  syn_out.attr("dimnames") = Rcpp::List::create(R_NilValue, column_names);
  if(n_syn > 0){
    sby_destandarize_synthetic(syn_scaled, n_syn, p, means, sds, syn_out);
  }

  Rcpp::IntegerVector target_synthetic(n_syn);
  for(int s = 0; s < n_syn; ++s){
    target_synthetic[s] = syn_codes[s];
  }

  return Rcpp::List::create(
    Rcpp::Named("sby_synthetic_rows")   = syn_out,
    Rcpp::Named("sby_target_synthetic") = target_synthetic,
    Rcpp::Named("sby_scaling_info")     = sby_build_scaling_info(means, sds, p)
  );
}


//' @title NearMiss-1 HPC
//' @description Retorna indices 1-based da maioria retida e metadados de escala.
//   Tibble montado na camada R.
// [[Rcpp::export]]
extern "C" SEXP sby_nearmiss_hpc_cpp(
    SEXP x_matrix, SEXP y_factor,
    SEXP k_nearmiss, SEXP under_ratio,
    SEXP max_threads, SEXP column_names,
    SEXP target_levels){

  sby_apply_native_threads(max_threads);

  Rcpp::NumericMatrix x(x_matrix);
  int n = x.nrow(), p = x.ncol();
  Rcpp::IntegerVector y_codes_in = sby_extract_factor_codes(y_factor);
  Rcpp::CharacterVector levels(target_levels);
  int n_levels  = levels.size();
  int k_under   = Rcpp::as<int>(k_nearmiss);
  double r_under = Rcpp::as<double>(under_ratio);

  int minority_code = sby_resolve_minority_role(y_codes_in, n_levels);

  sby_double_buffer means, sds;
  sby_zscore_population(x.begin(), n, p, means, sds);
  sby_float_buffer x_scaled;
  sby_apply_zscore(x.begin(), n, p, means, sds, x_scaled);

  std::vector<int> y_codes(y_codes_in.begin(), y_codes_in.end());

  std::vector<int> retained_maj_0based = sby_run_nearmiss_stage(
    x_scaled, n, p, y_codes, minority_code, k_under, r_under
  );
  int n_ret_maj = (int) retained_maj_0based.size();

  Rcpp::IntegerVector retained_majority_idx(n_ret_maj);
  for(int r = 0; r < n_ret_maj; ++r){
    retained_majority_idx[r] = retained_maj_0based[r] + 1;
  }

  return Rcpp::List::create(
    Rcpp::Named("sby_retained_majority_idx") = retained_majority_idx,
    Rcpp::Named("sby_scaling_info")          = sby_build_scaling_info(means, sds, p)
  );
}
