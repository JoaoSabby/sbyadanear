/*
 * sby_hpc_engine.cpp
 *
 * Motor HPC consolidado do pacote sbyadanear.
 *
 * Decisoes de desempenho:
 *   - estatisticas populacionais por kernel SIMD paralelo (Fortran).
 *   - matriz de distancias por D^2 = ||A||^2 + ||B||^2 - 2 A B^T com sgemm
 *     (interface Fortran do BLAS, portanto valida tanto no oneMKL quanto na
 *     BLAS de referencia do R).
 *   - blocagem bidimensional das chamadas sgemm com piso de tile, para nao
 *     degenerar em milhares de produtos finos.
 *   - os blocos de distancia saem orientados como (referencia x consulta), de
 *     modo que o laco interno de manutencao do top-k le memoria contigua.
 *   - interpolacao lambda do ADASYN com uniformes gerados por Rcpp::runif no espaco padronizado.
 *   - despadronizacao das sinteticas por FMA AVX-512 inteiramente no C++.
 *   - NearMiss processa a minoria em blocos (GEMM particionado) sem materializar
 *     a matriz ampliada (minoria + sinteticas).
 *   - reconstrucao final do tibble ocorre na camada R.
 *
 * Afinidade NUMA: os buffers quentes sao tocados pela mesma particao
 * schedule(static) que depois os consome, o que distribui as paginas entre os
 * dois sockets por first touch. Para que isso valha, o processo precisa rodar
 * com as threads fixadas, por exemplo com OMP_PROC_BIND=close, OMP_PLACES=cores
 * ou sob "numactl --interleave=all". O pacote nao impoe essa politica: fixar
 * threads e responsabilidade do ambiente de execucao.
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
#include <climits>
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

// Preenche o buffer inteiro com um valor. Use somente quando o valor inicial e
// de fato lido depois (por exemplo o infinito dos acumuladores de top-k).
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

// Reserva o buffer sem inicializar os valores, tocando apenas uma posicao por
// pagina de 4 KB. Isso preserva a distribuicao NUMA por first touch (que e o
// motivo real do laco paralelo) sem pagar a escrita completa de buffers que
// serao integralmente sobrescritos logo em seguida, como a saida do sgemm com
// beta = 0. O alocador alinhado faz default-init para tipos triviais, portanto
// resize() nao zera nada por conta propria.
template <typename Buffer>
static void sby_resize_first_touch_pages(Buffer& buffer, size_t n){
  buffer.resize(n);
  if(n == 0) return;
  using value_type = typename Buffer::value_type;
  value_type* data = buffer.data();
  const size_t stride = 4096u / sizeof(value_type) > 0 ? 4096u / sizeof(value_type) : 1u;
  const size_t pages = (n + stride - 1) / stride;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(size_t page = 0; page < pages; ++page){
    data[page * stride] = value_type();
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


// -------------------------------------------------------------------
// Blocagem das chamadas sgemm.
//
// A versao anterior dimensionava apenas a dimensao de referencia a partir do
// numero total de linhas de consulta. Em matrizes altas isso colapsava o bloco
// para poucas colunas e produzia milhares de sgemm finos, nos quais o custo de
// empacotamento do BLAS domina o produto. Agora as duas dimensoes sao
// blocadas, com um piso de tile que mantem cada chamada larga o suficiente
// para o MKL amortizar o empacotamento.
// -------------------------------------------------------------------
struct sby_gemm_tiling {
  int query_block;
  int ref_block;
};

static sby_gemm_tiling sby_resolve_gemm_tiling(int n_query, int n_ref, int p){
  const size_t dist_budget_elems = (size_t) 32 * 1024 * 1024;  // 128 MB em float
  const size_t copy_budget_elems = (size_t) 16 * 1024 * 1024;  //  64 MB em float
  const size_t tile_floor        = 1024;

  const size_t dim_query = (size_t) std::max(1, n_query);
  const size_t dim_ref   = (size_t) std::max(1, n_ref);
  const size_t width     = (size_t) std::max(1, p);

  // A copia dos blocos limita cada dimensao isoladamente.
  size_t ref_block = copy_budget_elems / width;
  if(ref_block < tile_floor) ref_block = tile_floor;
  if(ref_block > dim_ref)    ref_block = dim_ref;

  size_t query_block = copy_budget_elems / width;
  if(query_block < tile_floor) query_block = tile_floor;
  if(query_block > dim_query)  query_block = dim_query;

  // O bloco de distancias limita o produto das duas dimensoes. Encolhe primeiro
  // a consulta e so depois a referencia.
  if(query_block * ref_block > dist_budget_elems){
    size_t allowed = dist_budget_elems / ref_block;
    if(allowed < tile_floor) allowed = tile_floor;
    if(allowed < query_block) query_block = allowed;
  }
  if(query_block * ref_block > dist_budget_elems){
    size_t allowed = dist_budget_elems / query_block;
    if(allowed < 1) allowed = 1;
    if(allowed < ref_block) ref_block = allowed;
  }

  sby_gemm_tiling out;
  out.query_block = (int) std::max<size_t>(1, std::min(query_block, dim_query));
  out.ref_block   = (int) std::max<size_t>(1, std::min(ref_block, dim_ref));
  return out;
}

// Copia um bloco contiguo de linhas de uma matriz column major. O laco de copia
// e o proprio first touch das paginas, portanto nao ha pre-zeragem: todos os
// n_block * p elementos sao escritos aqui.
static void sby_copy_column_block(const sby_float_buffer& source, int n_source,
                                  int p, int row_start, int n_block,
                                  sby_float_buffer& block){
  block.resize((size_t) n_block * (size_t) p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* src = source.data() + (size_t) j * (size_t) n_source + (size_t) row_start;
    float* dst = block.data() + (size_t) j * (size_t) n_block;
    sby_parallel_copy(src, dst, (size_t) n_block);
  }
}

// Copia linhas arbitrarias (gather) de uma matriz column major para um bloco
// contiguo. Tambem faz o first touch pelo proprio laco de escrita.
static void sby_gather_row_block(const sby_float_buffer& source, int n_source,
                                 int p, const int* row_index, int n_block,
                                 sby_float_buffer& block){
  block.resize((size_t) n_block * (size_t) p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j = 0; j < p; ++j){
    const float* col = source.data() + (size_t) j * (size_t) n_source;
    float* dst = block.data() + (size_t) j * (size_t) n_block;
    for(int r = 0; r < n_block; ++r){
      dst[r] = col[row_index[r]];
    }
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

// omp_set_num_threads muda o estado global do runtime OpenMP e vale para todo o
// processo, nao apenas para a chamada corrente. Sem restauracao, um unico
// sby_max_threads baixo contaminava permanentemente qualquer outro codigo
// paralelo da sessao R. O guard salva o valor anterior e o repoe na saida,
// inclusive quando a chamada aborta por excecao.
class sby_omp_thread_guard {
public:
  explicit sby_omp_thread_guard(SEXP max_threads){
#ifdef _OPENMP
    int requested = sby_resolve_native_threads(max_threads);
    if(requested > 0){
      previous_ = omp_get_max_threads();
      omp_set_num_threads(requested);
      active_ = true;
    }
#else
    (void) max_threads;
#endif
  }

  ~sby_omp_thread_guard(){
#ifdef _OPENMP
    if(active_){
      omp_set_num_threads(previous_);
    }
#endif
  }

  sby_omp_thread_guard(const sby_omp_thread_guard&) = delete;
  sby_omp_thread_guard& operator=(const sby_omp_thread_guard&) = delete;

private:
#ifdef _OPENMP
  int  previous_ = 1;
  bool active_   = false;
#endif
};

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
  // O kernel Fortran escreve todos os n * p elementos; o first touch fica por
  // conta do toque por pagina, sem a escrita completa redundante.
  sby_resize_first_touch_pages(x_scaled, (size_t) n * (size_t) p);
  int status = 0;
  sby_apply_zscore_simd_f(x, n, p, means.data(), sds.data(), x_scaled.data(), &status);
  if(status != 0){
    Rcpp::stop("Falha na aplicacao de z-score (status=%d)", status);
  }
}

// -------------------------------------------------------------------
// sby_destandarize_synthetic
// Reverte o z-score das linhas sinteticas (float32) escrevendo direto no
// NumericMatrix de saida. Ambos sao column major com o mesmo n_syn x p, entao
// nao ha motivo para alocar um buffer intermediario e copia-lo em seguida.
// -------------------------------------------------------------------
static void sby_destandarize_synthetic(
    const sby_float_buffer& syn_scaled,   // column major n_syn x p
    int n_syn, int p,
    const sby_double_buffer& means,
    const sby_double_buffer& sds,
    Rcpp::NumericMatrix& out){

  int status = 0;
  sby_revert_zscore_fma_f(syn_scaled.data(), n_syn, p,
                          means.data(), sds.data(), out.begin(), &status);
  if(status != 0){
    Rcpp::stop("Falha na despadronizacao das sinteticas (status=%d)", status);
  }
}

// -------------------------------------------------------------------
// sby_knn_topk_against_reference
//
// Retorna, para cada linha de consulta, os indices 1-based (na referencia) dos
// out_k vizinhos mais proximos exatos, em ordem crescente de distancia.
//
// self_ref_index, quando nao nulo, da a linha 0-based da referencia que e a
// propria linha de consulta e deve ser descartada. Isso generaliza o antigo
// "drop_self" por igualdade de indice, que so valia quando consulta e
// referencia eram literalmente a mesma matriz.
//
// A consulta e a referencia sao blocadas em duas dimensoes e a matriz de
// distancias sai orientada como (referencia x consulta): assim o laco interno
// de manutencao do top-k percorre memoria contigua e o laco externo sobre as
// consultas e paralelizavel sem escrita compartilhada.
// -------------------------------------------------------------------
static void sby_knn_topk_against_reference(
    const sby_float_buffer& query, int n_query,
    const sby_float_buffer& reference, int n_ref, int p,
    int k, const int* self_ref_index,
    std::vector<int>& out_index,
    int& out_k){
  const bool drop_self = (self_ref_index != nullptr);
  const int  max_k     = n_ref - (drop_self ? 1 : 0);

  out_k = std::min(k, max_k);
  if(out_k < 1 || n_query < 1 || n_ref < 1){
    out_k = 0;
    out_index.clear();
    return;
  }
  // Sentinela negativa: nenhum indice valido pode permanecer sem escrita. A
  // versao anterior inicializava com 1 e, quando um bloco nao preenchia o
  // top-k, devolvia silenciosamente a primeira linha como vizinho.
  out_index.assign((size_t) n_query * (size_t) out_k, -1);

  int keep_k = out_k + (drop_self ? 1 : 0);
  if(keep_k > n_ref) keep_k = n_ref;

  const sby_gemm_tiling tiling = sby_resolve_gemm_tiling(n_query, n_ref, p);
  const double infinity = std::numeric_limits<double>::infinity();

  sby_float_buffer query_block;
  sby_float_buffer reference_block;
  sby_float_buffer dist_block;
  sby_double_buffer top_dist;
  sby_double_buffer worst_val;
  std::vector<int> top_index;
  std::vector<int> worst_pos;
  int status = 0;
  int incomplete = 0;

  for(int q_start = 0; q_start < n_query; q_start += tiling.query_block){
    const int n_q = std::min(tiling.query_block, n_query - q_start);
    sby_copy_column_block(query, n_query, p, q_start, n_q, query_block);

    sby_resize_first_touch(top_dist, (size_t) n_q * (size_t) keep_k, infinity);
    top_index.assign((size_t) n_q * (size_t) keep_k, -1);
    worst_pos.assign(n_q, 0);
    sby_resize_first_touch(worst_val, (size_t) n_q, infinity);

    for(int ref_start = 0; ref_start < n_ref; ref_start += tiling.ref_block){
      const int n_r = std::min(tiling.ref_block, n_ref - ref_start);
      sby_copy_column_block(reference, n_ref, p, ref_start, n_r, reference_block);
      // dist_block sai como (n_r x n_q): coluna i guarda as distancias da
      // consulta i para todas as referencias do bloco, em posicoes contiguas.
      sby_resize_first_touch_pages(dist_block, (size_t) n_r * (size_t) n_q);
      sby_pairwise_sqdist_sgemm_f(reference_block.data(), n_r,
                                  query_block.data(), n_q, p,
                                  dist_block.data(), &status);
      if(status != 0){
        Rcpp::stop("Falha no calculo blocado de distancias por sgemm (status=%d)", status);
      }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for(int i = 0; i < n_q; ++i){
        const float* col = dist_block.data() + (size_t) i * (size_t) n_r;
        const size_t row_offset = (size_t) i * (size_t) keep_k;
        double* row_dist = top_dist.data() + row_offset;
        int*    row_idx  = top_index.data() + row_offset;
        int    local_worst_pos = worst_pos[i];
        double local_worst_val = worst_val[i];
        for(int b = 0; b < n_r; ++b){
          const double candidate = (double) col[b];
          if(candidate >= local_worst_val) continue;
          row_dist[local_worst_pos] = candidate;
          row_idx[local_worst_pos]  = ref_start + b;

          int new_worst_pos = 0;
          double new_worst_val = row_dist[0];
          for(int t = 1; t < keep_k; ++t){
            if(row_dist[t] > new_worst_val){
              new_worst_val = row_dist[t];
              new_worst_pos = t;
            }
          }
          local_worst_pos = new_worst_pos;
          local_worst_val = new_worst_val;
        }
        worst_pos[i] = local_worst_pos;
        worst_val[i] = local_worst_val;
      }
    }

#ifdef _OPENMP
#pragma omp parallel reduction(|:incomplete)
#endif
    {
      // Buffer de ordenacao por thread: evita uma alocacao por linha de consulta.
      std::vector<int> order(keep_k);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for(int i = 0; i < n_q; ++i){
        const size_t row_offset = (size_t) i * (size_t) keep_k;
        const double* row_dist = top_dist.data() + row_offset;
        const int*    row_idx  = top_index.data() + row_offset;
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b){
          if(row_dist[a] != row_dist[b]) return row_dist[a] < row_dist[b];
          return row_idx[a] < row_idx[b];
        });

        const int self_row = drop_self ? self_ref_index[q_start + i] : -1;
        int written = 0;
        for(int pos : order){
          const int cand = row_idx[pos];
          if(cand < 0) continue;
          if(cand == self_row) continue;
          out_index[(size_t) (q_start + i) * (size_t) out_k + (size_t) written] = cand + 1;
          ++written;
          if(written >= out_k) break;
        }
        if(written < out_k) incomplete = 1;
      }
    }
  }

  if(incomplete){
    Rcpp::stop("Busca exata de vizinhos nao completou o top-k solicitado");
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
  sby_gather_row_block(x_scaled_orig, n, p, minority_index.data(), n_min, minority);

  // O ratio ADASYN e sempre somado a 1 para determinar o tamanho final da
  // classe rara. Somente a diferenca em relacao aos originais e sintetizada.
  double final_minority_double = std::floor((double) n_min * (1.0 + over_ratio));
  if(final_minority_double > (double) INT_MAX){
    Rcpp::stop("quantidade final ADASYN excede o limite suportado");
  }
  // Nao ha piso artificial de uma sintetica: com ratio 0 o contrato e nao
  // gerar nada, exatamente como na rota classica em R.
  int final_minority_count = (int) final_minority_double;
  if(final_minority_count < n_min) final_minority_count = n_min;
  int synthetic_count = final_minority_count - n_min;

  syn_codes_out.assign(synthetic_count, minority_code);
  if(synthetic_count < 1){
    syn_scaled_out.clear();
    return 0;
  }

  // -----------------------------------------------------------------
  // Ponderacao ADASYN por densidade local (He et al., 2008).
  //
  // Para cada linha rara conta-se a fracao de vizinhos que pertencem a classe
  // majoritaria (r_i). Quanto mais cercada de maioria, mais dificil o exemplo
  // e mais sinteticas ele recebe. A implementacao espelha a rota classica em
  // R/sby_generate_adasyn_samples.R, incluindo o fallback uniforme quando
  // nenhuma linha rara tem vizinho majoritario.
  // -----------------------------------------------------------------
  int global_k = k_neighbor;
  if(global_k > n - 1) global_k = n - 1;
  if(global_k < 1) global_k = 1;

  std::vector<int> global_neighbors;
  int global_neighbor_k = 0;
  sby_knn_topk_against_reference(minority, n_min, x_scaled_orig, n, p,
                                 global_k, minority_index.data(),
                                 global_neighbors, global_neighbor_k);

  std::vector<double> majority_ratio(n_min, 0.0);
  double ratio_sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:ratio_sum)
#endif
  for(int i = 0; i < n_min; ++i){
    int majority_hits = 0;
    for(int t = 0; t < global_neighbor_k; ++t){
      const int neighbor_row = global_neighbors[(size_t) i * (size_t) global_neighbor_k + (size_t) t] - 1;
      if(y_codes_orig[neighbor_row] != minority_code){
        ++majority_hits;
      }
    }
    const double ratio = (double) majority_hits / (double) global_neighbor_k;
    majority_ratio[i] = ratio;
    ratio_sum += ratio;
  }

  std::vector<double> raw_counts(n_min);
  std::vector<int> synthetic_per_row(n_min);
  long long assigned = 0;
  for(int i = 0; i < n_min; ++i){
    const double weight = (ratio_sum > 0.0)
      ? majority_ratio[i] / ratio_sum
      : 1.0 / (double) n_min;
    raw_counts[i] = (double) synthetic_count * weight;
    synthetic_per_row[i] = (int) std::floor(raw_counts[i]);
    assigned += synthetic_per_row[i];
  }

  // Distribui o resto pelos maiores residuos fracionarios. O desempate por
  // indice crescente reproduz a estabilidade de order(..., decreasing = TRUE)
  // do R, mantendo o resultado deterministico.
  long long remaining = (long long) synthetic_count - assigned;
  if(remaining > 0){
    std::vector<int> fractional_order(n_min);
    std::iota(fractional_order.begin(), fractional_order.end(), 0);
    std::sort(fractional_order.begin(), fractional_order.end(), [&](int a, int b){
      const double fa = raw_counts[a] - (double) synthetic_per_row[a];
      const double fb = raw_counts[b] - (double) synthetic_per_row[b];
      if(fa != fb) return fa > fb;
      return a < b;
    });
    for(long long t = 0; t < remaining; ++t){
      synthetic_per_row[fractional_order[(size_t) (t % (long long) n_min)]] += 1;
    }
  }

  int effective_k = k_neighbor;
  if(effective_k > n_min - 1) effective_k = n_min - 1;
  if(effective_k < 1) effective_k = 1;

  std::vector<int> minority_self(n_min);
  std::iota(minority_self.begin(), minority_self.end(), 0);
  std::vector<int> minority_neighbors;
  int minority_neighbor_k = 0;
  sby_knn_topk_against_reference(minority, n_min, minority, n_min, p,
                                 effective_k, minority_self.data(),
                                 minority_neighbors, minority_neighbor_k);

  std::vector<int> base_idx(synthetic_count);
  std::vector<int> nbr_idx(synthetic_count);
  sby_float_buffer lambda;
  lambda.resize((size_t) synthetic_count);

  Rcpp::NumericVector unif_lambda = Rcpp::runif(synthetic_count);
  Rcpp::NumericVector unif_pick   = Rcpp::runif(synthetic_count);
  int s = 0;
  for(int i = 0; i < n_min && s < synthetic_count; ++i){
    for(int rep = 0; rep < synthetic_per_row[i] && s < synthetic_count; ++rep, ++s){
      base_idx[s] = i + 1;
      int pick = (int) std::floor(unif_pick[s] * (double) minority_neighbor_k);
      if(pick < 0) pick = 0;
      if(pick >= minority_neighbor_k) pick = minority_neighbor_k - 1;
      nbr_idx[s] = minority_neighbors[(size_t) i * (size_t) minority_neighbor_k + (size_t) pick];
      lambda[s]  = (float) unif_lambda[s];
    }
  }
  if(s != synthetic_count){
    Rcpp::stop("Alocacao ponderada ADASYN nao cobriu todas as linhas sinteticas");
  }

  sby_resize_first_touch_pages(syn_scaled_out, (size_t) synthetic_count * (size_t) p);
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
  // Sem piso artificial: se o ratio zera a maioria, o contrato e o mesmo da
  // rota classica em R, que rejeita a configuracao em vez de reter uma linha
  // arbitraria para disfarcar o problema.
  if(retained_majority < 1 && n_maj > 0){
    Rcpp::stop("'sby_nearmiss_ratio' reteve zero linhas da classe majoritaria");
  }

  sby_float_buffer minority;
  minority.resize((size_t) n_min * (size_t) p);
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
  // A maioria e a dimensao de consulta e a minoria a de referencia. As duas sao
  // blocadas pelo mesmo orcamento, com piso de tile.
  const sby_gemm_tiling tiling = sby_resolve_gemm_tiling(n_maj, n_min, p);
  sby_float_buffer majority_block;
  sby_float_buffer minority_block;
  sby_float_buffer dist_block;
  int status = 0;

  for(int maj_start = 0; maj_start < n_maj; maj_start += tiling.query_block){
    int n_maj_block = std::min(tiling.query_block, n_maj - maj_start);
    sby_gather_row_block(x_scaled_orig, n, p, majority_index.data() + maj_start,
                         n_maj_block, majority_block);

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

    for(int min_start = 0; min_start < n_min; min_start += tiling.ref_block){
      int n_block = std::min(tiling.ref_block, n_min - min_start);
      sby_copy_column_block(minority, n_min, p, min_start, n_block, minority_block);
      // Orientacao (minoria x maioria): para uma linha majoritaria fixa, as
      // distancias contra o bloco minoritario ficam contiguas. A orientacao
      // anterior fazia o laco interno saltar n_maj_block floats por leitura,
      // desperdicando toda a linha de cache carregada.
      sby_resize_first_touch_pages(dist_block, (size_t) n_block * (size_t) n_maj_block);
      sby_pairwise_sqdist_sgemm_f(minority_block.data(), n_block,
                                  majority_block.data(), n_maj_block, p,
                                  dist_block.data(), &status);
      if(status != 0){
        Rcpp::stop("Falha no calculo blocado de distancias NearMiss por sgemm (status=%d)", status);
      }

      if(effective_k == n_min){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int m = 0; m < n_maj_block; ++m){
          const float* col = dist_block.data() + (size_t) m * (size_t) n_block;
          double acc = 0.0;
#ifdef _OPENMP
#pragma omp simd reduction(+:acc)
#endif
          for(int c = 0; c < n_block; ++c){
            acc += (double) col[c];
          }
          score_acc[m] += acc;
        }
      } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int m = 0; m < n_maj_block; ++m){
          const float* col = dist_block.data() + (size_t) m * (size_t) n_block;
          double* row_topk = topk.data() + (size_t) m * (size_t) effective_k;
          int local_worst_pos = worst_pos[m];
          double local_worst_val = worst_val[m];
          for(int c = 0; c < n_block; ++c){
            const double candidate = (double) col[c];
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
BEGIN_RCPP
  Rcpp::List out;
  // Derivado das macros que o compilador realmente define, e nao de um -D fixo
  // no Makevars: antes o relatorio anunciava suporte Cascade Lake mesmo quando
  // nenhuma flag de arquitetura chegava ao compilador.
#if defined(__AVX512F__) && defined(__AVX512CD__) && defined(__AVX512BW__) && \
    defined(__AVX512DQ__) && defined(__AVX512VL__)
  out["cascade_lake_native"] = true;
#else
  out["cascade_lake_native"] = false;
#endif
#if defined(SBYADANEAR_ONEAPI_MKL)
  out["mkl_linked"] = true;
#else
  out["mkl_linked"] = false;
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
END_RCPP
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
BEGIN_RCPP
  // Rcpp::runif so respeita set.seed() com o estado do RNG do R sincronizado.
  // Em um ponto de entrada .Call() puro isso nao acontece sozinho, e sem o
  // escopo abaixo a saida deixava de ser reprodutivel por sby_seed.
  Rcpp::RNGScope rng_scope;
  const sby_omp_thread_guard thread_guard(max_threads);

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
END_RCPP
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
BEGIN_RCPP
  Rcpp::RNGScope rng_scope;
  const sby_omp_thread_guard thread_guard(max_threads);

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
END_RCPP
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
BEGIN_RCPP
  const sby_omp_thread_guard thread_guard(max_threads);

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
END_RCPP
}
