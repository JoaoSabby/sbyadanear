/*
 * Projeto: adanear
 * Autor: Joao Batista Goncalves de Brito
 *
 * Justificativa Arquitetural para Implementacao em C:
 * Esta biblioteca nativa foi desenvolvida para transpor gargalos criticos de
 * performance inerentes a execucao de loops massivos no ambiente R puro.
 * Algoritmos de superamostragem como o ADASYN e rotinas de padronizacao
 * exigem manipulacao intensiva de tensores multidimensionais.
 *
 * A integracao em baixo nivel fornece controle absoluto sobre a alocacao da heap,
 * permitindo operacoes diretas em memoria via ponteiros escalares sem acionar
 * copias intermediarias de objetos. Alem de garantir tempos de execucao drasticamente
 * menores em pipelines de modelagem preditiva, a base em C blinda a estabilidade
 * matematica ao utilizar variaveis de precisao estendida no tratamento de
 * variancias proximas a zero, isolando as operacoes do coletor de lixo da maquina virtual R.
 *
 * Referencias Historicas Fundamentais:
 * 1. He, H., Bai, Y., Garcia, E. A., & Li, S. (2008). ADASYN: Adaptive synthetic
 *    sampling approach for imbalanced learning. IEEE International Joint Conference
 *    on Neural Networks. (Base fundacional do algoritmo de superamostragem implementado).
 *
 * 2. Knuth, D. E. (1997). The Art of Computer Programming, Volume 2: Seminumerical
 *    Algorithms. Addison Wesley. (Base fundacional historica para a estabilidade
 *    numerica e algoritmos seguros de calculo de variancia em ponto flutuante).
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <R_ext/Visibility.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <vector>

#if defined(SBYADANEAR_ONEAPI_MKL)
#include <mkl.h>
#endif

#if defined(SBYADANEAR_ONEDPL)
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#endif

/* FCONE: macro de R-API para passar comprimentos de strings Fortran nas
 * chamadas a BLAS. Definida em Rconfig.h em R >= 3.6; fornece um fallback
 * vazio em versoes mais antigas para manter portabilidade. */
#ifndef FCONE
# define FCONE
#endif

/* Distancias ao quadrado calculadas por ||q||^2 + ||r||^2 - 2 q.r sao
 * extremamente rapidas via BLAS, mas sofrem cancelamento catastrofico quando
 * os vetores sao quase identicos ou muito colineares em escala alta. Nesses
 * pares raros, recalculamos a distancia diretamente sobre as diferencas para
 * preservar a micro-ordenacao dos vizinhos em ADASYN/NearMiss. */
#define SBY_D2_RESCUE_EPS 1e-11

struct nearmiss_score {
  double mean;
  int index;
};

struct knn_neighbor {
  double dist2;
  int index;
};

static bool nearmiss_score_less(const nearmiss_score& lhs, const nearmiss_score& rhs){
  return lhs.mean < rhs.mean || (lhs.mean == rhs.mean && lhs.index < rhs.index);
}

static bool knn_neighbor_less(const knn_neighbor& lhs, const knn_neighbor& rhs){
  return lhs.dist2 < rhs.dist2 || (lhs.dist2 == rhs.dist2 && lhs.index < rhs.index);
}

static bool knn_neighbor_worse(const knn_neighbor& lhs, const knn_neighbor& rhs){
  return lhs.dist2 < rhs.dist2 || (lhs.dist2 == rhs.dist2 && lhs.index < rhs.index);
}

/* Selecao dos k menores: evita chamar std::nth_element caso k == n_ref (middle == last) */
template <typename iterator, typename compare>
static void sby_select_topk(iterator first, iterator middle, iterator last, compare comp){
  if (middle != last) {
    std::nth_element(first, middle, last, comp);
  }
}

/* Ordenacao final dos k primeiros: oneDPL quando habilitado e vantajoso. */
template <typename iterator, typename compare>
static void sby_sort_topk(iterator first, iterator last, compare comp,
                          std::size_t k, std::size_t n_ref){
#if defined(SBYADANEAR_ONEDPL)
  /* Limiar heuristico: ativa par_unseq somente quando a granularidade (k)
   * justifica o overhead de despacho do backend TBB/OpenMP. */
  if(k >= 256){
    oneapi::dpl::sort(oneapi::dpl::execution::par_unseq, first, last, comp);
  } else {
    std::sort(first, last, comp);
  }
#else
  std::sort(first, last, comp);
#endif
}

/* sby_partial_sort: mantido apenas para vetores de nearmiss_score. */
template <typename iterator, typename compare>
static void sby_partial_sort(iterator first, iterator middle, iterator last, compare comp){
#if defined(SBYADANEAR_ONEDPL)
  oneapi::dpl::partial_sort(oneapi::dpl::execution::par_unseq, first, middle, last, comp);
#else
  std::partial_sort(first, middle, last, comp);
#endif
}

static void sby_compute_column_norms2(const double *matrix, double *norms, R_xlen_t rows, R_xlen_t cols){
  std::fill(norms, norms + rows, 0.0);
  double *squared = (double *) R_alloc((size_t) rows, sizeof(double));
  for(R_xlen_t c = 0; c < cols; ++c){
    const double *col = matrix + c * rows;
#if defined(SBYADANEAR_ONEAPI_MKL)
    if(rows <= INT_MAX){
      vdSqr((MKL_INT) rows, col, squared);
      vdAdd((MKL_INT) rows, norms, squared, norms);
    } else
#endif
    {
#pragma omp simd
      for(R_xlen_t i = 0; i < rows; ++i){
        norms[i] += col[i] * col[i];
      }
    }
  }
}

static double euclidean_d2_with_rescue(
  double d2,
  const double *q_ref,
  const double *r_ref,
  R_xlen_t q_idx,
  R_xlen_t r_idx,
  R_xlen_t n_query,
  R_xlen_t n_ref,
  R_xlen_t col_count
){
  if(d2 < SBY_D2_RESCUE_EPS){
    long double exact_dist = 0.0L;
#pragma omp simd reduction(+:exact_dist)
    for(R_xlen_t c = 0; c < col_count; ++c){
      const long double diff = (long double) q_ref[q_idx + c * n_query] -
                               (long double) r_ref[r_idx + c * n_ref];
      exact_dist += diff * diff;
    }
    d2 = (double) exact_dist;
  }
  if(d2 < 0.0){
    d2 = 0.0;
  }
  return d2;
}

//' @title Verificação de interrupção do usuário
//' @description Encaminha a checagem cooperativa de interrupção para a API do R.
//' @return Valor nulo do R.
extern "C" SEXP check_user_interrupt_c(void){
  R_CheckUserInterrupt();
  return R_NilValue;
}

static void require_real_matrix(SEXP x, const char *name){
  if (!isReal(x) || !isMatrix(x)) {
    error("O parametro '%s' deve ser uma matrix double", name);
  }
}

static void require_finite_real_matrix_values(const double *values, R_xlen_t length, const char *name){
  for(R_xlen_t i = 0; i < length; ++i){
    if(!R_FINITE(values[i])){
      error("O parametro '%s' nao pode conter NA, NaN, Inf ou -Inf", name);
    }
  }
}

static void require_real_vector_length(SEXP x, R_xlen_t expected_length, const char *name){
  if(!isReal(x) || XLENGTH(x) != expected_length){
    error("O parametro '%s' deve ser um vetor double com comprimento compativel", name);
  }
}

//' @title Parâmetros nativos de Z-Score
//' @description Calcula médias e desvios-padrão amostrais por coluna com aceleração vetorial quando disponível.
//' @param x_matrix Matriz double de entrada.
//' @return Lista com centros e escalas.
extern "C" SEXP compute_z_score_params_c(SEXP x_matrix){
  require_real_matrix(x_matrix, "x_matrix");

  SEXP dims = getAttrib(x_matrix, R_DimSymbol);
  const R_xlen_t n = (R_xlen_t) INTEGER(dims)[0];
  const R_xlen_t p = (R_xlen_t) INTEGER(dims)[1];

  if(n < 2){
    error("'x_matrix' deve conter ao menos duas linhas");
  }

  SEXP centers = PROTECT(allocVector(REALSXP, p));
  SEXP scales  = PROTECT(allocVector(REALSXP, p));

  const double *x = REAL(x_matrix);
  double *mu = REAL(centers);
  double *sd = REAL(scales);

  for(R_xlen_t j = 0; j < p; ++j){
    if((j & 15) == 0){
      R_CheckUserInterrupt();
    }

    const double *col = x + j * n;
    long double sum = 0.0L;

#pragma omp simd reduction(+:sum)
    for(R_xlen_t i = 0; i < n; ++i){
      sum += (long double) col[i];
    }

    const long double mean = sum / (long double) n;
    mu[j] = (double) mean;

    long double ss = 0.0L;
    double *center_vec  = (double *) R_alloc((size_t) n, sizeof(double));
    double *centered_vec = (double *) R_alloc((size_t) n, sizeof(double));
    std::fill(center_vec, center_vec + n, (double) mean);
#if defined(SBYADANEAR_ONEAPI_MKL)
    if(n <= INT_MAX){
      vdSub((MKL_INT) n, col, center_vec, centered_vec);
      vdSqr((MKL_INT) n, centered_vec, centered_vec);
    } else
#endif
    {
#pragma omp simd
      for(R_xlen_t i = 0; i < n; ++i){
        const double centered = col[i] - (double) mean;
        centered_vec[i] = centered * centered;
      }
    }

#pragma omp simd reduction(+:ss)
    for(R_xlen_t i = 0; i < n; ++i){
      ss += (long double) centered_vec[i];
    }

    sd[j] = sqrt((double) (ss / (long double) (n - 1)));
  }

  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("centers"));
  SET_STRING_ELT(names, 1, mkChar("scales"));

  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, centers);
  SET_VECTOR_ELT(out, 1, scales);
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(4);
  return out;
}

//' @title Transformação Z-Score vetorizada
//' @description Aplica ou reverte a padronização Z-Score utilizando primitivas vetoriais oneMath quando habilitadas.
//' @param x_matrix Matriz double de entrada.
//' @param centers Vetor double de médias.
//' @param scales Vetor double de desvios-padrão.
//' @param reverse Indicador lógico de reversão.
//' @return Matriz double transformada.
extern "C" SEXP apply_z_score_c(SEXP x_matrix, SEXP centers, SEXP scales, SEXP reverse){
  require_real_matrix(x_matrix, "x_matrix");

  if(!isReal(centers) || !isReal(scales)){
    error("'centers' e 'scales' devem ser vetores double");
  }

  SEXP dims = getAttrib(x_matrix, R_DimSymbol);
  const R_xlen_t n = (R_xlen_t) INTEGER(dims)[0];
  const R_xlen_t p = (R_xlen_t) INTEGER(dims)[1];

  if(XLENGTH(centers) != p || XLENGTH(scales) != p){
    error("'centers' e 'scales' devem ter comprimento igual ao numero de colunas");
  }

  if(n > INT_MAX || p > INT_MAX){
    error("'x_matrix' excede o limite de dimensoes suportado por allocMatrix");
  }
  SEXP out = PROTECT(allocMatrix(REALSXP, (int) n, (int) p));

  const double *x  = REAL(x_matrix);
  const double *mu = REAL(centers);
  const double *sd = REAL(scales);
  double *y = REAL(out);

  const int do_reverse = asLogical(reverse);

  for(R_xlen_t j = 0; j < p; ++j){
    if((j & 15) == 0){
      R_CheckUserInterrupt();
    }

    const R_xlen_t offset = j * n;
    const double center = mu[j];
    const double scale  = sd[j];

    if(!do_reverse && (!R_FINITE(scale) || scale <= 0.0)){
      error("'scales[%lld]' deve ser positivo e finito", (long long) (j + 1));
    }
    if(!R_FINITE(center)){
      error("'centers[%lld]' deve finito", (long long) (j + 1));
    }

    if(do_reverse){
#if defined(SBYADANEAR_ONEAPI_MKL)
      if(n <= INT_MAX){
        double *scale_vec  = (double *) R_alloc((size_t) n, sizeof(double));
        double *center_vec = (double *) R_alloc((size_t) n, sizeof(double));
        std::fill(scale_vec,  scale_vec  + n, scale);
        std::fill(center_vec, center_vec + n, center);
        vdMul((MKL_INT) n, x + offset, scale_vec,  y + offset);
        vdAdd((MKL_INT) n, y + offset, center_vec, y + offset);
      } else
#endif
      {
#pragma omp simd
        for(R_xlen_t i = 0; i < n; ++i){
          y[offset + i] = x[offset + i] * scale + center;
        }
      }
    } else {
#if defined(SBYADANEAR_ONEAPI_MKL)
      if(n <= INT_MAX){
        double *center_vec = (double *) R_alloc((size_t) n, sizeof(double));
        double *scale_vec  = (double *) R_alloc((size_t) n, sizeof(double));
        std::fill(center_vec, center_vec + n, center);
        std::fill(scale_vec,  scale_vec  + n, scale);
        vdSub((MKL_INT) n, x + offset, center_vec, y + offset);
        vdDiv((MKL_INT) n, y + offset, scale_vec,  y + offset);
      } else
#endif
      {
#pragma omp simd
        for(R_xlen_t i = 0; i < n; ++i){
          y[offset + i] = (x[offset + i] - center) / scale;
        }
      }
    }
  }

  setAttrib(out, R_DimNamesSymbol, getAttrib(x_matrix, R_DimNamesSymbol));
  UNPROTECT(1);
  return out;
}

//' @title Síntese ADASYN nativa por linha
//' @description Gera observações sintéticas a partir da interpolação entre instâncias minoritárias e seus vizinhos.
//' @return Matriz double sintética.
extern "C" SEXP generate_synthetic_adasyn_c(SEXP minority_matrix, SEXP minority_neighbor_index, SEXP synthetic_per_row){
  require_real_matrix(minority_matrix, "minority_matrix");

  if(!isInteger(minority_neighbor_index) || !isMatrix(minority_neighbor_index)){
    error("'minority_neighbor_index' deve ser uma matrix integer");
  }
  if(!isInteger(synthetic_per_row)){
    error("'synthetic_per_row' deve ser um vetor integer");
  }

  SEXP minority_dims = getAttrib(minority_matrix, R_DimSymbol);
  SEXP neighbor_dims = getAttrib(minority_neighbor_index, R_DimSymbol);

  const int minority_rows  = INTEGER(minority_dims)[0];
  const int col_count      = INTEGER(minority_dims)[1];
  const int neighbor_rows  = INTEGER(neighbor_dims)[0];
  const int neighbor_count = INTEGER(neighbor_dims)[1];

  if(neighbor_rows != minority_rows || XLENGTH(synthetic_per_row) != minority_rows){
    error("Dimensoes inconsistentes para geracao ADASYN");
  }
  if(neighbor_count < 1){
    error("'minority_neighbor_index' deve conter ao menos uma coluna");
  }

  const double *minority = REAL(minority_matrix);
  const int    *neighbor = INTEGER(minority_neighbor_index);
  const int    *per_row  = INTEGER(synthetic_per_row);

  R_xlen_t total_synthetic = 0;
  for(int i = 0; i < minority_rows; ++i){
    if ((i & 8191) == 0){
      R_CheckUserInterrupt();
    }
    if(per_row[i] < 0){
      error("'synthetic_per_row' nao pode conter valores negativos");
    }
    total_synthetic += (R_xlen_t) per_row[i];
  }

  if(total_synthetic > INT_MAX){
    error("Numero de linhas sinteticas excede o limite suportado por matrix R");
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, (int) total_synthetic, col_count));
  double *synthetic = REAL(out);

  GetRNGstate();

  int write_row = 0;
  for(int i = 0; i < minority_rows; ++i){
    if((i & 255) == 0 && i > 0){
      PutRNGstate();
      R_CheckUserInterrupt();
      GetRNGstate();
    }

    const int row_count = per_row[i];
    for(int r = 0; r < row_count; ++r){
      int sampled_neighbor_column = (int) floor(unif_rand() * (double) neighbor_count);
      if(sampled_neighbor_column >= neighbor_count){
        sampled_neighbor_column = neighbor_count - 1;
      }

      const int selected_neighbor_row = neighbor[i + ((R_xlen_t) sampled_neighbor_column * minority_rows)] - 1;
      if(selected_neighbor_row < 0 || selected_neighbor_row >= minority_rows){
        PutRNGstate();
        error("'minority_neighbor_index' contem indice fora do intervalo");
      }

      const double weight = unif_rand();
      for(int j = 0; j < col_count; ++j){
        const double base_value     = minority[i + ((R_xlen_t) j * minority_rows)];
        const double neighbor_value = minority[selected_neighbor_row + ((R_xlen_t) j * minority_rows)];
        synthetic[write_row + ((R_xlen_t) j * total_synthetic)] = base_value + weight * (neighbor_value - base_value);
      }
      ++write_row;
    }
  }

  PutRNGstate();
  UNPROTECT(1);
  return out;
}

//' @title Seleção NearMiss acelerada
//' @description Seleciona exemplos majoritários por médias de distância usando ordenação parcial oneDPL quando habilitada.
//' @return Vetor inteiro de índices retidos.
extern "C" SEXP select_nearmiss_majority_c(SEXP nn_dist, SEXP majority_index, SEXP retained_majority_count){
  require_real_matrix(nn_dist, "nn_dist");

  if(!isInteger(majority_index)){
    error("'majority_index' deve ser integer");
  }
  if(!isInteger(retained_majority_count) || LENGTH(retained_majority_count) != 1){
    error("'retained_majority_count' deve ser integer escalar");
  }

  SEXP dims    = getAttrib(nn_dist, R_DimSymbol);
  const int n  = INTEGER(dims)[0];
  const int k  = INTEGER(dims)[1];
  const int retain = INTEGER(retained_majority_count)[0];

  if(LENGTH(majority_index) != n){
    error("'majority_index' deve ter comprimento igual a nrow(nn_dist)");
  }
  if(retain < 1 || retain > n){
    error("'retained_majority_count' deve estar entre 1 e nrow(nn_dist)");
  }
  if(k < 1){
    error("'nn_dist' deve conter ao menos uma coluna");
  }

  const double *dist   = REAL(nn_dist);
  const int    *majority = INTEGER(majority_index);
  std::vector<nearmiss_score> scores(static_cast<std::size_t>(n));

  for(int i = 0; i < n; ++i){
    if((i & 8191) == 0){
      R_CheckUserInterrupt();
    }
    long double sum = 0.0L;
#pragma omp simd reduction(+:sum)
    for(int j = 0; j < k; ++j){
      const double value = dist[i + ((R_xlen_t) j * n)];
      if(!R_FINITE(value)){
        error("'nn_dist' contem distancia ausente ou infinita");
      }
      sum += (long double) value;
    }
    scores[static_cast<std::size_t>(i)] = {
      (double) (sum / (long double) k),
      majority[i]
    };
  }

  /* sby_partial_sort mantido aqui: vetor de nearmiss_score, n_query elementos */
  sby_partial_sort(scores.begin(), scores.begin() + retain, scores.end(), nearmiss_score_less);

  SEXP out = PROTECT(allocVector(INTSXP, retain));
  for(int i = 0; i < retain; ++i){
    INTEGER(out)[i] = scores[static_cast<std::size_t>(i)].index;
  }

  UNPROTECT(1);
  return out;
}

//' @title Síntese ADASYN nativa por coluna
//' @description Gera observações sintéticas em layout favorável à localidade de colunas.
//' @return Matriz double sintética.
extern "C" SEXP generate_synthetic_adasyn_col_c(SEXP minority_matrix, SEXP minority_neighbor_index, SEXP synthetic_per_row){
  require_real_matrix(minority_matrix, "minority_matrix");

  if(!isInteger(minority_neighbor_index) || !isMatrix(minority_neighbor_index)){
    error("'minority_neighbor_index' deve ser uma matrix integer");
  }
  if(!isInteger(synthetic_per_row)){
    error("'synthetic_per_row' deve ser um vetor integer");
  }

  SEXP minority_dims = getAttrib(minority_matrix, R_DimSymbol);
  SEXP neighbor_dims = getAttrib(minority_neighbor_index, R_DimSymbol);
  const R_xlen_t minority_rows  = (R_xlen_t) INTEGER(minority_dims)[0];
  const R_xlen_t col_count      = (R_xlen_t) INTEGER(minority_dims)[1];
  const R_xlen_t neighbor_rows  = (R_xlen_t) INTEGER(neighbor_dims)[0];
  const R_xlen_t neighbor_count = (R_xlen_t) INTEGER(neighbor_dims)[1];

  if(neighbor_rows != minority_rows || XLENGTH(synthetic_per_row) != minority_rows){
    error("Dimensoes inconsistentes para geracao ADASYN");
  }
  if(neighbor_count < 1){
    error("'minority_neighbor_index' deve conter ao menos uma coluna");
  }

  const double *minority = REAL(minority_matrix);
  const int    *neighbor = INTEGER(minority_neighbor_index);
  const int    *per_row  = INTEGER(synthetic_per_row);

  R_xlen_t total_synthetic = 0;
  for(R_xlen_t i = 0; i < minority_rows; ++i){
    if(per_row[i] < 0){
      error("'synthetic_per_row' nao pode conter valores negativos");
    }
    total_synthetic += (R_xlen_t) per_row[i];
  }
  if(total_synthetic > INT_MAX){
    error("Numero de linhas sinteticas excede o limite suportado por matrix R");
  }
  if(col_count > INT_MAX){
    error("Numero de colunas excede o limite suportado por matrix R");
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, (int) total_synthetic, (int) col_count));
  double *synthetic = REAL(out);

  if(total_synthetic == 0){
    UNPROTECT(1);
    return out;
  }

  int    *base_row = (int    *) R_alloc((size_t) total_synthetic, sizeof(int));
  int    *nbr_row  = (int    *) R_alloc((size_t) total_synthetic, sizeof(int));
  double *wgt      = (double *) R_alloc((size_t) total_synthetic, sizeof(double));

  GetRNGstate();
  R_xlen_t s = 0;
  for(R_xlen_t i = 0; i < minority_rows; ++i){
    if((i & 8191) == 0){
      R_CheckUserInterrupt();
    }
    const int row_count = per_row[i];
    for(int r = 0; r < row_count; ++r){
      int sampled_col = (int) floor(unif_rand() * (double) neighbor_count);
      if(sampled_col >= (int) neighbor_count){
        sampled_col = (int) neighbor_count - 1;
      }
      const int selected_neighbor_row = neighbor[i + (R_xlen_t) sampled_col * minority_rows] - 1;
      if(selected_neighbor_row < 0 || selected_neighbor_row >= (int) minority_rows){
        PutRNGstate();
        error("'minority_neighbor_index' contem indice fora do intervalo");
      }
      base_row[s] = (int) i;
      nbr_row[s]  = selected_neighbor_row;
      wgt[s]      = unif_rand();
      ++s;
    }
  }
  PutRNGstate();

  for(R_xlen_t j = 0; j < col_count; ++j){
    if((j & 15) == 0){
      R_CheckUserInterrupt();
    }
    const R_xlen_t col_offset_min = j * minority_rows;
    const R_xlen_t col_offset_syn = j * total_synthetic;
    for(R_xlen_t t = 0; t < total_synthetic; ++t){
      const double bv = minority[base_row[t] + col_offset_min];
      const double nv = minority[nbr_row[t]  + col_offset_min];
      synthetic[t + col_offset_syn] = bv + wgt[t] * (nv - bv);
    }
  }

  UNPROTECT(1);
  return out;
}

//' @title Remoção do próprio vizinho
//' @description Remove o índice da própria observação de uma matriz KNN registrada.
//' @return Matriz inteira de vizinhos filtrados.
extern "C" SEXP drop_self_neighbor_c(SEXP nbr, SEXP self_index, SEXP desired_k){
  if(!isInteger(nbr) || !isMatrix(nbr)){
    error("'nbr' deve ser uma matrix integer");
  }
  if(!isInteger(self_index)){
    error("'self_index' deve ser integer");
  }
  if(!isInteger(desired_k) || LENGTH(desired_k) != 1){
    error("'desired_k' deve ser integer escalar");
  }

  SEXP dims             = getAttrib(nbr, R_DimSymbol);
  const R_xlen_t n      = (R_xlen_t) INTEGER(dims)[0];
  const R_xlen_t k_plus = (R_xlen_t) INTEGER(dims)[1];
  const int k           = INTEGER(desired_k)[0];

  if(XLENGTH(self_index) != n){
    error("'self_index' deve ter comprimento igual a nrow(nbr)");
  }
  if(k < 1){
    error("'desired_k' deve ser >= 1");
  }
  if((R_xlen_t) k > k_plus){
    error("Nao foi possivel remover o proprio ponto mantendo vizinhos suficientes");
  }
  if(n > INT_MAX){
    error("'nbr' excede o limite suportado por allocMatrix");
  }

  SEXP out     = PROTECT(allocMatrix(INTSXP, (int) n, k));
  int *outp    = INTEGER(out);
  const int *src = INTEGER(nbr);

  for(R_xlen_t i = 0; i < n; ++i){
    if((i & 8191) == 0){
      R_CheckUserInterrupt();
    }
    const int self = INTEGER(self_index)[i];
    R_xlen_t written = 0;
    R_xlen_t j = 0;
    while(written < k && j < k_plus){
      const int v = src[i + j * n];
      if(v != NA_INTEGER && v != self){
        outp[i + written * n] = v;
        ++written;
      }
      ++j;
    }
    if(written < k){
      UNPROTECT(1);
      error("Nao foi possivel remover o proprio ponto mantendo vizinhos suficientes");
    }
  }

  UNPROTECT(1);
  return out;
}

/**
 * @brief Top-k brute force KNN para distancia euclidiana usando BLAS.
 *
 * Calcula D^2[i,j] = ||q_i||^2 + ||r_j||^2 - 2 q_i . r_j em blocos de
 * consulta, via dgemm para o termo cruzado, e usa selecao top-k via
 * nth_element + sort para extrair os k vizinhos com menor distancia.
 */
static SEXP brute_force_knn_impl(SEXP data, SEXP query, SEXP k_value,
                                  int return_index, int return_dist,
                                  int query_is_data = 0,
                                  int exclude_self  = 0,
                                  R_xlen_t query_offset = 0){
  require_real_matrix(data,  "data");
  require_real_matrix(query, "query");
  if(!isInteger(k_value) || LENGTH(k_value) != 1){
    error("'k' deve ser integer escalar");
  }
  if(!return_index && !return_dist){
    error("Ao menos um componente de KNN deve ser solicitado");
  }

  SEXP d_dims = getAttrib(data,  R_DimSymbol);
  SEXP q_dims = getAttrib(query, R_DimSymbol);
  const R_xlen_t n_ref   = (R_xlen_t) INTEGER(d_dims)[0];
  const R_xlen_t p1      = (R_xlen_t) INTEGER(d_dims)[1];
  const R_xlen_t n_query = (R_xlen_t) INTEGER(q_dims)[0];
  const R_xlen_t p2      = (R_xlen_t) INTEGER(q_dims)[1];
  const int k            = INTEGER(k_value)[0];

  if(p1 != p2){
    error("'data' e 'query' devem ter o mesmo numero de colunas");
  }
  if(k < 1){
    error("'k' deve ser >= 1");
  }
  if(n_ref > INT_MAX){
    error("'data' excede o limite de indices inteiros suportado");
  }
  if((R_xlen_t) k > n_ref){
    error("'k' nao pode exceder o numero de linhas de 'data'");
  }
  if(exclude_self == 1){
    if(query_is_data != 1){
      error("'exclude_self' requer 'query_is_data = TRUE'");
    }
    if((R_xlen_t) k > n_ref - 1){
      error("'k' nao pode exceder nrow(data) - 1 quando 'exclude_self = TRUE'");
    }
  }
  if(n_query > INT_MAX){
    error("'query' excede o limite suportado por allocMatrix");
  }

  const double *r_ref = REAL(data);
  const double *q_ref = REAL(query);
  require_finite_real_matrix_values(r_ref, n_ref   * p1, "data");
  require_finite_real_matrix_values(q_ref, n_query * p2, "query");

  double *q_norm2 = (double *) R_alloc((size_t) n_query, sizeof(double));
  double *r_norm2 = (double *) R_alloc((size_t) n_ref,   sizeof(double));
  sby_compute_column_norms2(q_ref, q_norm2, n_query, p1);
  sby_compute_column_norms2(r_ref, r_norm2, n_ref,   p1);

  int protect_count = 0;
  SEXP out_idx = R_NilValue;
  SEXP out_dst = R_NilValue;
  int    *idx_out = NULL;
  double *dst_out = NULL;

  if(return_index){
    out_idx = PROTECT(allocMatrix(INTSXP,  (int) n_query, k)); ++protect_count;
    idx_out = INTEGER(out_idx);
  }
  if(return_dist){
    out_dst = PROTECT(allocMatrix(REALSXP, (int) n_query, k)); ++protect_count;
    dst_out = REAL(out_dst);
  }

  const R_xlen_t max_cells = 8 * 1024 * 1024;
  R_xlen_t block_q = max_cells / (n_ref > 0 ? n_ref : 1);
  if(block_q < 1)       block_q = 1;
  if(block_q > n_query) block_q = n_query;
  if(block_q > INT_MAX) block_q = INT_MAX;

  double *cross = (double *) R_alloc((size_t) (block_q * n_ref), sizeof(double));

  for(R_xlen_t q_start = 0; q_start < n_query; q_start += block_q){
    R_xlen_t q_end = q_start + block_q;
    if(q_end > n_query) q_end = n_query;
    const R_xlen_t cur_b = q_end - q_start;
    R_CheckUserInterrupt();

    const double alpha    = 1.0;
    const double beta     = 0.0;
    const int m_int   = (int) cur_b;
    const int n_int   = (int) n_ref;
    const int kk_int  = (int) p1;
    const int lda_int = (int) n_query;
    const int ldb_int = (int) n_ref;
    const int ldc_int = (int) cur_b;
    F77_CALL(dgemm)(
      "N", "T",
      &m_int, &n_int, &kk_int,
      &alpha,
      q_ref + q_start, &lda_int,
      r_ref,           &ldb_int,
      &beta,
      cross,           &ldc_int FCONE FCONE
    );

    std::vector<knn_neighbor> neighbors(static_cast<std::size_t>(n_ref));
    for(R_xlen_t bi = 0; bi < cur_b; ++bi){
      const R_xlen_t q_idx = q_start + bi;
      const double   qn2   = q_norm2[q_idx];

      for(R_xlen_t j = 0; j < n_ref; ++j){
        if((exclude_self == 1) && (query_offset + q_idx) == j){
          neighbors[static_cast<std::size_t>(j)] = {R_PosInf, (int) (j + 1)};
          continue;
        }
        double d2 = qn2 + r_norm2[j] - 2.0 * cross[bi + j * cur_b];
        d2 = euclidean_d2_with_rescue(d2, q_ref, r_ref, q_idx, j, n_query, n_ref, p1);
        neighbors[static_cast<std::size_t>(j)] = {d2, (int) (j + 1)};
      }

      auto first  = neighbors.begin();
      auto middle = neighbors.begin() + k;
      auto last   = neighbors.end();

      sby_select_topk(first, middle, last, knn_neighbor_less);
      sby_sort_topk(first, middle, knn_neighbor_less,
                    static_cast<std::size_t>(k),
                    static_cast<std::size_t>(n_ref));

      for(int kk = 0; kk < k; ++kk){
        const knn_neighbor& neighbor = neighbors[static_cast<std::size_t>(kk)];
        if(return_index){
          idx_out[q_idx + (R_xlen_t) kk * n_query] = neighbor.index;
        }
        if(return_dist){
          dst_out[q_idx + (R_xlen_t) kk * n_query] = sqrt(neighbor.dist2);
        }
      }
    }
  }

  const int out_length = return_index && return_dist ? 2 : 1;
  SEXP names = PROTECT(allocVector(STRSXP, out_length)); ++protect_count;
  SEXP out   = PROTECT(allocVector(VECSXP, out_length)); ++protect_count;

  int pos = 0;
  if(return_index){
    SET_STRING_ELT(names, pos, mkChar("nn.index"));
    SET_VECTOR_ELT(out,   pos, out_idx);
    ++pos;
  }
  if(return_dist){
    SET_STRING_ELT(names, pos, mkChar("nn.dist"));
    SET_VECTOR_ELT(out,   pos, out_dst);
  }
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(protect_count);
  return out;
}

//' @title KNN bruto nativo
//' @description Executa KNN euclidiano exato com distância vetorizada e ordenação parcial.
//' @return Lista com índices e distâncias.
extern "C" SEXP brute_force_knn_c(SEXP data, SEXP query, SEXP k_value){
  return brute_force_knn_impl(data, query, k_value, 1, 1);
}

//' @title Índices KNN brutos nativos
//' @description Executa KNN euclidiano exato retornando apenas índices.
//' @return Lista com matriz de índices.
extern "C" SEXP brute_force_knn_index_c(SEXP data, SEXP query, SEXP k_value){
  return brute_force_knn_impl(data, query, k_value, 1, 0);
}

//' @title Distâncias KNN brutas nativas
//' @description Executa KNN euclidiano exato retornando apenas distâncias.
//' @return Lista com matriz de distâncias.
extern "C" SEXP brute_force_knn_dist_c(SEXP data, SEXP query, SEXP k_value){
  return brute_force_knn_impl(data, query, k_value, 0, 1);
}

//' @title KNN bruto nativo parametrizado
//' @description Executa KNN euclidiano exato interno com controle de retorno e self-neighbor.
//' @return Lista com matrizes KNN solicitadas.
extern "C" SEXP brute_force_knn_native_c(SEXP data, SEXP query, SEXP k_value,
                                           SEXP return_code, SEXP query_is_data_value,
                                           SEXP exclude_self_value, SEXP query_offset_value){
  const int ret              = asInteger(return_code);
  const R_xlen_t query_offset = (R_xlen_t) asInteger(query_offset_value);
  const int need_index       = ret == 0 || ret == 1;
  const int need_dist        = ret == 0 || ret == 2;
  if(ret < 0 || ret > 2){
    error("'return_code' deve ser 0 (both), 1 (index) ou 2 (dist)");
  }
  if(query_offset < 0){
    error("'query_offset' deve ser >= 0");
  }
  return brute_force_knn_impl(
    data, query, k_value,
    need_index, need_dist,
    asLogical(query_is_data_value),
    asLogical(exclude_self_value),
    query_offset
  );
}

//' @title NearMiss bruto fundido
//' @description Calcula distâncias KNN e seleciona exemplos majoritários sem materializar a matriz intermediária no R.
//' @return Vetor inteiro de índices retidos.
extern "C" SEXP nearmiss_brute_select_c(SEXP minority_data, SEXP majority_query,
                                         SEXP majority_index, SEXP k_value,
                                         SEXP retained_majority_count){
  require_real_matrix(minority_data,  "minority_data");
  require_real_matrix(majority_query, "majority_query");
  if(!isInteger(majority_index)){
    error("'majority_index' deve ser integer");
  }
  if(!isInteger(k_value) || LENGTH(k_value) != 1){
    error("'k' deve ser integer escalar");
  }
  if(!isInteger(retained_majority_count) || LENGTH(retained_majority_count) != 1){
    error("'retained_majority_count' deve ser integer escalar");
  }

  SEXP d_dims = getAttrib(minority_data,  R_DimSymbol);
  SEXP q_dims = getAttrib(majority_query, R_DimSymbol);
  const R_xlen_t n_ref   = (R_xlen_t) INTEGER(d_dims)[0];
  const R_xlen_t p1      = (R_xlen_t) INTEGER(d_dims)[1];
  const R_xlen_t n_query = (R_xlen_t) INTEGER(q_dims)[0];
  const R_xlen_t p2      = (R_xlen_t) INTEGER(q_dims)[1];
  const int k            = INTEGER(k_value)[0];
  const int retain       = INTEGER(retained_majority_count)[0];

  if(p1 != p2){
    error("'minority_data' e 'majority_query' devem ter o mesmo numero de colunas");
  }
  if(k < 1 || (R_xlen_t) k > n_ref){
    error("'k' deve estar entre 1 e nrow(minority_data)");
  }
  if(retain < 1 || (R_xlen_t) retain > n_query){
    error("'retained_majority_count' deve estar entre 1 e nrow(majority_query)");
  }
  if(XLENGTH(majority_index) != n_query){
    error("'majority_index' deve ter comprimento igual a nrow(majority_query)");
  }
  if(n_query > INT_MAX){
    error("'majority_query' excede o limite suportado por allocVector");
  }

  const double *r_ref   = REAL(minority_data);
  const double *q_ref   = REAL(majority_query);
  const int    *majority = INTEGER(majority_index);

  double *q_norm2 = (double *) R_alloc((size_t) n_query, sizeof(double));
  double *r_norm2 = (double *) R_alloc((size_t) n_ref,   sizeof(double));
  sby_compute_column_norms2(q_ref, q_norm2, n_query, p1);
  sby_compute_column_norms2(r_ref, r_norm2, n_ref,   p1);

  std::vector<nearmiss_score> scores;
  scores.reserve(static_cast<std::size_t>(n_query));

  const R_xlen_t max_cells = 8 * 1024 * 1024;
  R_xlen_t block_q = max_cells / (n_ref > 0 ? n_ref : 1);
  if(block_q < 1)       block_q = 1;
  if(block_q > n_query) block_q = n_query;
  if(block_q > INT_MAX) block_q = INT_MAX;
  double *cross = (double *) R_alloc((size_t) (block_q * n_ref), sizeof(double));

  for(R_xlen_t q_start = 0; q_start < n_query; q_start += block_q){
    R_xlen_t q_end = q_start + block_q;
    if(q_end > n_query) q_end = n_query;
    const R_xlen_t cur_b = q_end - q_start;
    R_CheckUserInterrupt();

    const double alpha    = 1.0;
    const double beta     = 0.0;
    const int m_int   = (int) cur_b;
    const int n_int   = (int) n_ref;
    const int kk_int  = (int) p1;
    const int lda_int = (int) n_query;
    const int ldb_int = (int) n_ref;
    const int ldc_int = (int) cur_b;
    F77_CALL(dgemm)(
      "N", "T",
      &m_int, &n_int, &kk_int,
      &alpha,
      q_ref + q_start, &lda_int,
      r_ref,           &ldb_int,
      &beta,
      cross,           &ldc_int FCONE FCONE
    );

    std::vector<knn_neighbor> neighbors(static_cast<std::size_t>(n_ref));
    for(R_xlen_t bi = 0; bi < cur_b; ++bi){
      const R_xlen_t q_idx = q_start + bi;
      const double   qn2   = q_norm2[q_idx];

      for(R_xlen_t j = 0; j < n_ref; ++j){
        double d2 = qn2 + r_norm2[j] - 2.0 * cross[bi + j * cur_b];
        d2 = euclidean_d2_with_rescue(d2, q_ref, r_ref, q_idx, j, n_query, n_ref, p1);
        neighbors[static_cast<std::size_t>(j)] = {d2, (int) (j + 1)};
      }

      auto first  = neighbors.begin();
      auto middle = neighbors.begin() + k;
      auto last   = neighbors.end();

      sby_select_topk(first, middle, last, knn_neighbor_less);
      sby_sort_topk(first, middle, knn_neighbor_less,
                    static_cast<std::size_t>(k),
                    static_cast<std::size_t>(n_ref));

      long double sum_dist = 0.0L;
#pragma omp simd reduction(+:sum_dist)
      for(int kk = 0; kk < k; ++kk){
        sum_dist += (long double) sqrt(neighbors[static_cast<std::size_t>(kk)].dist2);
      }
      scores.push_back({
        (double) (sum_dist / (long double) k),
        majority[q_idx]
      });
    }
  }

  /* sby_partial_sort mantido: ordena vetor de nearmiss_score (n_query elementos) */
  sby_partial_sort(scores.begin(), scores.begin() + retain, scores.end(), nearmiss_score_less);

  SEXP out = PROTECT(allocVector(INTSXP, retain));
  for(int i = 0; i < retain; ++i){
    INTEGER(out)[i] = scores[static_cast<std::size_t>(i)].index;
  }

  UNPROTECT(1);
  return out;
}

//' @title Concatenação nativa de matrizes double
//' @description Empilha duas matrizes double por linhas preservando nomes de colunas quando disponíveis.
//' @return Matriz double concatenada.
extern "C" SEXP rbind_double_matrix_c(SEXP first_matrix, SEXP second_matrix){
  require_real_matrix(first_matrix,  "first_matrix");
  require_real_matrix(second_matrix, "second_matrix");

  SEXP first_dims  = getAttrib(first_matrix,  R_DimSymbol);
  SEXP second_dims = getAttrib(second_matrix, R_DimSymbol);
  const R_xlen_t first_rows  = (R_xlen_t) INTEGER(first_dims)[0];
  const R_xlen_t second_rows = (R_xlen_t) INTEGER(second_dims)[0];
  const R_xlen_t first_cols  = (R_xlen_t) INTEGER(first_dims)[1];
  const R_xlen_t second_cols = (R_xlen_t) INTEGER(second_dims)[1];
  const R_xlen_t out_rows    = first_rows + second_rows;

  if(first_cols != second_cols){
    error("'first_matrix' e 'second_matrix' devem ter o mesmo numero de colunas");
  }
  if(out_rows > INT_MAX || first_cols > INT_MAX){
    error("Dimensoes excedem o limite suportado por matrix R");
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, (int) out_rows, (int) first_cols));
  const double *first    = REAL(first_matrix);
  const double *second   = REAL(second_matrix);
  double       *expanded = REAL(out);

  for(R_xlen_t j = 0; j < first_cols; ++j){
    if((j & 31) == 0){
      R_CheckUserInterrupt();
    }
    double *out_col = expanded + j * out_rows;
    memcpy(out_col,             first  + j * first_rows,  (size_t) first_rows  * sizeof(double));
    memcpy(out_col + first_rows, second + j * second_rows, (size_t) second_rows * sizeof(double));
  }

  SEXP first_dimnames = getAttrib(first_matrix, R_DimNamesSymbol);
  if(!isNull(first_dimnames) && XLENGTH(first_dimnames) >= 2){
    SEXP first_colnames = VECTOR_ELT(first_dimnames, 1);
    if(!isNull(first_colnames)){
      SEXP out_dimnames = PROTECT(allocVector(VECSXP, 2));
      SET_VECTOR_ELT(out_dimnames, 0, R_NilValue);
      SET_VECTOR_ELT(out_dimnames, 1, first_colnames);
      setAttrib(out, R_DimNamesSymbol, out_dimnames);
      UNPROTECT(1);
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP brute_force_knn_rcpp_parallel_c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP brute_force_knn_native_parallel_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP rcpp_parallel_uses_tbb_c(void);

extern "C" {
  void sby_compute_zscore_population_f(const double *x, int n, int p,
                                       double *means, double *sds, int *status);
  void sby_apply_zscore_f(const double *x, int n, int p,
                          const double *means, const double *sds,
                          double *x_out, int *status);
  void sby_revert_zscore_f(const double *x, int n, int p,
                           const double *means, const double *sds,
                           double *x_out, int *status);
  void sby_rbind_matrix_f(const double *a, int n1, int p,
                          const double *b, int n2,
                          double *c_out, int *status);
}

//' @title Z-Score populacional nativo Fortran
//' @description Calcula media e desvio padrao populacionais via motor Fortran AVX-512.
//' @param x_matrix Matriz double em layout R n x p.
//' @return Lista com vetores means e sds.
extern "C" SEXP compute_zscore_population_fortran_c(SEXP x_matrix){
  require_real_matrix(x_matrix, "x_matrix");
  SEXP dims   = getAttrib(x_matrix, R_DimSymbol);
  const int n = INTEGER(dims)[0];
  const int p = INTEGER(dims)[1];

  SEXP means_sexp = PROTECT(allocVector(REALSXP, p));
  SEXP sds_sexp   = PROTECT(allocVector(REALSXP, p));

  int status = 0;
  sby_compute_zscore_population_f(REAL(x_matrix), n, p,
                                  REAL(means_sexp), REAL(sds_sexp), &status);
  if(status != 0){
    UNPROTECT(2);
    error("Erro no calculo de z-score populacional nativo (status=%d)", status);
  }

  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("means"));
  SET_STRING_ELT(names, 1, mkChar("sds"));

  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, means_sexp);
  SET_VECTOR_ELT(out, 1, sds_sexp);
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(4);
  return out;
}

//' @title Aplicacao de Z-Score nativo Fortran
//' @description Aplica padronizacao Z-Score via motor Fortran AVX-512.
//' @param x_matrix Matriz double em layout R n x p.
//' @param means Vetor double de medias.
//' @param sds Vetor double de desvios padrao.
//' @return Matriz double padronizada.
extern "C" SEXP apply_zscore_fortran_c(SEXP x_matrix, SEXP means, SEXP sds){
  require_real_matrix(x_matrix, "x_matrix");
  SEXP dims   = getAttrib(x_matrix, R_DimSymbol);
  const int n = INTEGER(dims)[0];
  const int p = INTEGER(dims)[1];
  require_real_vector_length(means, p, "means");
  require_real_vector_length(sds, p, "sds");

  SEXP out   = PROTECT(allocMatrix(REALSXP, n, p));
  int status = 0;
  sby_apply_zscore_f(REAL(x_matrix), n, p, REAL(means), REAL(sds),
                     REAL(out), &status);
  if(status != 0){
    UNPROTECT(1);
    error("Erro na aplicacao de z-score nativo (status=%d)", status);
  }

  UNPROTECT(1);
  return out;
}

//' @title Reversao de Z-Score nativo Fortran
//' @description Reverte padronizacao Z-Score via motor Fortran AVX-512.
//' @param x_matrix Matriz double em layout R n x p.
//' @param means Vetor double de medias.
//' @param sds Vetor double de desvios padrao.
//' @return Matriz double na escala original.
extern "C" SEXP revert_zscore_fortran_c(SEXP x_matrix, SEXP means, SEXP sds){
  require_real_matrix(x_matrix, "x_matrix");
  SEXP dims   = getAttrib(x_matrix, R_DimSymbol);
  const int n = INTEGER(dims)[0];
  const int p = INTEGER(dims)[1];
  require_real_vector_length(means, p, "means");
  require_real_vector_length(sds, p, "sds");

  SEXP out   = PROTECT(allocMatrix(REALSXP, n, p));
  int status = 0;
  sby_revert_zscore_f(REAL(x_matrix), n, p, REAL(means), REAL(sds),
                      REAL(out), &status);
  if(status != 0){
    UNPROTECT(1);
    error("Erro na reversao de z-score nativo (status=%d)", status);
  }

  UNPROTECT(1);
  return out;
}

//' @title Concatenacao de matrizes nativa Fortran
//' @description Empilha duas matrizes n x p por linhas via motor Fortran com OpenMP.
//' @param a_matrix Primeira matriz double (N1 x P).
//' @param b_matrix Segunda matriz double (N2 x P).
//' @return Matriz double ((N1+N2) x P).
extern "C" SEXP rbind_matrix_fortran_c(SEXP a_matrix, SEXP b_matrix){
  require_real_matrix(a_matrix, "a_matrix");
  require_real_matrix(b_matrix, "b_matrix");

  SEXP a_dims    = getAttrib(a_matrix, R_DimSymbol);
  SEXP b_dims    = getAttrib(b_matrix, R_DimSymbol);
  const int n1   = INTEGER(a_dims)[0];
  const int p    = INTEGER(a_dims)[1];
  const int n2   = INTEGER(b_dims)[0];
  const int p2   = INTEGER(b_dims)[1];

  if(p != p2){
    error("As matrizes devem ter o mesmo numero de colunas (features)");
  }

  SEXP out   = PROTECT(allocMatrix(REALSXP, n1 + n2, p));
  int status = 0;
  sby_rbind_matrix_f(REAL(a_matrix), n1, p, REAL(b_matrix), n2,
                     REAL(out), &status);
  if(status != 0){
    UNPROTECT(1);
    error("Erro na concatenacao nativa (status=%d)", status);
  }

  UNPROTECT(1);
  return out;
}

static const R_CallMethodDef call_entries[] = {
  {"check_user_interrupt_c",              (DL_FUNC) &check_user_interrupt_c,              0},
  {"compute_z_score_params_c",            (DL_FUNC) &compute_z_score_params_c,            1},
  {"apply_z_score_c",                     (DL_FUNC) &apply_z_score_c,                     4},
  {"generate_synthetic_adasyn_c",         (DL_FUNC) &generate_synthetic_adasyn_c,         3},
  {"generate_synthetic_adasyn_col_c",     (DL_FUNC) &generate_synthetic_adasyn_col_c,     3},
  {"select_nearmiss_majority_c",          (DL_FUNC) &select_nearmiss_majority_c,          3},
  {"drop_self_neighbor_c",                (DL_FUNC) &drop_self_neighbor_c,                3},
  {"brute_force_knn_c",                   (DL_FUNC) &brute_force_knn_c,                   3},
  {"brute_force_knn_index_c",             (DL_FUNC) &brute_force_knn_index_c,             3},
  {"brute_force_knn_dist_c",              (DL_FUNC) &brute_force_knn_dist_c,              3},
  {"brute_force_knn_rcpp_parallel_c",     (DL_FUNC) &brute_force_knn_rcpp_parallel_c,     5},
  {"brute_force_knn_native_c",            (DL_FUNC) &brute_force_knn_native_c,            7},
  {"brute_force_knn_native_parallel_c",   (DL_FUNC) &brute_force_knn_native_parallel_c,   8},
  {"rcpp_parallel_uses_tbb_c",            (DL_FUNC) &rcpp_parallel_uses_tbb_c,            0},
  {"nearmiss_brute_select_c",             (DL_FUNC) &nearmiss_brute_select_c,             5},
  {"rbind_double_matrix_c",               (DL_FUNC) &rbind_double_matrix_c,               2},
  {"compute_zscore_population_fortran_c", (DL_FUNC) &compute_zscore_population_fortran_c, 1},
  {"apply_zscore_fortran_c",              (DL_FUNC) &apply_zscore_fortran_c,              3},
  {"revert_zscore_fortran_c",             (DL_FUNC) &revert_zscore_fortran_c,             3},
  {"rbind_matrix_fortran_c",              (DL_FUNC) &rbind_matrix_fortran_c,              2},
  {NULL, NULL, 0}
};

//' @title Registro nativo do pacote
//' @description Registra as rotinas .Call e desativa resolução dinâmica de símbolos.
//' @param dll Ponteiro para a biblioteca dinâmica carregada pelo R.
extern "C" void attribute_visible R_init_sbyadanear(DllInfo *dll){
  R_registerRoutines(dll, NULL, call_entries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
