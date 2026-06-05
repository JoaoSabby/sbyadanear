#include <Rcpp.h>
#include <RcppParallel.h>
#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Constants.h>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <climits>
#include <vector>

namespace {

struct sby_neighbor {
  double dist2;
  int index;
};

inline bool sby_neighbor_less(const sby_neighbor& lhs, const sby_neighbor& rhs){
  return lhs.dist2 < rhs.dist2 || (lhs.dist2 == rhs.dist2 && lhs.index < rhs.index);
}

struct sby_brute_force_knn_worker : public RcppParallel::Worker {
  RcppParallel::RMatrix<double> data;
  RcppParallel::RMatrix<double> query;
  RcppParallel::RMatrix<int> index_out;
  RcppParallel::RMatrix<double> dist_out;
  const int k;
  const bool need_index;
  const bool need_dist;
  const bool query_is_data;
  const bool exclude_self;
  const std::size_t query_offset;

  sby_brute_force_knn_worker(
    Rcpp::NumericMatrix data_,
    Rcpp::NumericMatrix query_,
    Rcpp::IntegerMatrix index_out_,
    Rcpp::NumericMatrix dist_out_,
    int k_,
    bool need_index_,
    bool need_dist_,
    bool query_is_data_,
    bool exclude_self_,
    std::size_t query_offset_
  ) : data(data_),
      query(query_),
      index_out(index_out_),
      dist_out(dist_out_),
      k(k_),
      need_index(need_index_),
      need_dist(need_dist_),
      query_is_data(query_is_data_),
      exclude_self(exclude_self_),
      query_offset(query_offset_) {}

  void operator()(std::size_t begin, std::size_t end){
    const std::size_t n_ref = data.nrow();
    const std::size_t n_cols = data.ncol();

    std::vector<sby_neighbor> neighbors(n_ref);
    for(std::size_t q = begin; q < end; ++q){
      for(std::size_t r = 0; r < n_ref; ++r){
        if(exclude_self && query_is_data && (query_offset + q) == r){
          neighbors[r] = {R_PosInf, static_cast<int>(r + 1)};
          continue;
        }
        long double dist2_ld = 0.0L;
#pragma omp simd reduction(+:dist2_ld)
        for(std::size_t c = 0; c < n_cols; ++c){
          const long double diff = static_cast<long double>(query(q, c)) -
            static_cast<long double>(data(r, c));
          dist2_ld += diff * diff;
        }
        double dist2 = static_cast<double>(dist2_ld);
        if(dist2 < 0.0){
          dist2 = 0.0;
        }
        neighbors[r] = {dist2, static_cast<int>(r + 1)};
      }

      std::partial_sort(
        neighbors.begin(),
        neighbors.begin() + k,
        neighbors.end(),
        sby_neighbor_less
      );
      for(int j = 0; j < k; ++j){
        const sby_neighbor& neighbor = neighbors[static_cast<std::size_t>(j)];
        if(need_index){
          index_out(q, static_cast<std::size_t>(j)) = neighbor.index;
        }
        if(need_dist){
          dist_out(q, static_cast<std::size_t>(j)) = std::sqrt(neighbor.dist2);
        }
      }
    }
  }
};

void sby_require_real_matrix(SEXP x, const char *name){
  if(TYPEOF(x) != REALSXP || !Rf_isMatrix(x)){
    Rf_error("O parametro '%s' deve ser uma matrix double", name);
  }
}

void sby_require_finite_matrix(SEXP x, const char *name){
  const double *values = REAL(x);
  const R_xlen_t n = Rf_xlength(x);
  for(R_xlen_t i = 0; i < n; ++i){
    if(!std::isfinite(values[i])){
      Rf_error("O parametro '%s' nao pode conter NA, NaN, Inf ou -Inf", name);
    }
  }
}

} // namespace

//' @title Detecção de backend TBB em RcppParallel
//' @description Informa se o RcppParallel foi compilado com suporte TBB.
//' @return Valor lógico indicado pelo backend nativo.
extern "C" SEXP rcpp_parallel_uses_tbb_c(){
#if RCPP_PARALLEL_USE_TBB
  return Rf_ScalarLogical(1);
#else
  return Rf_ScalarLogical(0);
#endif
}

//' @title KNN bruto paralelo com RcppParallel
//' @description Executa KNN euclidiano exato distribuindo consultas por trabalhadores nativos.
//' @return Lista R com índices e distâncias conforme solicitado.
extern "C" SEXP brute_force_knn_rcpp_parallel_c(
  SEXP data_matrix,
  SEXP query_matrix,
  SEXP k_value,
  SEXP return_code,
  SEXP workers_value
){
  sby_require_real_matrix(data_matrix, "data_matrix");
  sby_require_real_matrix(query_matrix, "query_matrix");

  SEXP data_dims = Rf_getAttrib(data_matrix, R_DimSymbol);
  SEXP query_dims = Rf_getAttrib(query_matrix, R_DimSymbol);
  const int n_ref = INTEGER(data_dims)[0];
  const int data_cols = INTEGER(data_dims)[1];
  const int n_query = INTEGER(query_dims)[0];
  const int query_cols = INTEGER(query_dims)[1];
  const int k = Rf_asInteger(k_value);
  const int ret = Rf_asInteger(return_code);
  int workers = Rf_asInteger(workers_value);

  if(data_cols != query_cols){
    Rf_error("'data_matrix' e 'query_matrix' devem ter o mesmo numero de colunas");
  }
  if(k < 1 || k > n_ref){
    Rf_error("'k_value' deve estar entre 1 e nrow(data_matrix)");
  }
  if(workers < 1){
    workers = 1;
  }
  if(ret < 0 || ret > 2){
    Rf_error("'return_code' deve ser 0 (both), 1 (index) ou 2 (dist)");
  }

  const bool need_index = ret == 0 || ret == 1;
  const bool need_dist = ret == 0 || ret == 2;
  SEXP index_out = R_NilValue;
  SEXP dist_out = R_NilValue;

  if(need_index){
    index_out = PROTECT(Rf_allocMatrix(INTSXP, n_query, k));
  }else{
    index_out = PROTECT(Rf_allocMatrix(INTSXP, 0, 0));
  }
  if(need_dist){
    dist_out = PROTECT(Rf_allocMatrix(REALSXP, n_query, k));
  }else{
    dist_out = PROTECT(Rf_allocMatrix(REALSXP, 0, 0));
  }

  Rcpp::NumericMatrix data_rcpp(data_matrix);
  Rcpp::NumericMatrix query_rcpp(query_matrix);
  Rcpp::IntegerMatrix index_out_rcpp(index_out);
  Rcpp::NumericMatrix dist_out_rcpp(dist_out);

  sby_brute_force_knn_worker worker(
    data_rcpp,
    query_rcpp,
    index_out_rcpp,
    dist_out_rcpp,
    k,
    need_index,
    need_dist,
    false,
    false,
    0
  );
  RcppParallel::parallelFor(
    0,
    static_cast<std::size_t>(n_query),
    worker,
    1,
    workers
  );

  const int out_length = need_index && need_dist ? 2 : 1;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, out_length));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, out_length));
  int pos = 0;
  if(need_index){
    SET_STRING_ELT(names, pos, Rf_mkChar("nn.index"));
    SET_VECTOR_ELT(out, pos, index_out);
    ++pos;
  }
  if(need_dist){
    SET_STRING_ELT(names, pos, Rf_mkChar("nn.dist"));
    SET_VECTOR_ELT(out, pos, dist_out);
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(4);
  return out;
}

//' @title KNN bruto nativo paralelo parametrizado
//' @description Executa KNN euclidiano exato distribuindo consultas por trabalhadores nativos, com controle de retorno e self-neighbor.
//' @return Lista R com índices e distâncias conforme solicitado.
extern "C" SEXP brute_force_knn_native_parallel_c(
  SEXP data_matrix,
  SEXP query_matrix,
  SEXP k_value,
  SEXP return_code,
  SEXP workers_value,
  SEXP query_is_data_value,
  SEXP exclude_self_value,
  SEXP query_offset_value
){
  sby_require_real_matrix(data_matrix, "data_matrix");
  sby_require_real_matrix(query_matrix, "query_matrix");
  sby_require_finite_matrix(data_matrix, "data_matrix");
  sby_require_finite_matrix(query_matrix, "query_matrix");

  SEXP data_dims = Rf_getAttrib(data_matrix, R_DimSymbol);
  SEXP query_dims = Rf_getAttrib(query_matrix, R_DimSymbol);
  const int n_ref = INTEGER(data_dims)[0];
  const int data_cols = INTEGER(data_dims)[1];
  const int n_query = INTEGER(query_dims)[0];
  const int query_cols = INTEGER(query_dims)[1];
  const int k = Rf_asInteger(k_value);
  const int ret = Rf_asInteger(return_code);
  int workers = Rf_asInteger(workers_value);
  const bool query_is_data = Rf_asLogical(query_is_data_value) == TRUE;
  const bool exclude_self = Rf_asLogical(exclude_self_value) == TRUE;
  const int query_offset_int = Rf_asInteger(query_offset_value);
  if(query_offset_int < 0){
    Rf_error("'query_offset' deve ser >= 0");
  }
  const std::size_t query_offset = static_cast<std::size_t>(query_offset_int);

  if(data_cols != query_cols){
    Rf_error("'data_matrix' e 'query_matrix' devem ter o mesmo numero de colunas");
  }
  if(exclude_self){
    if(!query_is_data){
      Rf_error("'exclude_self' requer 'query_is_data = TRUE'");
    }
    if(k < 1 || k > n_ref - 1){
      Rf_error("'k_value' deve estar entre 1 e nrow(data_matrix) - 1 quando 'exclude_self = TRUE'");
    }
  }else if(k < 1 || k > n_ref){
    Rf_error("'k_value' deve estar entre 1 e nrow(data_matrix)");
  }
  if(workers < 1){
    workers = 1;
  }
  if(ret < 0 || ret > 2){
    Rf_error("'return_code' deve ser 0 (both), 1 (index) ou 2 (dist)");
  }

  const bool need_index = ret == 0 || ret == 1;
  const bool need_dist = ret == 0 || ret == 2;
  SEXP index_out = R_NilValue;
  SEXP dist_out = R_NilValue;

  if(need_index){
    index_out = PROTECT(Rf_allocMatrix(INTSXP, n_query, k));
  }else{
    index_out = PROTECT(Rf_allocMatrix(INTSXP, 0, 0));
  }
  if(need_dist){
    dist_out = PROTECT(Rf_allocMatrix(REALSXP, n_query, k));
  }else{
    dist_out = PROTECT(Rf_allocMatrix(REALSXP, 0, 0));
  }

  Rcpp::NumericMatrix data_rcpp(data_matrix);
  Rcpp::NumericMatrix query_rcpp(query_matrix);
  Rcpp::IntegerMatrix index_out_rcpp(index_out);
  Rcpp::NumericMatrix dist_out_rcpp(dist_out);

  sby_brute_force_knn_worker worker(
    data_rcpp,
    query_rcpp,
    index_out_rcpp,
    dist_out_rcpp,
    k,
    need_index,
    need_dist,
    query_is_data,
    exclude_self,
    query_offset
  );
  RcppParallel::parallelFor(
    0,
    static_cast<std::size_t>(n_query),
    worker,
    1,
    workers
  );

  const int out_length = need_index && need_dist ? 2 : 1;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, out_length));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, out_length));
  int pos = 0;
  if(need_index){
    SET_STRING_ELT(names, pos, Rf_mkChar("nn.index"));
    SET_VECTOR_ELT(out, pos, index_out);
    ++pos;
  }
  if(need_dist){
    SET_STRING_ELT(names, pos, Rf_mkChar("nn.dist"));
    SET_VECTOR_ELT(out, pos, dist_out);
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(4);
  return out;
}
