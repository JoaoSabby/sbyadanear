#include <Rcpp.h>
#include <RcppParallel.h>
#include <Rinternals.h>
#include <R_ext/Error.h>
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

  sby_brute_force_knn_worker(
    Rcpp::NumericMatrix data_,
    Rcpp::NumericMatrix query_,
    Rcpp::IntegerMatrix index_out_,
    Rcpp::NumericMatrix dist_out_,
    int k_,
    bool need_index_,
    bool need_dist_
  ) : data(data_),
      query(query_),
      index_out(index_out_),
      dist_out(dist_out_),
      k(k_),
      need_index(need_index_),
      need_dist(need_dist_) {}

  void operator()(std::size_t begin, std::size_t end){
    const std::size_t n_ref = data.nrow();
    const std::size_t n_cols = data.ncol();

    for(std::size_t q = begin; q < end; ++q){
      std::vector<sby_neighbor> neighbors(n_ref);
      for(std::size_t r = 0; r < n_ref; ++r){
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

  const bool need_index = ret == 0 || ret == 1;
  const bool need_dist = ret == 0 || ret == 2;
  SEXP index_out = R_NilValue;
  SEXP dist_out = R_NilValue;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("nn.index"));
  SET_STRING_ELT(names, 1, Rf_mkChar("nn.dist"));

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
    need_dist
  );
  RcppParallel::parallelFor(
    0,
    static_cast<std::size_t>(n_query),
    worker,
    1,
    workers
  );

  if(need_index){
    SET_VECTOR_ELT(out, 0, index_out);
  }
  if(need_dist){
    SET_VECTOR_ELT(out, 1, dist_out);
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(4);
  return out;
}
