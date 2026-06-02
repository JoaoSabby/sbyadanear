#include <RcppParallel.h>
#include <Rinternals.h>
#include <R_ext/Error.h>
#include <cmath>
#include <cfloat>
#include <climits>
#include <vector>

namespace {

struct SbyNeighbor {
  double dist2;
  int index;
};

inline bool sby_worse_neighbor(const SbyNeighbor& lhs, const SbyNeighbor& rhs){
  return lhs.dist2 > rhs.dist2 || (lhs.dist2 == rhs.dist2 && lhs.index > rhs.index);
}

inline bool sby_better_neighbor(const double dist2, const int index, const SbyNeighbor& rhs){
  return dist2 < rhs.dist2 || (dist2 == rhs.dist2 && index < rhs.index);
}

inline void sby_heap_sift_up(std::vector<SbyNeighbor>& heap, int pos){
  while(pos > 0){
    const int parent = (pos - 1) / 2;
    if(!sby_worse_neighbor(heap[pos], heap[parent])){
      break;
    }
    std::swap(heap[pos], heap[parent]);
    pos = parent;
  }
}

inline void sby_heap_sift_down(std::vector<SbyNeighbor>& heap, int heap_size, int pos){
  while(true){
    const int left = 2 * pos + 1;
    const int right = left + 1;
    int worst = pos;
    if(left < heap_size && sby_worse_neighbor(heap[left], heap[worst])){
      worst = left;
    }
    if(right < heap_size && sby_worse_neighbor(heap[right], heap[worst])){
      worst = right;
    }
    if(worst == pos){
      break;
    }
    std::swap(heap[pos], heap[worst]);
    pos = worst;
  }
}

inline void sby_sort_neighbors(std::vector<SbyNeighbor>& heap, int heap_size){
  for(int i = 1; i < heap_size; ++i){
    SbyNeighbor key = heap[i];
    int j = i - 1;
    while(j >= 0 && sby_worse_neighbor(heap[j], key)){
      heap[j + 1] = heap[j];
      --j;
    }
    heap[j + 1] = key;
  }
}

struct SbyBruteForceKnnWorker : public RcppParallel::Worker {
  RcppParallel::RMatrix<double> data;
  RcppParallel::RMatrix<double> query;
  RcppParallel::RMatrix<int> index_out;
  RcppParallel::RMatrix<double> dist_out;
  const int k;
  const bool need_index;
  const bool need_dist;

  SbyBruteForceKnnWorker(
    SEXP data_,
    SEXP query_,
    SEXP index_out_,
    SEXP dist_out_,
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
    const std::size_t n_query = query.nrow();
    const std::size_t n_cols = data.ncol();

    std::vector<SbyNeighbor> heap(static_cast<std::size_t>(k));

    for(std::size_t q = begin; q < end; ++q){
      int heap_size = 0;
      for(std::size_t r = 0; r < n_ref; ++r){
        long double dist2_ld = 0.0L;
        for(std::size_t c = 0; c < n_cols; ++c){
          const long double diff = static_cast<long double>(query(q, c)) -
            static_cast<long double>(data(r, c));
          dist2_ld += diff * diff;
        }
        double dist2 = static_cast<double>(dist2_ld);
        if(dist2 < 0.0){
          dist2 = 0.0;
        }
        const int one_based_index = static_cast<int>(r + 1);

        if(heap_size < k){
          heap[static_cast<std::size_t>(heap_size)] = {dist2, one_based_index};
          sby_heap_sift_up(heap, heap_size);
          ++heap_size;
        }else if(sby_better_neighbor(dist2, one_based_index, heap[0])){
          heap[0] = {dist2, one_based_index};
          sby_heap_sift_down(heap, heap_size, 0);
        }
      }

      sby_sort_neighbors(heap, heap_size);
      for(int j = 0; j < k; ++j){
        if(need_index){
          index_out(q, static_cast<std::size_t>(j)) = heap[static_cast<std::size_t>(j)].index;
        }
        if(need_dist){
          dist_out(q, static_cast<std::size_t>(j)) = std::sqrt(heap[static_cast<std::size_t>(j)].dist2);
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

extern "C" SEXP OU_RcppParallelUsesTbbC(){
#if RCPP_PARALLEL_USE_TBB
  return Rf_ScalarLogical(1);
#else
  return Rf_ScalarLogical(0);
#endif
}

extern "C" SEXP OU_BruteForceKnnRcppParallelC(
  SEXP dataMatrix,
  SEXP queryMatrix,
  SEXP kValue,
  SEXP returnCode,
  SEXP workersValue
){
  sby_require_real_matrix(dataMatrix, "dataMatrix");
  sby_require_real_matrix(queryMatrix, "queryMatrix");

  SEXP dataDims = Rf_getAttrib(dataMatrix, R_DimSymbol);
  SEXP queryDims = Rf_getAttrib(queryMatrix, R_DimSymbol);
  const int nRef = INTEGER(dataDims)[0];
  const int dataCols = INTEGER(dataDims)[1];
  const int nQuery = INTEGER(queryDims)[0];
  const int queryCols = INTEGER(queryDims)[1];
  const int k = Rf_asInteger(kValue);
  const int ret = Rf_asInteger(returnCode);
  int workers = Rf_asInteger(workersValue);

  if(dataCols != queryCols){
    Rf_error("'dataMatrix' e 'queryMatrix' devem ter o mesmo numero de colunas");
  }
  if(k < 1 || k > nRef){
    Rf_error("'kValue' deve estar entre 1 e nrow(dataMatrix)");
  }
  if(workers < 1){
    workers = 1;
  }

  const bool needIndex = ret == 0 || ret == 1;
  const bool needDist = ret == 0 || ret == 2;
  SEXP indexOut = R_NilValue;
  SEXP distOut = R_NilValue;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("nn.index"));
  SET_STRING_ELT(names, 1, Rf_mkChar("nn.dist"));

  if(needIndex){
    indexOut = PROTECT(Rf_allocMatrix(INTSXP, nQuery, k));
  }else{
    indexOut = PROTECT(Rf_allocMatrix(INTSXP, 0, 0));
  }
  if(needDist){
    distOut = PROTECT(Rf_allocMatrix(REALSXP, nQuery, k));
  }else{
    distOut = PROTECT(Rf_allocMatrix(REALSXP, 0, 0));
  }

  SbyBruteForceKnnWorker worker(dataMatrix, queryMatrix, indexOut, distOut, k, needIndex, needDist);
  RcppParallel::parallelFor(
    0,
    static_cast<std::size_t>(nQuery),
    worker,
    1,
    workers
  );

  if(needIndex){
    SET_VECTOR_ELT(out, 0, indexOut);
  }
  if(needDist){
    SET_VECTOR_ELT(out, 1, distOut);
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(4);
  return out;
}
