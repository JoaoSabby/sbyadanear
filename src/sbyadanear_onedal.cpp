#include <Rinternals.h>
#include <R_ext/Error.h>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <vector>

#if defined(SBYADANEAR_ONEDAL)
#include "oneapi/dal/algo/knn.hpp"
#include "oneapi/dal/table/homogen.hpp"
#include "oneapi/dal/table/row_accessor.hpp"
#include <tbb/global_control.h>
#endif

extern "C" SEXP OU_BruteForceKnnRcppParallelC(
  SEXP dataMatrix,
  SEXP queryMatrix,
  SEXP kValue,
  SEXP returnCode,
  SEXP workersValue
);

namespace {

void sby_require_real_matrix(SEXP x, const char *name){
  if(TYPEOF(x) != REALSXP || !Rf_isMatrix(x)){
    Rf_error("O parametro '%s' deve ser uma matrix double", name);
  }
}

SEXP sby_named_knn_result(SEXP index_out, SEXP dist_out){
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("nn.index"));
  SET_STRING_ELT(names, 1, Rf_mkChar("nn.dist"));
  SET_VECTOR_ELT(out, 0, index_out);
  SET_VECTOR_ELT(out, 1, dist_out);
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

#if defined(SBYADANEAR_ONEDAL)

void sby_copy_r_matrix_to_row_major(
  const double *src,
  double *dst,
  const int n_rows,
  const int n_cols
){
  for(int i = 0; i < n_rows; ++i){
    for(int j = 0; j < n_cols; ++j){
      dst[static_cast<std::size_t>(i) * static_cast<std::size_t>(n_cols) +
        static_cast<std::size_t>(j)] =
        src[static_cast<std::size_t>(i) +
          static_cast<std::size_t>(j) * static_cast<std::size_t>(n_rows)];
    }
  }
}

#endif

} // namespace

//' @title Rotina de busca KNN acelerada via Intel oneDAL
//' @description
//' Implementa a busca de k-vizinhos mais proximos utilizando o backend
//' Intel oneDAL. Delega o processamento pesado para as instrucoes
//' vetorizadas do hardware e utiliza TBB para paralelismo.
//' @param data_matrix Matriz de dados de referencia (SEXP matriz double).
//' @param query_matrix Matriz de dados de consulta (SEXP matriz double).
//' @param k_value Numero inteiro de vizinhos a buscar.
//' @param workers_value Quantidade de threads alocadas para o TBB.
//' @return Uma lista R com duas matrizes nn.index e nn.dist.
extern "C" SEXP run_one_dal_knn_c(
  SEXP data_matrix,
  SEXP query_matrix,
  SEXP k_value,
  SEXP workers_value
) {
  // valida entradas R antes de delegar para o backend numerico
  sby_require_real_matrix(data_matrix, "data_matrix");
  sby_require_real_matrix(query_matrix, "query_matrix");

  // extrai dimensoes respeitando matrizes column_major do R
  SEXP data_dims = Rf_getAttrib(data_matrix, R_DimSymbol);
  SEXP query_dims = Rf_getAttrib(query_matrix, R_DimSymbol);
  const int n_ref = INTEGER(data_dims)[0];
  const int data_cols = INTEGER(data_dims)[1];
  const int n_query = INTEGER(query_dims)[0];
  const int query_cols = INTEGER(query_dims)[1];
  const int k = Rf_asInteger(k_value);
  int workers = Rf_asInteger(workers_value);

  // bloqueia combinacoes invalidas antes de alocar saidas
  if(data_cols != query_cols){
    Rf_error("'data_matrix' e 'query_matrix' devem ter o mesmo numero de colunas");
  }
  if(k < 1 || k > n_ref){
    Rf_error("'k_value' deve estar entre 1 e nrow(data_matrix)");
  }
  if(workers < 1){
    workers = 1;
  }

#if defined(SBYADANEAR_ONEDAL)
  // copia dados para buffers row_major consumidos diretamente pelo oneDAL
  const std::size_t data_size = static_cast<std::size_t>(n_ref) *
    static_cast<std::size_t>(data_cols);
  const std::size_t query_size = static_cast<std::size_t>(n_query) *
    static_cast<std::size_t>(query_cols);
  std::vector<double> data_row_major(data_size);
  std::vector<double> query_row_major(query_size);
  sby_copy_r_matrix_to_row_major(REAL(data_matrix), data_row_major.data(), n_ref, data_cols);
  sby_copy_r_matrix_to_row_major(REAL(query_matrix), query_row_major.data(), n_query, query_cols);

  namespace dal = oneapi::dal;
  namespace knn = dal::knn;

  // limita paralelismo TBB conforme solicitado pelo chamador R
  tbb::global_control tbb_control(
    tbb::global_control::max_allowed_parallelism,
    static_cast<std::size_t>(workers)
  );

  // executa busca brute_force exata com distancia euclidiana padrao do oneDAL
  const auto data_table = dal::homogen_table::wrap(data_row_major.data(), n_ref, data_cols);
  const auto query_table = dal::homogen_table::wrap(query_row_major.data(), n_query, query_cols);
  const auto knn_desc = knn::descriptor<double, knn::method::brute_force, knn::task::search>(k)
    .set_result_options(knn::result_options::indices | knn::result_options::distances);
  const auto train_result = dal::train(knn_desc, data_table);
  const auto infer_result = dal::infer(knn_desc, query_table, train_result.get_model());

  // materializa tabelas oneDAL em matrizes R com indices um_based
  SEXP index_out = PROTECT(Rf_allocMatrix(INTSXP, n_query, k));
  SEXP dist_out = PROTECT(Rf_allocMatrix(REALSXP, n_query, k));
  const auto indices_block = dal::row_accessor<const std::int64_t>(infer_result.get_indices()).pull();
  const auto distances_block = dal::row_accessor<const double>(infer_result.get_distances()).pull();
  const std::int64_t *indices = indices_block.get_data();
  const double *distances = distances_block.get_data();
  int *index_ptr = INTEGER(index_out);
  double *dist_ptr = REAL(dist_out);

  for(int q = 0; q < n_query; ++q){
    for(int j = 0; j < k; ++j){
      const std::size_t src_pos = static_cast<std::size_t>(q) * static_cast<std::size_t>(k) +
        static_cast<std::size_t>(j);
      const std::size_t dst_pos = static_cast<std::size_t>(q) +
        static_cast<std::size_t>(j) * static_cast<std::size_t>(n_query);
      if(indices[src_pos] < 0 || indices[src_pos] >= n_ref || indices[src_pos] >= INT_MAX){
        Rf_error("oneDAL retornou indice de vizinho fora do intervalo suportado");
      }
      index_ptr[dst_pos] = static_cast<int>(indices[src_pos]) + 1;
      dist_ptr[dst_pos] = distances[src_pos];
    }
  }

  SEXP out = PROTECT(sby_named_knn_result(index_out, dist_out));
  UNPROTECT(3);
  return out;
#else
  // usa fallback nativo paralelo quando oneDAL nao foi habilitado na compilacao
  SEXP return_code = PROTECT(Rf_ScalarInteger(0));
  SEXP out = PROTECT(OU_BruteForceKnnRcppParallelC(
    data_matrix,
    query_matrix,
    k_value,
    return_code,
    workers_value
  ));
  UNPROTECT(2);
  return out;
#endif
}

extern "C" SEXP OU_OneDalAvailableC(){
#if defined(SBYADANEAR_ONEDAL)
  return Rf_ScalarLogical(1);
#else
  return Rf_ScalarLogical(0);
#endif
}
