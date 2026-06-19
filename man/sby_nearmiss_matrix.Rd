\name{sby_nearmiss_matrix}
\alias{sby_nearmiss_matrix}
\title{Aplicar NearMiss-1 sobre matrix double densa}
\usage{
sby_nearmiss_matrix(sby_x_matrix, sby_y_vector, sby_under_ratio = 1,
  sby_knn_under_k = 5L, sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE, sby_audit_level = c("none", "light", "full"),
  sby_return_index = TRUE, sby_return_scaled = FALSE,
  sby_return_original_scale = TRUE, sby_scaling_info = NULL,
  sby_input_already_scaled = FALSE, sby_fixed_minority_label = NULL,
  sby_fixed_majority_label = NULL,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "native", "FNN", "RcppHNSW", "KernelKnn", "bigKNN"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L, sby_knn_parallel_backend = c("parallel", "RcppParallel"), sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L, sby_knn_query_chunk_size = 1000L, sby_memory_guard = TRUE)
}
\arguments{
\item{sby_x_matrix}{Matrix double densa de preditores.}
\item{sby_y_vector}{Fator binário alinhado às linhas.}
\item{sby_under_ratio}{Fração alvo para retenção da maioria em relação à minoria.}
\item{sby_knn_under_k}{Número de vizinhos minoritários usados no NearMiss-1.}
\item{sby_seed}{Semente reprodutível.}
\item{sby_audit, sby_audit_level}{Controle de auditoria.}
\item{sby_return_index}{Retorna índices retidos.}
\item{sby_return_scaled}{Retorna matriz final em z-score.}
\item{sby_return_original_scale}{Reverte a saída para escala original.}
\item{sby_scaling_info, sby_input_already_scaled}{Permitem reutilizar z-score já aplicado.}
\item{sby_fixed_minority_label, sby_fixed_majority_label}{Papéis de classe fixos usados por ADANEAR.}
\item{sby_knn_algorithm, sby_knn_engine, sby_knn_distance_metric, sby_knn_workers, sby_knn_parallel_backend, sby_knn_hnsw_m, sby_knn_hnsw_ef, sby_knn_query_chunk_size}{Configuração KNN.}
\item{sby_memory_guard}{Ativa proteções de memória.}
}
\description{Implementa NearMiss-1, retendo majoritários mais próximos da minoria.}
\details{Quando \code{sby_knn_parallel_backend = "RcppParallel"}, o runtime concreto pode ser TBB/oneTBB ou TinyThread conforme a instalação de \pkg{RcppParallel}; os diagnósticos auditáveis incluem \code{sby_knn_parallel_runtime}. NearMiss-1 não remove majoritários por redundância interna da classe majoritária. A função reaproveita a matriz reduzida escalada calculada por \code{sby_nearmiss_index()} para evitar recálculo de z-score.}
\value{Lista com matriz e alvo reduzidos, distribuições, razões, diagnósticos e opcionalmente índices/escala.}
\examples{
x <- matrix(rnorm(40), ncol = 2); storage.mode(x) <- "double"
y <- factor(c(rep("minor", 5), rep("major", 15)))
sby_nearmiss_matrix(x, y, sby_seed = 1, sby_knn_engine = "FNN")
}
\seealso{\code{\link{sby_nearmiss_index}}, \code{\link{sby_adasyn_matrix}}}
