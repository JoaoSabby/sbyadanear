\name{sby_adanear_matrix}
\alias{sby_adanear_matrix}
\title{Aplicar ADASYN seguido de NearMiss-1}
\usage{
sby_adanear_matrix(sby_x_matrix, sby_y_vector, sby_over_ratio = 0.2,
  sby_under_ratio = 0.5, sby_knn_over_k = 5L, sby_knn_under_k = 5L,
  sby_seed = sample.int(10L^5L, 1L), sby_audit = FALSE,
  sby_audit_level = c("none", "light", "full"), sby_return_scaled = FALSE,
  sby_return_original_scale = TRUE,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "native", "FNN", "RcppHNSW", "KernelKnn", "bigKNN"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L, sby_knn_parallel_backend = c("parallel", "RcppParallel"), sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L, sby_knn_query_chunk_size = 1000L, sby_memory_guard = TRUE,
  sby_max_output_rows = Inf, sby_max_dense_gb = Inf)
}
\arguments{
\item{sby_x_matrix}{Matrix double densa de preditores.}
\item{sby_y_vector}{Fator binário.}
\item{sby_over_ratio, sby_under_ratio}{Controles de expansão ADASYN e retenção NearMiss.}
\item{sby_knn_over_k, sby_knn_under_k}{Vizinhos usados nas duas etapas.}
\item{sby_seed}{Semente reprodutível.}
\item{sby_audit, sby_audit_level}{Auditoria none/light/full.}
\item{sby_return_scaled, sby_return_original_scale}{Controle de escala retornada.}
\item{sby_knn_algorithm, sby_knn_engine, sby_knn_distance_metric, sby_knn_workers, sby_knn_parallel_backend, sby_knn_hnsw_m, sby_knn_hnsw_ef, sby_knn_query_chunk_size}{Configuração KNN.}
\item{sby_memory_guard, sby_max_output_rows, sby_max_dense_gb}{Controles de memória.}
}
\description{Aplica sobreamostragem ADASYN e depois subamostragem NearMiss-1 reutilizando a mesma escala.}
\details{Quando \code{sby_knn_parallel_backend = "RcppParallel"}, o runtime concreto pode ser TBB/oneTBB ou TinyThread conforme a instalação de \pkg{RcppParallel}; os diagnósticos auditáveis incluem \code{sby_knn_parallel_runtime}. Os papéis originais de minoria e maioria são preservados entre etapas. Auditoria padrão não guarda objetos intermediários pesados; \code{full} preserva artefatos completos.}
\value{Lista com matriz, alvo, razões, distribuições e diagnósticos; opcionalmente escala, índices e intermediários.}
\examples{
x <- matrix(rnorm(60), ncol = 2); storage.mode(x) <- "double"
y <- factor(c(rep("minor", 8), rep("major", 22)))
sby_adanear_matrix(x, y, sby_seed = 1, sby_knn_engine = "FNN")
}
\seealso{\code{\link{sby_adasyn_matrix}}, \code{\link{sby_nearmiss_matrix}}}
