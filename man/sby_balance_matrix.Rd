\name{sby_balance_matrix}
\alias{sby_balance_matrix}
\title{Roteador matricial industrial de estratégias de balanceamento}
\usage{
sby_balance_matrix(sby_x_matrix, sby_y_vector,
  sby_strategy = c("none", "weight", "adasyn", "nearmiss", "adanear", "adanearWeight"),
  sby_over_ratio = 0.2, sby_under_ratio = 0.5, sby_knn_over_k = 5L,
  sby_knn_under_k = 5L, sby_seed = sample.int(10L^5L, 1L),
  sby_audit = FALSE, sby_audit_level = c("none", "light", "full"),
  sby_return_index = TRUE, sby_return_scaled = FALSE,
  sby_return_original_scale = TRUE,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "FNN", "RcppHNSW", "KernelKnn", "bigKNN"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L, sby_knn_parallel_backend = c("parallel", "RcppParallel"), sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L, sby_knn_query_chunk_size = 1000L, sby_memory_guard = TRUE,
  sby_max_output_rows = Inf, sby_max_dense_gb = Inf)
}
\arguments{
\item{sby_x_matrix}{Objeto matricial. Para \code{none} e \code{weight} pode ser esparso; para balanceamento físico deve ser matrix double densa.}
\item{sby_y_vector}{Fator binário alinhado às linhas.}
\item{sby_strategy}{Estratégia: \code{none}, \code{weight}, \code{adasyn}, \code{nearmiss}, \code{adanear} ou \code{adanearWeight}.}
\item{sby_over_ratio, sby_under_ratio}{Controles de sobreamostragem e subamostragem.}
\item{sby_knn_over_k, sby_knn_under_k}{Números de vizinhos por etapa.}
\item{sby_seed}{Semente reprodutível.}
\item{sby_audit, sby_audit_level}{Auditoria none/light/full.}
\item{sby_return_index, sby_return_scaled, sby_return_original_scale}{Componentes e escala retornados.}
\item{sby_knn_algorithm, sby_knn_engine, sby_knn_distance_metric, sby_knn_workers, sby_knn_parallel_backend, sby_knn_hnsw_m, sby_knn_hnsw_ef, sby_knn_query_chunk_size}{Configuração KNN.}
\item{sby_memory_guard, sby_max_output_rows, sby_max_dense_gb}{Controles de memória.}
}
\description{Seleciona a estratégia de balanceamento adequada para pipelines industriais.}
\details{Quando \code{sby_knn_parallel_backend = "RcppParallel"}, o runtime concreto pode ser TBB/oneTBB ou TinyThread conforme a instalação de \pkg{RcppParallel}; os diagnósticos auditáveis incluem \code{sby_knn_parallel_runtime}. \code{none} e \code{weight} não alteram dados nem densificam matrizes. \code{adasyn}, \code{nearmiss}, \code{adanear} e \code{adanearWeight} fazem balanceamento físico e exigem matrix double densa. Em \code{adanearWeight}, use \code{sby_class_ratio_output} para \code{scale_pos_weight} pós-balanceamento.}
\value{Lista padronizada com \code{sby_x_matrix}, \code{sby_y_vector}, razões/distribuições e diagnósticos.}
\examples{
x <- matrix(rnorm(40), ncol = 2); storage.mode(x) <- "double"
y <- factor(c(rep("minor", 5), rep("major", 15)))
sby_balance_matrix(x, y, sby_strategy = "none")
}
\seealso{\code{\link{sby_adasyn_matrix}}, \code{\link{sby_nearmiss_index}}}
