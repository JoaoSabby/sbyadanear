\name{sby_nearmiss_index}
\alias{sby_nearmiss_index}
\title{Selecionar índices NearMiss-1}
\usage{
sby_nearmiss_index(sby_x_matrix, sby_y_vector, sby_under_ratio = 0.5,
  sby_knn_under_k = 5L, sby_seed = sample.int(10L^5L, 1L),
  sby_scaling_info = NULL, sby_input_already_scaled = FALSE,
  sby_fixed_minority_label = NULL, sby_fixed_majority_label = NULL,
  sby_knn_algorithm = c("auto", "kd_tree", "cover_tree", "brute"),
  sby_knn_engine = c("auto", "FNN", "RcppHNSW"),
  sby_knn_distance_metric = c("euclidean", "ip", "cosine"),
  sby_knn_workers = 1L, sby_knn_hnsw_m = 16L,
  sby_knn_hnsw_ef = 200L, sby_audit = FALSE,
  sby_audit_level = c("none", "light", "full"),
  sby_return_scaling_info = FALSE, sby_return_reduced_scaled = FALSE)
}
\arguments{
\item{sby_x_matrix}{Matrix double densa usada para calcular os índices.}
\item{sby_y_vector}{Fator binário.}
\item{sby_under_ratio, sby_knn_under_k, sby_seed}{Controles NearMiss.}
\item{sby_scaling_info, sby_input_already_scaled}{Reuso de escala.}
\item{sby_fixed_minority_label, sby_fixed_majority_label}{Papéis fixos de classe.}
\item{sby_knn_algorithm, sby_knn_engine, sby_knn_distance_metric, sby_knn_workers, sby_knn_hnsw_m, sby_knn_hnsw_ef}{Configuração KNN.}
\item{sby_audit, sby_audit_level}{Auditoria.}
\item{sby_return_scaling_info}{Retorna parâmetros de escala sem exigir auditoria full.}
\item{sby_return_reduced_scaled}{Retorna matriz reduzida em z-score sem expor a matriz escalada completa.}
}
\description{Retorna índices de retenção NearMiss-1 para que pipelines possam aplicar a seleção em matrizes externas.}
\details{Útil para preservar matrizes esparsas fora do pacote quando a estratégia for apenas NearMiss. A seleção nativa retém menores distâncias médias; empates favorecem o menor índice original.}
\value{Lista com \code{sby_retained_index}, distribuições e diagnósticos; opcionalmente escala e matriz reduzida escalada.}
\examples{
x <- matrix(rnorm(40), ncol = 2); storage.mode(x) <- "double"
y <- factor(c(rep("minor", 5), rep("major", 15)))
sby_nearmiss_index(x, y, sby_seed = 1, sby_knn_engine = "FNN")
}
\seealso{\code{\link{sby_nearmiss_matrix}}, \code{\link{sby_balance_matrix}}}
