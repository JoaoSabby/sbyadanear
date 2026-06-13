\name{sby_nearmiss_hpc}
\alias{sby_nearmiss_hpc}
\title{Atalho HPC do undersampling NearMiss-1}
\usage{
sby_nearmiss_hpc(.data, formula, sby_k_neighbor_nearmiss = 7,
  sby_under_ratio = 0.5, sby_config_max_threads = -1)
}
\arguments{
\item{.data}{Data frame, tibble ou matriz com a coluna de desfecho e os preditores numericos.}
\item{formula}{Formula no formato alvo ~ preditores.}
\item{sby_k_neighbor_nearmiss}{Numero inteiro positivo de vizinhos da etapa NearMiss-1. O padrao e 7.}
\item{sby_under_ratio}{Razao minima desejada entre minoria e maioria apos o NearMiss-1. O padrao e 0.5.}
\item{sby_config_max_threads}{Numero inteiro de threads do motor HPC. O padrao -1 detecta os nucleos fisicos disponiveis.}
}
\description{
Atalho de alto desempenho do undersampling NearMiss-1. Rankeia e retem as linhas
majoritarias estritamente no espaco padronizado por z-score e monta o tibble final
por zero-copy diretamente em C++.
}
\details{
A funcao substitui internamente a rota denominada "native" como caminho rapido do
NearMiss-1. Quando o motor HPC consolidado nao esta disponivel, recorre de forma
transparente a \code{sby_nearmiss()} com \code{sby_knn_engine = "native"}. Apenas \code{MKL_NUM_THREADS} e \code{OMP_NUM_THREADS} sao controladas
temporariamente e restauradas por um bloco \code{on.exit()} inflexivel.
}
\value{Tibble balanceado com classe \code{c("tbl_df", "tbl", "data.frame")}.}
\seealso{\code{\link{sby_nearmiss}}, \code{\link{sby_adanear_hpc}}, \code{\link{sby_adasyn_hpc}}}
