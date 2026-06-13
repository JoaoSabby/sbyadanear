\name{sby_adasyn_hpc}
\alias{sby_adasyn_hpc}
\title{Atalho HPC do oversampling ADASYN}
\usage{
sby_adasyn_hpc(.data, formula, sby_k_neighbor_adanear = 3,
  sby_over_ratio = 0.2, sby_config_max_threads = -1)
}
\arguments{
\item{.data}{Data frame, tibble ou matriz com a coluna de desfecho e os preditores numericos.}
\item{formula}{Formula no formato alvo ~ preditores.}
\item{sby_k_neighbor_adanear}{Numero inteiro positivo de vizinhos da etapa ADASYN. O padrao e 3.}
\item{sby_over_ratio}{Fator relativo de expansao da classe minoritaria. O padrao e 0.2.}
\item{sby_config_max_threads}{Numero inteiro de threads do motor HPC. O padrao -1 detecta os nucleos fisicos disponiveis.}
}
\description{
Atalho de alto desempenho do oversampling ADASYN. Gera as linhas sinteticas
estritamente no espaco padronizado por z-score e monta o tibble final por zero-copy
diretamente em C++.
}
\details{
A funcao substitui internamente a rota denominada "native" como caminho rapido do
ADASYN. Quando o motor HPC consolidado nao esta disponivel, recorre de forma
transparente a \code{sby_adasyn()} com \code{sby_knn_engine = "native"}. Apenas \code{MKL_NUM_THREADS}, \code{OMP_NUM_THREADS} e \code{MKL_NUM_STRIPES} sao controladas
temporariamente e restauradas por um bloco \code{on.exit()} inflexivel.
}
\value{Tibble balanceado com classe \code{c("tbl_df", "tbl", "data.frame")}.}
\seealso{\code{\link{sby_adasyn}}, \code{\link{sby_adanear_hpc}}, \code{\link{sby_nearmiss_hpc}}}
