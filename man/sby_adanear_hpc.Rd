\name{sby_adanear_hpc}
\alias{sby_adanear_hpc}
\title{Atalho HPC do pipeline combinado ADASYN e NearMiss-1}
\usage{
sby_adanear_hpc(.data, formula, sby_k_neighbor_adanear = 3,
  sby_k_neighbor_nearmiss = 7, sby_config_max_threads = -1)
}
\arguments{
\item{.data}{Data frame, tibble ou matriz com a coluna de desfecho e os preditores numericos.}
\item{formula}{Formula no formato alvo ~ preditores.}
\item{sby_k_neighbor_adanear}{Numero inteiro positivo de vizinhos da etapa ADASYN. O padrao e 3.}
\item{sby_k_neighbor_nearmiss}{Numero inteiro positivo de vizinhos da etapa NearMiss-1. O padrao e 7.}
\item{sby_config_max_threads}{Numero inteiro de threads do motor HPC. O padrao -1 detecta os nucleos fisicos disponiveis.}
}
\description{
Atalho de alto desempenho do pipeline combinado de balanceamento binario. Consolida
ADASYN e NearMiss-1 em uma unica passada no espaco padronizado, eliminando a dupla
normalizacao e montando o tibble final por zero-copy diretamente em C++.
}
\details{
A funcao substitui internamente a rota denominada "native" como caminho rapido.
Quando o motor HPC consolidado nao esta disponivel, recorre de forma transparente
a \code{sby_adanear()} com \code{sby_knn_engine = "native"}. As variaveis NUMA e MKL
sao isoladas e restauradas por um bloco \code{on.exit()} inflexivel. As estatisticas
iniciais usam a Vector Statistics Library, a matriz de distancias usa
\code{cblas_dgemm} e a reversao final do z-score aproveita as unidades FMA do hardware.
}
\value{Tibble balanceado com classe \code{c("tbl_df", "tbl", "data.frame")}.}
\seealso{\code{\link{sby_adanear}}, \code{\link{sby_adasyn_hpc}}, \code{\link{sby_nearmiss_hpc}}}
