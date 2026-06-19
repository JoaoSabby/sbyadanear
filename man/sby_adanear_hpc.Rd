\name{sby_adanear_hpc}
\alias{sby_adanear_hpc}
\title{Atalho HPC do pipeline combinado ADASYN e NearMiss-1}
\usage{
sby_adanear_hpc(.data, formula, sby_k_neighbor_adanear = 3,
  sby_k_neighbor_nearmiss = 7, sby_config_max_threads = -1,
  sby_seed = sample.int(10L^5L, 1L), sby_over_ratio = 0.2,
  sby_under_ratio = 1)
}
\arguments{
\item{.data}{Data frame, tibble ou matriz com a coluna de desfecho e os preditores numericos.}
\item{formula}{Formula no formato alvo ~ preditores.}
\item{sby_k_neighbor_adanear}{Numero inteiro positivo de vizinhos da etapa ADASYN. O padrao e 3.}
\item{sby_k_neighbor_nearmiss}{Numero inteiro positivo de vizinhos da etapa NearMiss-1. O padrao e 7.}
\item{sby_config_max_threads}{Numero inteiro de threads do motor HPC. O padrao -1 detecta os nucleos fisicos disponiveis.}
\item{sby_seed}{Valor numerico inteiro utilizado para inicializar o gerador de numeros pseudoaleatorios da etapa ADASYN no motor HPC. O padrao e \code{sample.int(10L^5L, 1L)}, gerando uma semente aleatoria quando o usuario nao informa valor. Informe uma semente fixa para tornar reproduziveis as escolhas de vizinhos sinteticos e pesos de interpolacao.}
\item{sby_over_ratio}{Fator relativo de expansao da classe minoritaria pela etapa ADASYN. O padrao e 0.2.}
\item{sby_under_ratio}{Razao minima desejada entre minoria e maioria apos o NearMiss-1. O padrao e 1.}
}
\description{
Atalho de alto desempenho do pipeline combinado de balanceamento binario. Consolida
ADASYN e NearMiss-1 em uma unica passada no espaco padronizado, eliminando a dupla
normalizacao e montando o tibble final por zero-copy diretamente em C++.
}
\details{
A funcao substitui internamente a rota denominada "native" como caminho rapido.
Quando o motor HPC consolidado nao esta disponivel, recorre de forma transparente
a \code{sby_adanear()} com \code{sby_knn_engine = "native"}. Apenas \code{MKL_NUM_THREADS}, \code{OMP_NUM_THREADS} e \code{MKL_NUM_STRIPES}
sao controladas temporariamente e restauradas por um bloco \code{on.exit()} inflexivel. As estatisticas
iniciais usam a Vector Statistics Library, a matriz de distancias usa
\code{cblas_sgemm} e a reversao final do z-score aproveita as unidades FMA do hardware.
}
\value{Tibble balanceado com classe \code{c("tbl_df", "tbl", "data.frame")}.}
\seealso{\code{\link{sby_adanear}}, \code{\link{sby_adasyn_hpc}}, \code{\link{sby_nearmiss_hpc}}}
