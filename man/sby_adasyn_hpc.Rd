\name{sby_adasyn_hpc}
\alias{sby_adasyn_hpc}
\title{Atalho HPC do oversampling ADASYN}
\usage{
sby_adasyn_hpc(.data, formula, sby_k_adasyn = 3,
  sby_over_ratio = 0.2, sby_config_max_threads = -1,
  sby_seed = sample.int(10L^5L, 1L))
}
\arguments{
\item{.data}{Data frame, tibble ou matriz com a coluna de desfecho e os preditores numericos.}
\item{formula}{Formula no formato alvo ~ preditores.}
\item{sby_k_adasyn}{Numero inteiro positivo de vizinhos da etapa ADASYN. O padrao e 3.}
\item{sby_over_ratio}{Fator relativo de expansao da classe minoritaria. O padrao e 0.2.}
\item{sby_config_max_threads}{Numero inteiro de threads do motor HPC. O padrao -1 detecta os nucleos fisicos disponiveis.}
\item{sby_seed}{Semente inteira para o gerador de numeros pseudoaleatorios. A semente e aplicada em escopo local e o estado RNG global do chamador e restaurado ao final. O padrao e \code{sample.int(10L^5L, 1L)}.}
}
\description{
Atalho de alto desempenho do oversampling ADASYN. Gera as linhas sinteticas no espaco padronizado usando estatisticas via MKL VSL, distancias por \code{cblas_sgemm} e pesos de interpolacao gerados por \code{Rcpp::runif()} sob controle da semente local.
}
\details{
A funcao substitui internamente a rota denominada "native" como caminho rapido do ADASYN. O estado de RNG do chamador e preservado: a funcao aplica \code{sby_seed} apenas durante a chamada nativa e restaura \code{RNGkind()} e \code{.Random.seed} ao final. O runtime MKL/OpenMP deve ser configurado externamente.
}
\value{Tibble balanceado com classe \code{c("tbl_df", "tbl", "data.frame")}.}
\seealso{\code{\link{sby_adasyn}}, \code{\link{sby_adanear_hpc}}, \code{\link{sby_nearmiss_hpc}}}
