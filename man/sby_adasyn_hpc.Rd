\name{sby_adasyn_hpc}
\alias{sby_adasyn_hpc}
\title{Atalho HPC do oversampling ADASYN}
\usage{
sby_adasyn_hpc(.data, formula, sby_adasyn_k = 3,
  sby_adasyn_ratio = 0.2, sby_config_max_threads = -1,
  sby_seed = sample.int(10L^5L, 1L))
}
\arguments{
\item{.data}{Data frame, tibble ou matriz com a coluna de desfecho e os preditores numericos.}
\item{formula}{Formula no formato alvo ~ preditores.}
\item{sby_adasyn_k}{Numero inteiro positivo de vizinhos da etapa ADASYN. O padrao e 3.}
\item{sby_adasyn_ratio}{Acréscimo relativo sobre a quantidade original da classe rara. O cálculo sempre soma \code{1} ao valor informado: \code{n_rara_final = floor(n_rara_original * (1 + sby_adasyn_ratio))}. Portanto, \code{0.4} usa o multiplicador final \code{1.4}, nunca \code{0.4}, e preserva todos os raros originais.}
\item{sby_config_max_threads}{Numero inteiro de threads do motor HPC. O padrao -1 detecta os nucleos fisicos disponiveis.}
\item{sby_seed}{Semente inteira para o gerador de numeros pseudoaleatorios. A semente e aplicada em escopo local e o estado RNG global do chamador e restaurado ao final. O padrao e \code{sample.int(10L^5L, 1L)}.}
}
\description{
Atalho de alto desempenho do oversampling ADASYN. Gera as linhas sinteticas no espaco padronizado usando estatisticas populacionais por laco SIMD paralelo, distancias exatas por \code{sgemm} (oneMKL quando ligado, BLAS do R caso contrario) e pesos de interpolacao gerados por \code{Rcpp::runif()} sob controle da semente local. A alocacao das sinteticas e ponderada pela densidade majoritaria local de cada ponto raro, como no ADASYN classico.
}
\details{
A funcao substitui internamente a rota denominada "native" como caminho rapido do ADASYN. O estado de RNG do chamador e preservado: a funcao aplica \code{sby_seed} apenas durante a chamada nativa e restaura \code{RNGkind()} e \code{.Random.seed} ao final. O runtime MKL/OpenMP deve ser configurado externamente.
}
\value{Tibble balanceado com classe \code{c("tbl_df", "tbl", "data.frame")}.}
\seealso{\code{\link{sby_adasyn}}, \code{\link{sby_adanear_hpc}}, \code{\link{sby_nearmiss_hpc}}}
