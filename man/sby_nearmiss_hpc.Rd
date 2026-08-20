\name{sby_nearmiss_hpc}
\alias{sby_nearmiss_hpc}
\title{Atalho HPC do undersampling NearMiss-1}
\usage{
sby_nearmiss_hpc(.data, formula, sby_nearmiss_k = 7,
  sby_nearmiss_ratio = 1, sby_config_max_threads = -1,
  sby_seed = sample.int(10L^5L, 1L))
}
\arguments{
\item{.data}{Data frame, tibble ou matriz com a coluna de desfecho e os preditores numericos.}
\item{formula}{Formula no formato alvo ~ preditores.}
\item{sby_nearmiss_k}{Numero inteiro positivo de vizinhos da etapa NearMiss-1. O padrao e 7.}
\item{sby_nearmiss_ratio}{Razao nao negativa de retencao da classe majoritaria em relacao ao tamanho da classe rara. O alvo e \code{floor(n_minoria * sby_nearmiss_ratio)}, limitado a maioria disponivel. Zero desativa o undersampling e o padrao e 1.}
\item{sby_config_max_threads}{Numero inteiro de threads do motor HPC. O padrao -1 detecta os nucleos fisicos disponiveis. O limite vale para OpenMP e para a oneMKL da chamada corrente, quando ligada.}
\item{sby_seed}{Semente inteira validada para reproducibilidade da interface. O padrao e \code{sample.int(10L^5L, 1L)}.}
}
\description{
Atalho de alto desempenho do undersampling NearMiss-1. Rankeia e retem as linhas
majoritarias estritamente no espaco padronizado por z-score e monta o tibble final
por zero-copy diretamente em C++.
}
\details{
A funcao substitui internamente a rota denominada "native" como caminho rapido do
NearMiss-1. Quando o motor HPC consolidado nao esta disponivel, recorre de forma
transparente a \code{sby_nearmiss()} com \code{sby_knn_engine = "native"}. As variaveis de
ambiente do runtime MKL/OpenMP nao sao alteradas; o numero de threads vale
apenas para a chamada corrente, com os limites OpenMP e local da oneMKL salvos e
restaurados em torno do kernel nativo.
}
\value{Tibble balanceado com classe \code{c("tbl_df", "tbl", "data.frame")}.}
\seealso{\code{\link{sby_nearmiss}}, \code{\link{sby_adanear_hpc}}, \code{\link{sby_adasyn_hpc}}}
