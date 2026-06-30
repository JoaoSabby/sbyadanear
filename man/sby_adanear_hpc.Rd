\name{sby_adanear_hpc}
\alias{sby_adanear_hpc}
\title{Atalho HPC do pipeline combinado ADASYN e NearMiss-1}
\usage{
sby_adanear_hpc(.data, formula, sby_k_adasyn = 3,
  sby_k_nearmiss = 7, sby_config_max_threads = -1,
  sby_seed = sample.int(10L^5L, 1L), sby_over_ratio = 0.2,
  sby_under_ratio = 1)
}
\arguments{
\item{.data}{Data frame, tibble ou matriz com a coluna de desfecho e os preditores numericos.}
\item{formula}{Formula no formato alvo ~ preditores.}
\item{sby_k_adasyn}{Numero inteiro positivo de vizinhos da etapa ADASYN. O padrao e 3.}
\item{sby_k_nearmiss}{Numero inteiro positivo de vizinhos da etapa NearMiss-1. O padrao e 7.}
\item{sby_config_max_threads}{Numero inteiro de threads do motor HPC. O padrao -1 detecta os nucleos fisicos disponiveis.}
\item{sby_seed}{Valor numerico inteiro utilizado para inicializar o gerador de numeros pseudoaleatorios da etapa ADASYN no motor HPC. A semente e aplicada em escopo local e o estado RNG global do chamador e restaurado ao final. O padrao e \code{sample.int(10L^5L, 1L)}.}
\item{sby_over_ratio}{Razao nao negativa de aumento da classe rara. Valores positivos executam ADASYN; zero desativa ADASYN nesta rotina hibrida. O padrao e 0.2.}
\item{sby_under_ratio}{Razao nao negativa de retencao da classe majoritaria em relacao ao tamanho final da classe rara. Valores positivos executam NearMiss-1; zero desativa NearMiss-1 nesta rotina hibrida. O padrao e 1.}
}
\description{
Atalho de alto desempenho do pipeline combinado de balanceamento binario. Consolida ADASYN e NearMiss-1 em uma unica passada no espaco padronizado, usando estatisticas via MKL VSL, distancias por \code{cblas_sgemm} e pesos de interpolacao gerados por \code{Rcpp::runif()} sob controle da semente local.
}
\details{
A funcao substitui internamente a rota denominada "native" como caminho rapido. O estado de RNG do chamador e preservado: a funcao aplica \code{sby_seed} apenas durante a chamada nativa e restaura \code{RNGkind()} e \code{.Random.seed} ao final. O runtime MKL/OpenMP deve ser configurado externamente.

Regras formais das razoes de reamostragem:

Sejam \eqn{r_o = sby\_over\_ratio}, \eqn{r_u = sby\_under\_ratio}, \eqn{n_{min}^{(0)}} o numero original de registros da classe rara, \eqn{n_{maj}^{(0)}} o numero original de registros da classe majoritaria, \eqn{n_{syn}} o numero de registros sinteticos gerados pelo ADASYN, \eqn{n_{min}^{(1)}} o numero de registros da classe rara disponiveis para o NearMiss-1, \eqn{n_{maj}^{disp}} o numero de registros majoritarios disponiveis antes do NearMiss-1 e \eqn{n_{maj}^{ret}} o numero-alvo de registros majoritarios retidos. O dominio das razoes nesta rotina hibrida e:
\deqn{r_o, r_u \in [0, \infty).}
Valores negativos sao invalidos.

A etapa ADASYN e executada se, e somente se:
\deqn{r_o > 0.}
Formalmente:
\deqn{\operatorname{ADASYN}(r_o) =
\begin{cases}
\text{nao executado}, & \text{se } r_o = 0 \\
\text{executado}, & \text{se } r_o > 0
\end{cases}}
Quando \eqn{r_o > 0}, o numero de registros sinteticos e aproximadamente:
\deqn{n_{syn} = \left\lfloor n_{min}^{(0)} r_o \right\rfloor}
e a classe rara apos ADASYN tem aproximadamente:
\deqn{n_{min}^{(1)} = n_{min}^{(0)} + n_{syn} \approx n_{min}^{(0)}(1 + r_o).}
Em bases pequenas, valores positivos podem gerar ao menos uma linha sintetica, conforme a politica interna de arredondamento.

A etapa NearMiss-1 e executada se, e somente se:
\deqn{r_u > 0.}
Formalmente:
\deqn{\operatorname{NearMiss}(r_u) =
\begin{cases}
\text{nao executado}, & \text{se } r_u = 0 \\
\text{executado}, & \text{se } r_u > 0
\end{cases}}
Quando \eqn{r_u > 0}, a quantidade-alvo de registros majoritarios retidos e:
\deqn{n_{maj}^{ret} = \min\left(n_{maj}^{disp}, \left\lfloor n_{min}^{(1)} r_u \right\rfloor\right).}
Assim, \code{sby_under_ratio = 0.5} retem ate metade do tamanho final da classe rara, \code{sby_under_ratio = 1} retem ate a mesma quantidade da classe rara e \code{sby_under_ratio = 2} retem ate duas vezes o tamanho da classe rara, limitado a maioria disponivel.

A logica do pipeline hibrido e:
\itemize{
  \item \eqn{r_o = 0 \land r_u = 0}: ADASYN e NearMiss-1 nao executam; os dados originais sao retornados preservados.
  \item \eqn{r_o > 0 \land r_u = 0}: apenas ADASYN executa.
  \item \eqn{r_o = 0 \land r_u > 0}: apenas NearMiss-1 executa.
  \item \eqn{r_o > 0 \land r_u > 0}: ADASYN executa seguido de NearMiss-1.
}
}
\value{Tibble balanceado com classe \code{c("tbl_df", "tbl", "data.frame")}.}
\seealso{\code{\link{sby_adanear}}, \code{\link{sby_adasyn_hpc}}, \code{\link{sby_nearmiss_hpc}}}
