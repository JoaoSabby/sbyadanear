# Gera vetor numerico com distribuicoes aleatorias controladas
#
# Funcao auxiliar para produzir vetores com distribuicoes
# estatisticas variadas, simulando ruido de fundo em features
# sinteticas. Inclui retornos antecipados para execucao otimizada.
#
# @param num_rows integer. Numero de linhas a gerar.
# @param data_type character. Tipo de dado a simular:
#   'monetary', 'ratio', 'integer' ou 'binary'.
# @param negative logical. Reservado para uso futuro.
#
# @return Vetor numerico, inteiro ou logico correspondente ao tipo solicitado.
# @keywords internal
sby_syn_build_vector <- function(
    num_rows = 4000000L,
    data_type = c("monetary", "ratio", "integer", "binary"),
    negative = FALSE
){
  data_type <- match.arg(data_type)

  if (data_type == "monetary") {
    dist_options <- c("uniform", "normal", "lognormal", "gamma", "weibull", "exponential", "chisq")
    chosen_dist <- sample(dist_options, 1L)

    result <- switch(chosen_dist,
      "uniform"     = runif(num_rows, 1000L, 40000L),
      "normal"      = abs(rnorm(num_rows, mean = 20000L, sd = 15000L)),
      "lognormal"   = rlnorm(num_rows, meanlog = 9.5, sdlog = 1),
      "gamma"       = rgamma(num_rows, shape = 2L, scale = 10000L),
      "weibull"     = rweibull(num_rows, shape = 1.5, scale = 20000L),
      "exponential" = rexp(num_rows, rate = 1 / 20000),
      rchisq(num_rows, df = 3L) * 5000L
    )

    result <- pmin(pmax(result, 1000L), 40000L)
    sparsity_mask <- runif(num_rows) > 0.2
    result[sparsity_mask] <- 0
    return(as.double(result))
  }

  if (data_type == "ratio") {
    dist_options <- c("uniform", "beta", "exponential", "folded_normal")
    chosen_dist <- sample(dist_options, 1L)

    result <- switch(chosen_dist,
      "uniform"     = runif(num_rows, 0, 1L),
      "beta"        = rbeta(num_rows, shape1 = 2L, shape2 = 5L),
      "exponential" = rexp(num_rows, rate = 5L),
      abs(rnorm(num_rows, mean = 0.2, sd = 0.3))
    )

    result <- pmin(result, 1L)
    sparsity_mask <- runif(num_rows) > 0.3
    result[sparsity_mask] <- 0
    return(as.double(result))
  }

  if (data_type == "integer") {
    dist_options <- c("poisson", "geometric", "negative_binomial", "binomial")
    chosen_dist <- sample(dist_options, 1L)

    result <- switch(chosen_dist,
      "poisson"            = rpois(num_rows, lambda = 2L),
      "geometric"          = rgeom(num_rows, prob = 0.3),
      "negative_binomial"  = rnbinom(num_rows, size = 1L, prob = 0.2),
      rbinom(num_rows, size = 10L, prob = 0.1)
    )

    sparsity_mask <- runif(num_rows) > 0.4
    result[sparsity_mask] <- 0L
    return(as.integer(result))
  }

  if (data_type == "binary") {
    prob_binomial <- runif(1L, 0.01, 0.15)
    result_raw <- runif(num_rows) < prob_binomial
    return(as.logical(result_raw))
  }

  stop("Tipo de dado invalido. Use 'monetary', 'ratio', 'integer' ou 'binary'.")
}
