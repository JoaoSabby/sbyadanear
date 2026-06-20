#' Executar expressao longa em processo filho interrompivel
#'
#' @details
#' A funcao isola chamadas nativas externas que nao consultam diretamente a flag
#' de interrupcao do R. No Unix, a expressao roda em fork e o processo pai fica
#' responsivo a Ctrl+C enquanto aguarda o resultado. Em plataformas sem fork, a
#' expressao e avaliada diretamente para preservar compatibilidade.
#'
#' @param sby_expr Expressao R a ser executada
#'
#' @return Resultado produzido por `sby_expr`
#'
#' @noRd
sby_run_interruptible_fork <- function(sby_expr){
  # Captura expressao e ambiente chamador antes de decidir o modo de execucao
  sby_quoted_expr <- substitute(sby_expr)
  sby_parent_env  <- parent.frame()

  # Usa execucao direta fora de Unix, onde fork nao esta disponivel de forma segura
  if(!identical(.Platform$OS.type, "unix")){
    return(eval(
      expr = sby_quoted_expr,
      envir = sby_parent_env
    ))
  }

  # Garante que uma interrupcao pendente seja tratada antes de iniciar o filho
  sby_adanear_check_user_interrupt()

  # Executa a chamada potencialmente bloqueante em processo filho
  sby_job <- parallel::mcparallel(
    expr = eval(
      expr = sby_quoted_expr,
      envir = sby_parent_env
    ),
    silent = TRUE
  )

  # Registra se o resultado ja foi coletado para evitar sinalizar processo encerrado
  sby_job_done <- FALSE

  # Finaliza o processo filho caso o usuario interrompa a espera no processo pai
  on.exit({
    if(!isTRUE(sby_job_done)){
      sby_remaining <- parallel::mccollect(
        jobs = sby_job,
        wait = FALSE
      )
      if(is.null(sby_remaining)){
        tools::pskill(
          pid = sby_job$pid
        )
        Sys.sleep(
          time = 0.05
        )
        tools::pskill(
          pid = sby_job$pid,
          signal = 9L
        )
        parallel::mccollect(
          jobs = sby_job,
          wait = FALSE
        )
      }
    }
  }, add = TRUE)

  # Aguarda cooperativamente, permitindo que Ctrl+C seja observado pelo processo pai
  repeat{
    sby_result <- parallel::mccollect(
      jobs = sby_job,
      wait = FALSE
    )

    if(!is.null(sby_result)){
      sby_job_done <- TRUE
      sby_value <- sby_result[[1L]]

      if(inherits(sby_value, "try-error")){
        stop(
          conditionMessage(attr(sby_value, "condition")),
          call. = FALSE
        )
      }

      return(sby_value)
    }

    sby_adanear_check_user_interrupt()
    Sys.sleep(
      time = 0.05
    )
  }
}
####
## Fim
#
