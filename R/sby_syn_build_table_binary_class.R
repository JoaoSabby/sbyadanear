# Gera tabela sintetica densa para classificacao binaria
#
# Orquestra a geracao de um dataset completo e escalavel para
# testes de balanceamento. Retorna um tibble denso com ID,
# TARGET (factor binario) e features numericas, pronto para
# uso direto com sby_adanear(TARGET ~ ., data).
#
# @param num_rows integer. Numero de observacoes a gerar.
# @param target_prev numeric. Prevalencia da classe rara (0, 1).
#
# @return tibble com ID, TARGET e 250 features densas.
# @keywords internal
sby_syn_build_table_binary_class <- function(
    num_rows = 4000000L,
    target_prev = 0.0013
){
  message("Iniciando geracao de dados. Alocando memoria para ", num_rows, " linhas...")

  monetary_vars <- c(
    "TOTAL_BANK_CREDIT", "BANK_CREDIT_CARDS", "BANK_OVERDRAFT", "BANK_LOANS",
    "BANK_FINANCING", "BANK_MORTGAGES", "BANK_CREDIT_LIMIT_USED", "BANK_OTHER_CREDIT",
    "NAT_SYS_CREDIT_CARDS", "NAT_SYS_OVERDRAFT", "NAT_SYS_LOANS", "NAT_SYS_FINANCING",
    "NAT_SYS_MORTGAGES", "NAT_SYS_CREDIT_LIMIT_USED", "NAT_SYS_OTHER_CREDIT",
    "PAYCHECK_VAL", "TOTAL_INCOME", "ADJ_BANK_CREDIT_LIMIT", "BANK_CARD_LIMIT",
    "BANK_TOTAL_LIMIT", "MAX_BANK_LIMIT", "TOTAL_UNPAID_CREDIT", "TOTAL_NAT_SYS_CREDIT",
    "BANK_OVERDRAFT_LIMIT", "TOTAL_INVEST_VALUE", "TOTAL_CREDIT_VALUE",
    "CHECK_ACC_BALANCE", "SAVINGS_ACC_BALANCE", "MONTHLY_TRANSACTED_VALUE",
    "EFFECTIVE_MONTHLY_TRANCT", "TOTAL_CONTRIBUTION_MARGIN", "NET_SALARY_PORT",
    "INFLOW_SALARY_PORT", "OUTFLOW_SALARY_PORT", "NET_CREDIT_PORT", "INFLOW_CREDIT_PORT",
    "OUTFLOW_CREDIT_PORT"
  )

  ratio_vars <- c(
    "CREDIT_RATIO", "WALLET_SHARE_CC", "WALLET_SHARE_OD", "WALLET_SHARE_LOANS",
    "WALLET_SHARE_FINANCING", "WALLET_SHARE_TOTAL", "WALLET_SHARE_MORTGAGE",
    "WALLET_SHARE_CREDIT_LIMIT", "WALLET_SHARE_OTHER_CREDIT", "REL_PURCHASE_FREQ",
    "CONT_CREDIT_SCORE"
  )

  integer_vars <- c(
    "LIFE_CYCLE_IBGE", "COHORT", "TERMINATED_INVESTMENTS", "RECENT_INVESTMENTS",
    "TERMINATED_CREDITS", "RECENT_CREDIT_PURCHASE", "PURCHASE_FREQ", "MONTHLY_TRANSACTIONS",
    "PURCHASE_RECENCY", "DISTINCT_PRODUCTS", "DIRECT_DEBITS_COUNT", "MOBILE_INTERACTIONS",
    "MOBILE_TRANSACTIONS"
  )

  binary_vars <- c(
    "NATIONAL_CREDIT", "MARRIED", "DIVORCED", "SEPARATED", "WIDOWED", "CIVIL_UNION",
    "SINGLE", "PARTLY_DIVORCED", "MARITAL_UNKNOWN", "NO_DIGITAL_ACTIVITY", "DIGITAL_QUERIES",
    "DIGITAL_PAYMENTS", "DIGITAL_INVEST_CREDIT", "HS_INCOMPLETE", "HIGH_SCHOOL",
    "PROF_TECH_INCOMPLETE", "PROF_TECH", "GRAD_INCOMPLETE", "GRADUATE", "POST_GRADUATE",
    "MASTERS", "PHD", "EDU_OTHER", "ILLITERATE", "FUNC_LITERATE", "ELEM_SCHOOL_INC",
    "ELEM_SCHOOL", "EDU_UNKNOWN", "CUST_SEG_B_AA", "CUST_SEG_B_AI", "CUST_SEG_B_AN",
    "CUST_SEG_B_AT", "CUST_SEG_B_ATS", "CUST_SEG_B_INACTIVE", "CUST_SEG_B_UNDEF",
    "CUST_SEG_B_PA", "CUST_SEG_B_PI", "CUST_SEG_B_PN", "CUST_SEG_B_PT", "CUST_SEG_B_PTS",
    "CUST_SEG_B_VA", "CUST_SEG_B_VI", "CUST_SEG_B_VN", "CUST_SEG_B_VT", "CUST_SEG_B_VTS",
    "CUST_SEG_A_1", "CUST_SEG_A_2", "CUST_SEG_A_3", "CUST_SEG_C_HIGH_VAL",
    "CUST_SEG_C_MED_VAL", "CUST_SEG_C_UNDEF", "CUST_SEG_C_LOW_VAL", "FEMALE", "MALE",
    "GENDER_UNKNOWN", "DEDICATED_ACC_MGR", "CREDIT_SCORE_1", "CREDIT_SCORE_10",
    "CREDIT_SCORE_11", "CREDIT_SCORE_12", "CREDIT_SCORE_2", "CREDIT_SCORE_3",
    "CREDIT_SCORE_4", "CREDIT_SCORE_5", "CREDIT_SCORE_6", "CREDIT_SCORE_7",
    "CREDIT_SCORE_8", "CREDIT_SCORE_9", "CREDIT_SCORE_99", "CREDIT_SCORE_UNDEF",
    "CIVIL_SERVANT", "PROF_PROFILE_1001", "PROF_PROFILE_1002", "PROF_PROFILE_1003",
    "PROF_PROFILE_1004", "PROF_PROFILE_1005", "PROF_PROFILE_1006", "PROF_PROFILE_1007",
    "PROF_PROFILE_1008", "PROF_PROFILE_1009", "PROF_PROFILE_1010", "PROF_PROFILE_1011",
    "PROF_PROFILE_1015", "PROF_PROFILE_1017", "PROF_PROFILE_1018", "PROF_PROFILE_1019",
    "PROF_PROFILE_1020", "PROF_PROFILE_1021", "PROF_PROFILE_1022", "PROF_PROFILE_1023",
    "PROF_PROFILE_1025", "PROF_PROFILE_1026", "PROF_PROFILE_1029", "PROF_PROFILE_1030",
    "PROF_PROFILE_1031", "PROF_PROFILE_1032", "PROF_PROFILE_1034", "PROF_PROFILE_2001",
    "PROF_PROFILE_2002", "PROF_PROFILE_2003", "PROF_PROFILE_2004", "PROF_PROFILE_2005",
    "PROF_PROFILE_2006", "PROF_PROFILE_2007", "PROF_PROFILE_2008", "PROF_PROFILE_2009",
    "PROF_PROFILE_2010", "PROF_PROFILE_2011", "PROF_PROFILE_2015", "PROF_PROFILE_2017",
    "PROF_PROFILE_2018", "PROF_PROFILE_2019", "PROF_PROFILE_2020", "PROF_PROFILE_2021",
    "PROF_PROFILE_2022", "PROF_PROFILE_2023", "PROF_PROFILE_2025", "PROF_PROFILE_2026",
    "PROF_PROFILE_2028", "PROF_PROFILE_2029", "PROF_PROFILE_2030", "PROF_PROFILE_2031",
    "PROF_PROFILE_2032", "PROF_PROFILE_2034", "PROF_PROFILE_NONE", "HAS_CREDIT",
    "SALARY_WITH_BANK", "OVERDUE_CREDIT", "INVESTMENT_INCIDENCE", "CREDIT_PURCHASE_INCIDENCE",
    "AUTO_DEBIT_USAGE", "PURCHASE_INCIDENCE", "DEBT_RENEGOTIATION", "HAS_CREDIT_CARD",
    "USED_CREDIT_CARD", "BRANCH_CATEGORY_A", "BRANCH_CATEGORY_B", "BRANCH_CATEGORY_C",
    "BRANCH_CATEGORY_D", "BRANCH_CATEGORY_E", "BRANCH_CATEGORY_F", "BRANCH_CATEGORY_NONE",
    "BRANCH_REGION_0", "BRANCH_REGION_1", "BRANCH_REGION_14", "BRANCH_REGION_15",
    "BRANCH_REGION_2", "BRANCH_REGION_3", "BRANCH_REGION_34", "BRANCH_REGION_37",
    "BRANCH_REGION_38", "BRANCH_REGION_4", "BRANCH_REGION_45", "BRANCH_REGION_46",
    "BRANCH_REGION_5", "BRANCH_REGION_6", "BRANCH_REGION_7", "BRANCH_REGION_8",
    "BRANCH_REGION_9", "MOBILE_USE_INCIDENCE", "SALARY_PORT_OUTFLOW", "SALARY_PORT_INFLOW",
    "NO_SALARY_PORT", "PREV_CREDIT_PORT", "CREDIT_PORT_INFLOW", "ANY_CREDIT_PORT",
    "NO_CREDIT_PORT", "MONTH_CREDIT_PORT",
    "MARTIAL_ART_KUNGFU", "MARTIAL_ART_TAEKWONDO", "MARTIAL_ART_SUMO",
    "OWNER_CAT_BENGAL", "OWNER_CAT_BIRMAN", "LISTENER_LED_ZEPPELIN",
    "LISTENER_ROBERT_PLANT", "CINEPHILE_INGMAR_BERGMAN", "CINEPHILE_JOHN_CARPENTER",
    "CINEPHILE_ELISEO_SUBIELA", "CONSUMPTION_COENZYME_Q10", "CONSUMPTION_OMEGA_3",
    "CONSUMPTION_MAGNESIUM", "INVESTOR_REIT", "PURCHASE_DUSTER_RADIATOR", "TRIP_XANGRILA"
  )

  all_vars <- c(monetary_vars, ratio_vars, integer_vars, binary_vars)

  # Garante exatamente 250 features preenchendo com binarias extras
  if (length(all_vars) < 250L) {
    missing_count <- 250L - length(all_vars)
    padding_vars <- paste0("RANDOM_FEATURE_", sprintf("%03d", seq_len(missing_count)))
    binary_vars <- c(binary_vars, padding_vars)
    all_vars <- c(all_vars, padding_vars)
  } else if (length(all_vars) > 250L) {
    all_vars <- all_vars[1:250L]
  }

  set.seed(42L)
  shuffled_vars <- sample(all_vars)

  # Seleciona 20 features que terao sinal real com o TARGET
  signal_vars <- sample(shuffled_vars, 20L)

  message("Gerando vetores densos com distribuicoes aleatorias controladas...")

  # Gera cada coluna como vetor denso usando base R (sem purrr/sparsevctrs)
  col_data <- vector("list", length(shuffled_vars))
  names(col_data) <- shuffled_vars

  for (idx in seq_along(shuffled_vars)) {
    var_name <- shuffled_vars[idx]
    is_signal <- var_name %in% signal_vars

    if (var_name %in% monetary_vars) {
      if (is_signal) {
        result_raw <- runif(num_rows, 1000L, 40000L)
        sparsity_mask <- runif(num_rows) > 0.2
        result_raw[sparsity_mask] <- 0
        col_data[[idx]] <- as.double(result_raw)
      } else {
        col_data[[idx]] <- sby_syn_build_vector(num_rows, "monetary")
      }
      next
    }

    if (var_name %in% ratio_vars) {
      if (is_signal) {
        result_raw <- runif(num_rows, 0, 1L)
        sparsity_mask <- runif(num_rows) > 0.3
        result_raw[sparsity_mask] <- 0
        col_data[[idx]] <- as.double(result_raw)
      } else {
        col_data[[idx]] <- sby_syn_build_vector(num_rows, "ratio")
      }
      next
    }

    if (var_name %in% integer_vars) {
      if (is_signal) {
        result_raw <- rpois(num_rows, lambda = 2L)
        sparsity_mask <- runif(num_rows) > 0.4
        result_raw[sparsity_mask] <- 0L
        col_data[[idx]] <- as.integer(result_raw)
      } else {
        col_data[[idx]] <- sby_syn_build_vector(num_rows, "integer")
      }
      next
    }

    # binary_vars (incluindo padding)
    if (is_signal) {
      col_data[[idx]] <- as.logical(runif(num_rows) < 0.05)
    } else {
      col_data[[idx]] <- sby_syn_build_vector(num_rows, "binary")
    }
  }

  message("Calculando relacoes nao-lineares para o TARGET...")

  latent_score <- numeric(num_rows)

  for (i in seq(1L, 19L, by = 2L)) {
    var1 <- as.numeric(col_data[[ signal_vars[i] ]])
    var2 <- as.numeric(col_data[[ signal_vars[i + 1L] ]])

    # Usa collapse::fmean e collapse::fsd para normalizacao rapida
    v1_mean <- collapse::fmean(var1)
    v1_sd   <- collapse::fsd(var1)
    v2_mean <- collapse::fmean(var2)
    v2_sd   <- collapse::fsd(var2)

    v1_norm <- (var1 - v1_mean) / (v1_sd + 1e-6)
    v2_norm <- (var2 - v2_mean) / (v2_sd + 1e-6)

    latent_score <- latent_score + v1_norm * v2_norm
  }

  gaussian_noise <- rnorm(num_rows, 0, collapse::fsd(latent_score) * 0.4)
  latent_score <- latent_score + gaussian_noise

  cutoff_value <- stats::quantile(latent_score, probs = 1 - target_prev)

  # Usa kit::iif para geracao condicional vetorizada sem overhead de ifelse
  target_vector_raw <- kit::iif(latent_score >= cutoff_value, "Yes", "No")
  target_vector <- factor(target_vector_raw, levels = c("Yes", "No"))

  rm(latent_score, gaussian_noise, cutoff_value, target_vector_raw)

  message("Consolidando o tibble denso final...")

  # Converte todas as colunas para double denso (sem sparsevctrs)
  for (idx in seq_along(col_data)) {
    col_data[[idx]] <- as.double(col_data[[idx]])
  }

  features_tibble <- tibble::as_tibble(col_data)

  final_data <- tibble::tibble(
    ID = seq_len(num_rows),
    TARGET = target_vector
  )
  final_data <- cbind(final_data, features_tibble)
  final_data <- tibble::as_tibble(final_data)

  message("Geracao concluida. ", collapse::fnrow(final_data), " linhas x ",
          collapse::fncol(final_data), " colunas.")

  return(final_data)
}
