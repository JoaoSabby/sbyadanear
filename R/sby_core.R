# Shared implementation for every public sampling interface. There is no
# alternative geometry, single-precision buffer, or serial native fallback.
sby_native <- function(name, ...) .Call(name, ..., PACKAGE = "sbyadanear")

sby_scalar <- function(x, name, integer = FALSE, lower = 0, upper = Inf) {
  if (!(is.integer(x) || is.double(x)) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x)))
    stop(name, " deve ser um número ", if (integer) "inteiro " else "",
         "entre ", lower, " e ", upper, ".", call. = FALSE)
  if (integer) as.integer(x) else as.double(x)
}
sby_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " deve ser TRUE ou FALSE.", call. = FALSE)
  x
}
sby_model <- function(x) {
  if (!(is.integer(x) || is.double(x)) || length(x) != 1L || is.na(x) ||
      !(x %in% 1:3) || x != floor(x))
    stop("nearmiss_model deve ser 1L (NearMiss-1), 2L (NearMiss-2) ou 3L (NearMiss-3).", call. = FALSE)
  as.integer(x)
}
sby_quota <- function(root="/sys/fs/cgroup", membership="/proc/self/cgroup") {
  read <- function(path) {
    if (!file.exists(path)) return(character())
    tryCatch(readLines(path, warn=FALSE),error=function(e) character())
  }
  # Check the current cgroup and its ancestors. The mount root alone may be
  # unlimited while a systemd service or a nested container is constrained.
  roots <- c(root,file.path(root,"cpu"),file.path(root,"cpu,cpuacct"))
  entries <- strsplit(read(membership),":",fixed=TRUE)
  paths <- list(character(),character(),character())
  for(v in entries) if(length(v)==3L) {
    if(v[2L]=="") paths[[1L]] <- c(paths[[1L]],v[3L])
    else if("cpu" %in% strsplit(v[2L],",",fixed=TRUE)[[1L]]) {
      paths[[2L]] <- c(paths[[2L]],v[3L])
      paths[[3L]] <- c(paths[[3L]],v[3L])
    }
  }
  quotas <- numeric()
  for(i in seq_along(roots)) {
    base <- normalizePath(roots[i],winslash="/",mustWork=FALSE)
    candidates <- base
    for(node in paths[[i]]) {
      current <- normalizePath(file.path(base,sub("^/","",node)),winslash="/",mustWork=FALSE)
      while(startsWith(current,paste0(base,"/"))) {
        candidates <- c(candidates,current)
        parent <- dirname(current)
        if(identical(parent,current)) break
        current <- parent
      }
    }
    for(path in unique(candidates)) {
      if(i==1L) {
        v <- strsplit(paste(read(file.path(path,"cpu.max")),collapse="")," +")[[1L]]
        if(length(v)!=2L || v[1L]=="max") next
        q <- suppressWarnings(as.double(v))
      } else q <- suppressWarnings(as.double(c(
        read(file.path(path,"cpu.cfs_quota_us")),read(file.path(path,"cpu.cfs_period_us")))))
      if(length(q)==2L && all(is.finite(q)) && all(q>0)) quotas <- c(quotas,q[1L]/q[2L])
    }
  }
  if(length(quotas)) min(quotas) else NA_real_
}
sby_thread_plan <- function(requested, runtime) {
  requested <- sby_scalar(requested, "sby_config_max_threads", TRUE, -1, .Machine$integer.max)
  if (requested == 0L) stop("sby_config_max_threads deve ser -1L ou um inteiro positivo.", call. = FALSE)
  physical <- parallel::detectCores(logical = FALSE)
  logical <- parallel::detectCores(logical = TRUE)
  affinity <- length(runtime$affinity_cpus)
  quota <- sby_quota()
  limits <- c(if (requested > 0) requested, physical, if (affinity > 0) affinity,
              quota, runtime$openmp_thread_limit)
  limits <- limits[is.finite(limits) & limits > 0]
  effective <- if (length(limits)) max(1L,as.integer(floor(min(limits)))) else 1L
  list(requested = requested, effective = effective, physical_cores = physical,
       logical_cores = logical, affinity_cpus = runtime$affinity_cpus,
       cgroup_quota = quota, openmp_thread_limit=runtime$openmp_thread_limit)
}
sby_domains <- function(x) {
  lapply(x, function(v) {
    lo <- min(v); hi <- max(v)
    kind <- if (all(v %in% c(0, 1))) "binary" else if (all(v == round(v))) "integer" else "continuous"
    list(kind = kind, min = lo, max = hi, storage = typeof(v))
  })
}
sby_restore_domains <- function(syn, domains) {
  out <- as.data.frame(syn)
  for (j in seq_along(domains)) {
    d <- domains[[j]]; v <- out[[j]]
    if (d$kind == "binary") v <- as.double(v >= 0.5)
    if (d$kind == "integer") v <- round(v)
    v <- pmax(d$min, pmin(d$max, v))
    if (d$storage == "integer") v <- as.integer(v)
    out[[j]] <- v
  }
  out
}
sby_formula_input <- function(data, formula) {
  if (!is.data.frame(data) || anyDuplicated(names(data)))
    stop(".data deve ser um data.frame com nomes de colunas únicos.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 3L || !is.symbol(formula[[2L]]))
    stop("Use uma fórmula com um desfecho simples: alvo ~ preditores.", call. = FALSE)
  target <- as.character(formula[[2L]])
  tt <- stats::terms(formula, data = data)
  predictors <- attr(tt, "term.labels")
  # Formula transformations/interactions cannot be copied back as original
  # columns. Reject them instead of silently returning transformed originals.
  predictors <- vapply(predictors, function(v) {
    ex <- tryCatch(str2lang(v), error = function(e) NULL)
    if (!is.symbol(ex)) stop("A fórmula deve selecionar colunas, sem transformações ou interações.", call. = FALSE)
    as.character(ex)
  }, character(1))
  predictors <- setdiff(predictors, target)
  if (!length(predictors) || !all(c(target, predictors) %in% names(data)))
    stop("A fórmula deve selecionar um desfecho existente e ao menos um preditor.", call. = FALSE)
  list(x = as.data.frame(data[predictors]), y = data[[target]], target = target,
       order = names(data)[names(data) %in% c(target, predictors)])
}
sby_local_rng <- function(seed, expr) {
  old_kind <- RNGkind()
  exists_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (exists_seed) old_seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (exists_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(expr)
}
sby_attach <- function(x, summary, detailed) {
  attr(x, "sbyaudit") <- summary
  attr(x, "sby") <- list(synthetic_rows = if (length(summary$adasyn)) summary$adasyn else 0L)
  if (!is.null(detailed)) attr(x, "audit") <- detailed
  x
}
sby_memory_snapshot <- function() {
  if (!file.exists("/proc/self/status")) return(c(rss=NA_real_,peak=NA_real_))
  lines <- tryCatch(readLines("/proc/self/status",warn=FALSE),error=function(e) character())
  get <- function(key) {
    v <- grep(paste0("^",key,":"),lines,value=TRUE)
    if (length(v)!=1L) return(NA_real_)
    as.double(strsplit(trimws(sub("^[^:]+:","",v))," +")[[1L]][1L])*1024
  }
  c(rss=get("VmRSS"),peak=get("VmHWM"))
}
sby_run <- function(x, y, method, params, function_name, ratio_supplied = FALSE) {
  start <- proc.time()[[3L]]
  memory_start <- if(isTRUE(params$sby_audit)) sby_memory_snapshot() else NULL
  if (!is.data.frame(x)) x <- as.data.frame(x)
  if (!ncol(x) || !nrow(x) || !all(vapply(x, function(v) (is.double(v) || is.integer(v)) && !is.object(v), logical(1))))
    stop("Preditores devem ser colunas numéricas simples, não vazias.", call. = FALSE)
  matrix_input <- as.matrix(x)
  if (any(!is.finite(matrix_input))) stop("Preditores não podem conter NA, NaN ou Inf.", call. = FALSE)
  if (length(y) != nrow(x) || anyNA(y) || is.list(y)) stop("Desfecho inválido ou com valores ausentes.", call. = FALSE)
  labels <- if (is.factor(y)) levels(droplevels(y)) else unique(y)
  if (length(labels) != 2L) stop("São necessárias exatamente duas classes observadas.", call. = FALSE)
  codes <- match(y, labels)
  counts <- tabulate(codes, 2L)
  # Stable tie: the first factor level / first observed non-factor label wins.
  rare <- which.min(counts)
  fixed <- params$sby_fixed_minority_label
  if (!is.null(fixed)) {
    rare <- match(fixed, labels)
    if (length(rare) != 1L || is.na(rare)) stop("sby_fixed_minority_label não pertence ao desfecho.", call. = FALSE)
  }
  if (!is.null(params$sby_fixed_majority_label) &&
      !identical(match(params$sby_fixed_majority_label, labels), 3L-rare))
    stop("sby_fixed_majority_label não corresponde à outra classe.", call. = FALSE)
  model <- sby_model(params$nearmiss_model)
  ko <- sby_scalar(params$sby_adasyn_k, "sby_adasyn_k", TRUE, 1, .Machine$integer.max)
  ku <- sby_scalar(params$sby_nearmiss_k, "sby_nearmiss_k", TRUE, 1, .Machine$integer.max)
  km <- sby_scalar(params$sby_nearmiss_m, "sby_nearmiss_m", TRUE, 1, .Machine$integer.max)
  over <- sby_scalar(params$sby_adasyn_ratio, "sby_adasyn_ratio")
  under <- sby_scalar(params$sby_nearmiss_ratio, "sby_nearmiss_ratio")
  seed <- sby_scalar(params$sby_seed, "sby_seed", TRUE, 0, .Machine$integer.max)
  audit <- sby_flag(params$sby_audit, "sby_audit")
  dth <- sby_scalar(params$sby_adasyn_d_th, "sby_adasyn_d_th", FALSE, 0, 1)
  fallback <- match.arg(params$sby_adasyn_zero_difficulty, c("error", "uniform"))
  beta <- params$sby_adasyn_beta
  nmin <- as.double(counts[rare]); nmaj <- as.double(counts[3L-rare])
  if (!is.null(beta)) {
    beta <- sby_scalar(beta, "sby_adasyn_beta", FALSE, 0, 1)
    if (ratio_supplied) stop("Informe sby_adasyn_beta OU sby_adasyn_ratio, não ambos.", call. = FALSE)
  }
  generating <- method %in% c("adasyn", "adanear")
  reducing <- method %in% c("nearmiss", "adanear")
  g <- if (!generating || nmin/nmaj >= dth) 0 else if (!is.null(beta)) floor((nmaj-nmin)*beta) else floor(nmin*over)
  if (!is.finite(g) || g > .Machine$integer.max - nrow(x))
    stop("O número solicitado de sintéticos excede o limite de linhas; reduza sby_adasyn_ratio.", call. = FALSE)
  g <- as.integer(g)
  if (!reducing) under <- 0
  max_rows <- params$sby_max_output_rows
  if (is.null(max_rows)) max_rows <- Inf
  if (!(is.double(max_rows) || is.integer(max_rows)) || length(max_rows) != 1L || is.na(max_rows) || max_rows < 1)
    stop("sby_max_output_rows deve ser positivo ou Inf.", call. = FALSE)
  # Conservative guard considers the expanded input, not merely final retention.
  if (nrow(x) + g > max_rows) stop("sby_max_output_rows seria excedido pela expansão.", call. = FALSE)
  tile <- params$sby_knn_query_chunk_size
  if (is.null(tile)) tile <- 128L
  tile <- sby_scalar(tile, "sby_knn_query_chunk_size", TRUE, 1, .Machine$integer.max)
  tile <- min(tile, 128L)
  runtime <- sby_native("sby_runtime_cpp")
  plan <- sby_thread_plan(params$sby_config_max_threads, runtime)
  if (!is.null(params$sby_memory_guard)) sby_flag(params$sby_memory_guard, "sby_memory_guard")
  gb <- params$sby_max_dense_gb
  if (is.null(gb)) gb <- Inf
  if (!(is.double(gb) || is.integer(gb)) || length(gb) != 1L || is.na(gb) || gb <= 0)
    stop("sby_max_dense_gb deve ser positivo ou Inf.", call. = FALSE)
  neighbor_slots <- max(nmin*min(ko,nrow(x)-1), nmin*max(0,min(ko,nmin-1)),
                        (nmin+g)*min(km,nmaj), nmaj*min(ku,nmin+g))
  scratch_slots <- plan$effective*min(max(ko,ku,km),nrow(x)+g)
  estimate <- 8 * (8 * (nrow(x)+g)*ncol(x) + tile*(nrow(x)+g)) +
              16 * (neighbor_slots + scratch_slots)
  if (!identical(params$sby_memory_guard, FALSE) && estimate > gb*1024^3)
    stop("Estimativa dos buffers densos excede sby_max_dense_gb; reduza a expansão ou aumente o orçamento.", call. = FALSE)
  scale <- params$sby_scaling_info
  if (is.null(scale)) scale <- params$sby_precomputed_scaling
  if (!is.null(scale)) {
    if (!is.list(scale)) stop("Scaling information must be a list of centers and scales.",call.=FALSE)
    if (is.null(scale$sby_center) && !is.null(scale$centers)) scale$sby_center <- scale$centers
    if (is.null(scale$sby_scale) && !is.null(scale$scales)) scale$sby_scale <- scale$scales
    if (is.null(scale$sby_center) && !is.null(scale$means)) scale$sby_center <- scale$means
    if (is.null(scale$sby_scale) && !is.null(scale$sds)) scale$sby_scale <- scale$sds
    if (!(is.double(scale$sby_center) || is.integer(scale$sby_center)) ||
        !(is.double(scale$sby_scale) || is.integer(scale$sby_scale)) ||
        !is.null(dim(scale$sby_center)) || !is.null(dim(scale$sby_scale)) ||
        length(scale$sby_center) != ncol(x) || length(scale$sby_scale) != ncol(x) ||
        any(!is.finite(scale$sby_center)) || any(!is.finite(scale$sby_scale)) || any(scale$sby_scale <= 0))
      stop("Escala fornecida requer sby_center e sby_scale finitos, um valor por coluna, com desvios positivos.", call. = FALSE)
  }
  if (!is.null(params$sby_input_already_scaled)) {
    sby_flag(params$sby_input_already_scaled, "sby_input_already_scaled")
    if (params$sby_input_already_scaled) {
      if (is.null(scale)) stop("Dados já padronizados requerem sby_scaling_info ou sby_precomputed_scaling.", call. = FALSE)
      matrix_input <- sweep(sweep(matrix_input, 2, scale$sby_scale, `*`), 2, scale$sby_center, `+`)
      x <- as.data.frame(matrix_input)
    }
  }
  domains <- if(g>0L || audit) sby_domains(x) else NULL
  opts <- list(threads=plan$effective, rare=as.integer(rare), synthetics=g,
               model=model, ko=ko, ku=ku, km=km, tile=as.integer(tile), audit=audit,
               uniform=fallback == "uniform", under=under)
  if (!is.null(scale)) {
    opts$centers <- as.double(scale$sby_center); opts$scales <- as.double(scale$sby_scale)
  }
  storage.mode(matrix_input) <- "double"
  memory_native_start <- if(audit) sby_memory_snapshot() else NULL
  native_start <- proc.time()[[3L]]
  result <- sby_local_rng(seed, sby_native("sby_pipeline_cpp", matrix_input, as.integer(codes), opts))
  native_end <- proc.time()[[3L]]
  memory_native_end <- if(audit) sby_memory_snapshot() else NULL
  if (result$candidate_shortage)
    warning("NearMiss-3: candidatos insuficientes para a meta; todos os candidatos foram mantidos, sem completar com não candidatos.", call. = FALSE)
  colnames(result$synthetic) <- names(x)
  if (audit) colnames(result$synthetic_z) <- names(x)
  syn <- if(g>0L) sby_restore_domains(result$synthetic, domains) else x[FALSE,,drop=FALSE]
  names(syn) <- names(x)
  # Output consists only of copied originals and domain-restored synthetics.
  orig <- sort(c(result$minority_indices, result$majority_indices))
  xx <- rbind(x[orig, , drop=FALSE], syn)
  rownames(xx) <- NULL
  yy <- c(y[orig], rep(y[result$minority_indices[1L]], g))
  synth_rows <- if (g) seq.int(length(orig)+1L, length.out=g) else integer()
  majority_rows <- match(result$majority_indices, orig)
  actual <- params[!names(params) %in% c(".data", "sby_data", "sby_x_matrix", "sby_y_vector")]
  actual$input <- list(rows=nrow(x), predictor_columns=names(x),
                       target_class=class(y), class_labels=labels,
                       target_levels=if(is.factor(y)) levels(y) else NULL)
  actual$sby_seed <- seed
  actual$sby_adasyn_zero_difficulty <- fallback
  actual$effective <- list(synthetic_count=g, initial_minority_label=labels[rare],
    initial_majority_label=labels[3L-rare], expansion_ratio=if(nmin) g/nmin else 0,
    beta=if(nmaj>nmin) g/(nmaj-nmin) else 0,
    neighbor_counts=stats::setNames(result$effective_k, c("difficulty", "interpolation", "nearmiss", "preselection")),
    thread_limit=plan$effective, tile_rows=tile)
  summary <- list(adasyn=as.integer(synth_rows), nearmiss=as.integer(majority_rows),
    package="sbyadanear", function_name=function_name,
    original_majority_indices=result$majority_indices,
    original_minority_indices=result$minority_indices)
  summary[[function_name]] <- actual
  detailed <- NULL
  scale_out <- list(sby_center=stats::setNames(result$centers,names(x)),
                    sby_scale=stats::setNames(result$scales,names(x)), denominator=nrow(x))
  if (audit) {
    telemetry <- tibble::as_tibble(result$telemetry)
    overhead <- function(stage, elapsed, before, after) tibble::tibble(stage=stage, elapsed_seconds=elapsed,
      rss_before_bytes=before[["rss"]], rss_after_bytes=after[["rss"]], process_peak_rss_bytes=after[["peak"]],
      omp_observed_threads=1L, thread_limit=plan$effective,
      mkl_configured_threads=NA_integer_, mkl_observed_threads=NA_integer_)
    transfer <- max(0,native_end-native_start-sum(telemetry$elapsed_seconds))
    memory_transfer <- c(rss=utils::tail(telemetry$rss_after_bytes,1),
                          peak=utils::tail(telemetry$process_peak_rss_bytes,1))
    telemetry <- rbind(overhead("validation_and_preparation",native_start-start,memory_start,memory_native_start),
      telemetry, overhead("native_transfer_and_return",transfer,memory_transfer,memory_native_end),
      overhead("domain_restoration_and_assembly",proc.time()[[3L]]-native_end,memory_native_end,sby_memory_snapshot()))
    telemetry$physical_cores <- plan$physical_cores
    telemetry$logical_cores <- plan$logical_cores
    telemetry$affinity_cpu_count <- length(plan$affinity_cpus)
    telemetry$cgroup_cpu_quota <- plan$cgroup_quota
    telemetry$estimated_dense_buffer_bytes <- estimate
    telemetry$output_object_bytes <- as.double(utils::object.size(xx) + utils::object.size(yy))
    detailed <- list(package="sbyadanear",function_name=function_name,initial_scale=scale_out, domains=domains,
      original_indices=orig, original_minority_indices=result$minority_indices,
      original_majority_indices=result$majority_indices,
      parameters=actual, counts=list(input=nrow(x), minority=nmin, majority=nmaj,
        synthetic=g, retained_majority=length(result$majority_indices), output=nrow(xx)),
      nearmiss=list(model=model, executed=result$nearmiss_executed,
                    candidates=result$candidates, scores=result$scores,
                    candidate_shortage=result$candidate_shortage),
      adasyn=list(executed=result$adasyn_executed,
                  difficulty_hits=result$difficulty_hits, quotas=result$quotas,
                  uniform_fallback=result$uniform_fallback),
      synthetics=tibble::tibble(output_row=synth_rows, parent=result$parents,
                                neighbor=result$partners, lambda=result$lambda),
      synthetic_continuous=result$synthetic, synthetic_standardized=result$synthetic_z,
      threads=c(plan,list(runtime=runtime)), telemetry=telemetry)
  }
  list(x=xx, y=yy, summary=summary, audit=detailed, scaling=scale_out,
       original_indices=orig, initial_counts=counts, rare=rare, labels=labels)
}

# Legacy KNN selectors are accepted only where they name the exact geometry.
# Approximate or alternate distances cannot silently redefine these methods.
sby_legacy_options <- function(p) {
  for(name in c("sby_return_scaled","sby_return_original_scale","sby_return_index",
                "sby_return_scaling_info","sby_return_reduced_scaled",
                "sby_input_already_scaled","sby_memory_guard","sby_audit"))
    if (!is.null(p[[name]])) sby_flag(p[[name]],name)
  if (!is.null(p$sby_knn_workers)) {
    workers <- sby_scalar(p$sby_knn_workers,"sby_knn_workers",TRUE,-1,.Machine$integer.max)
    if(workers==0L) stop("sby_knn_workers must be -1L or a positive integer.",call.=FALSE)
  }
  one <- function(v) if (length(v)>1L) v[1L] else v
  if (!is.null(p$sby_knn_engine) && !one(p$sby_knn_engine) %in% c("auto","native"))
    stop("Esta versão requer sby_knn_engine = 'auto' ou 'native' (Intel oneAPI exato).", call. = FALSE)
  if (!is.null(p$sby_knn_algorithm) && !one(p$sby_knn_algorithm) %in% c("auto","brute"))
    stop("Use sby_knn_algorithm = 'auto' ou 'brute' para o motor exato oneAPI.", call. = FALSE)
  if (!is.null(p$sby_knn_distance_metric) && one(p$sby_knn_distance_metric) != "euclidean")
    stop("ADASYN/NearMiss neste pacote requerem distância euclidiana.", call. = FALSE)
  if (!is.null(p$sby_type_info))
    stop("sby_type_info externo não é aceito; forneça dados originais para inferir os domínios.", call. = FALSE)
  if (!is.null(p$sby_restore_types) && !sby_flag(p$sby_restore_types,"sby_restore_types"))
    stop("sby_restore_types deve ser TRUE para preservar os domínios originais.", call. = FALSE)
  if (!is.null(p$sby_knn_parallel_backend) && !one(p$sby_knn_parallel_backend) %in% c("parallel","RcppParallel"))
    stop("Backend legado inválido; o processamento utiliza Intel OpenMP/oneMKL.", call. = FALSE)
  if (!is.null(p$sby_knn_hnsw_m) && p$sby_knn_hnsw_m != 16L ||
      !is.null(p$sby_knn_hnsw_ef) && p$sby_knn_hnsw_ef != 200L)
    stop("Parâmetros HNSW não se aplicam ao motor exato; remova as configurações HNSW.", call. = FALSE)
}
sby_dispatch <- function(method, p, function_name, ratio_supplied) {
  sby_legacy_options(p)
  hpc <- grepl("_hpc$", function_name)
  matrix_api <- grepl("_matrix$|_index$", function_name)
  index_api <- grepl("_index$", function_name)
  if (!is.null(p$sby_strategy)) {
    strategy <- match.arg(p$sby_strategy,c("none","weight","adasyn","nearmiss","adanear","adanearWeight"))
    method <- switch(strategy, none="none", weight="none", adanearWeight="adanear", strategy)
  }
  if (matrix_api) {
    if (!is.matrix(p$sby_x_matrix) || !is.numeric(p$sby_x_matrix))
      stop("sby_x_matrix deve ser uma matriz numérica densa.", call. = FALSE)
    input <- list(x=as.data.frame(p$sby_x_matrix),y=p$sby_y_vector)
  } else {
    input <- sby_formula_input(if(hpc) p$.data else p$sby_data,
                               if(hpc) p$formula else p$sby_formula)
  }
  defaults <- list(sby_adasyn_ratio=0, sby_nearmiss_ratio=0, sby_seed=1L,
                   sby_adasyn_k=5L,sby_nearmiss_k=5L,sby_nearmiss_m=3L,
                   nearmiss_model=3L,sby_config_max_threads=-1L,sby_audit=FALSE,
                   sby_adasyn_beta=NULL,sby_adasyn_d_th=1,sby_adasyn_zero_difficulty="error")
  for (name in names(defaults)) if (!name %in% names(p)) p[name] <- defaults[name]
  p$sby_config_max_threads <- sby_scalar(p$sby_config_max_threads,"sby_config_max_threads",TRUE,-1,.Machine$integer.max)
  if (!hpc && p$sby_config_max_threads == -1L) p$sby_config_max_threads <- p$sby_knn_workers
  if (is.null(p$sby_config_max_threads)) p$sby_config_max_threads <- -1L
  if (!is.null(p$sby_knn_over_k)) p$sby_adasyn_k <- p$sby_knn_over_k
  if (!is.null(p$sby_knn_under_k)) p$sby_nearmiss_k <- p$sby_knn_under_k
  if (!is.null(p$sby_audit_level)) {
    level <- match.arg(p$sby_audit_level,c("none","light","full"))
    if (level != "none") p$sby_audit <- TRUE
  }
  out <- sby_run(input$x,input$y,method,p,function_name,ratio_supplied)
  scaled <- function(xx) sweep(sweep(as.matrix(xx),2,out$scaling$sby_center,`-`),2,out$scaling$sby_scale,`/`)
  if (matrix_api) {
    if (index_api) {
      value <- list(sby_retained_index=out$original_indices,
        sby_selected_majority_index=out$summary$original_majority_indices)
      if (isTRUE(p$sby_return_scaling_info)) value$sby_scaling_info <- out$scaling
      if (isTRUE(p$sby_return_reduced_scaled)) value$sby_reduced_scaled <- scaled(out$x)
    } else {
      value <- list(sby_x_matrix=if(identical(p$sby_return_original_scale,FALSE)) scaled(out$x) else as.matrix(out$x),
                    sby_y_vector=out$y,sby_scaling_info=out$scaling)
      if (isTRUE(p$sby_return_index)) value$sby_retained_index <- out$original_indices
      if (isTRUE(p$sby_return_scaled)) value$sby_x_scaled <- scaled(out$x)
    }
    final_counts <- tabulate(match(out$y,out$labels),2L)
    value$sby_input_class_distribution <- stats::setNames(out$initial_counts,out$labels)
    value$sby_output_class_distribution <- stats::setNames(final_counts,out$labels)
    value$sby_class_ratio_input <- out$initial_counts[3L-out$rare]/out$initial_counts[out$rare]
    value$sby_class_ratio_output <- final_counts[3L-out$rare]/final_counts[out$rare]
    value$sby_diagnostics <- list(sby_method=method,sby_input_rows=nrow(input$x),
                                  sby_output_rows=nrow(out$x),sby_knn_engine="Intel oneAPI exact")
    return(sby_attach(value,out$summary,out$audit))
  }
  value <- tibble::as_tibble(out$x)
  value[[input$target]] <- out$y
  value <- value[input$order]
  value <- sby_attach(value,out$summary,out$audit)
  if (hpc || (!isTRUE(p$sby_audit) && !isTRUE(p$sby_return_scaled))) return(value)
  # Existing tabular audit return remains a list; recipes always return data.
  wrapper <- list(sby_balanced_data=value, sby_scaling_info=out$scaling)
  if (isTRUE(p$sby_return_scaled)) wrapper$sby_x_scaled <- scaled(out$x)
  sby_attach(wrapper,out$summary,out$audit)
}
