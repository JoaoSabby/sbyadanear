# Recipe state uses the recipes trained/skip/id fields. Audit never changes
# the bake return type, so downstream recipes steps still receive data.
sby_step_add <- function(recipe, terms, parameters, public, method, role,
                         trained, columns, skip, id, ratio_supplied) {
  sby_flag(skip, "skip"); sby_flag(trained,"trained")
  sby_legacy_options(parameters)
  if (!is.null(parameters$nearmiss_model)) sby_model(parameters$nearmiss_model)
  if (!is.null(parameters$sby_adasyn_beta) && ratio_supplied)
    stop("Informe sby_adasyn_beta OU sby_adasyn_ratio, não ambos.",call.=FALSE)
  # Force the seed now: each training/bake call uses the stored integer.
  parameters$sby_seed <- sby_scalar(parameters$sby_seed,"sby_seed",TRUE,0,.Machine$integer.max)
  object <- recipes::step(subclass=public,
    terms=terms,role=role,trained=trained,columns=columns,
    parameters=parameters,public=public,method=method,
    ratio_supplied=ratio_supplied,audit_log=new.env(parent=emptyenv()),skip=skip,id=id)
  recipes::add_step(recipe,object)
}
sby_recipe_prep <- function(x, training, info, ...) {
  selected <- recipes::recipes_eval_select(x$terms, training, info)
  if (length(selected) != 1L) stop("A etapa deve selecionar exatamente um desfecho.",call.=FALSE)
  x$columns <- names(selected)
  x$trained <- TRUE
  x
}
sby_recipe_bake <- function(object, new_data, ...) {
  if (!object$trained) stop("Use prep() antes de bake().",call.=FALSE)
  if (!object$columns %in% names(new_data)) stop("Desfecho ausente em new_data.",call.=FALSE)
  hpc <- grepl("_hpc$",object$public)
  p <- object$parameters
  formula <- stats::as.formula(call("~",as.name(object$columns),quote(.)))
  if (hpc) { p$.data <- new_data; p$formula <- formula }
  else { p$sby_data <- new_data; p$sby_formula <- formula }
  name <- sub("_step_","_",object$public,fixed=TRUE)
  result <- sby_dispatch(object$method,p,name,object$ratio_supplied)
  if (is.list(result) && !is.data.frame(result)) result <- result$sby_balanced_data
  # Preserve public step provenance rather than advertise only its delegate.
  sm <- attr(result,"sbyaudit")
  sm[[object$public]] <- sm[[sm$function_name]]
  sm[[sm$function_name]] <- NULL
  sm$function_name <- object$public
  attr(result,"sbyaudit") <- sm
  if (!is.null(attr(result,"audit"))) {
    a <- attr(result,"audit")
    a$function_name <- object$public
    attr(result,"audit") <- a
  }
  # recipes::bake(recipe) may strip attributes in its final column selection.
  # Keep the last invocation's audit accessible through the stored step too.
  object$audit_log$last <- list(sbyaudit=sm,audit=attr(result,"audit"))
  result
}
sby_recipe_print <- function(x, width=max(20, getOption("width")-30), ...) {
  cat(x$public, ": ", if(x$trained) paste(x$columns,collapse=", ") else "não treinado", "\n",sep="")
  invisible(x)
}
sby_recipe_tidy <- function(x, ...) tibble::tibble(terms=if(x$trained) x$columns else vapply(x$terms,rlang::as_label,character(1)), id=x$id)
sby_recipe_required <- function(x, ...) c("sbyadanear","recipes","tibble")
