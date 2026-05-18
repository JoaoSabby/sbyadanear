#' Resolver nivel de auditoria matricial
#'
#' @param sby_audit Indicador legado de auditoria completa
#' @param sby_audit_level Nivel none, light ou full
#'
#' @return String com nivel efetivo de auditoria
#' @noRd
sby_resolve_audit_level <- function(sby_audit, sby_audit_level = c("none", "light", "full")){
  sby_audit <- sby_validate_logical_scalar(sby_audit, "sby_audit")
  sby_audit_level <- match.arg(sby_audit_level)
  if(isTRUE(sby_audit)){
    return("full")
  }
  return(sby_audit_level)
}
