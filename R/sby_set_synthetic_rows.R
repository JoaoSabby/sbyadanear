# Attach metadata identifying the synthetic observations in a returned object.
sby_set_synthetic_rows <- function(sby_result, sby_synthetic_rows = integer()){
  sby_synthetic_rows <- as.integer(sby_synthetic_rows)
  if(length(sby_synthetic_rows) == 0L){
    sby_synthetic_rows <- 0L
  }
  attr(sby_result, "sby") <- list(synthetic_rows = sby_synthetic_rows)
  sby_result
}
