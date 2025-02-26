#' Akaike Information Criterion
#' @inherit default_params_doc
#' @return The AIC value
#' @author Rampal S. Etienne

AIC <- function(LogLik,k){
  aic <- (2 * k) - (2 * LogLik)
  return(aic)
}
