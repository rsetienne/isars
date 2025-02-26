#' Corrected Akaike Information Criterion
#' @inherit default_params_doc
#' @return The AICc value
#' @author Rampal S. Etienne

AICc <- function(LogLik,k,n){
  aicc <- AIC(LogLik,k) + ((2 * k * (k + 1))/(n - k - 1))
  return(aicc)
}
