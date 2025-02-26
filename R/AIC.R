AIC <- function(LogLik,k){
  aic <- (2 * k) - ( 2 * LogLik)
  return(aic)
}
