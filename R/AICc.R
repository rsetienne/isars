AICc <- function(LogLik,k,n){
  aicc <- AIC(LogLik,k) + ((2 * k * (k + 1))/(n - k - 1))
  return(aicc)
}
