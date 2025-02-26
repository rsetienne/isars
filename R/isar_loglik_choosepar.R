isar_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix, ...)
{
  trpars1 <- rep(0,length(idparsopt) + length(idparsfix))
  trpars1[idparsopt] <- trparsopt
  if(length(idparsfix) != 0)
  {
    trpars1[idparsfix] <- trparsfix
  }
  if(max(trpars1) > 1 || min(trpars1) < 0)
  {
    loglik <- -Inf
  } else {
    pars1 <- trpars1/(1 - trpars1)
    loglik <- isar_loglik(pars = pars1, ...)
    if(is.nan(loglik) || is.na(loglik))
    {
      warning("Parameter values have been used that cause numerical problems.\n")
      loglik <- -Inf
    }
  }
  return(loglik)
}
