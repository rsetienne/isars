#' Maximum Likelihood of isar model
#' @name isar_M
#' @title Computes the loglikelihood of a parameter set of a specified island
#' species-area relationships (isars) for a given data set of values of
#' area and richness.
#' @description Computes the logarithm of the likelihood of a model for the
#' island species-area with a given set of parameters for a given data set with
#' values of area and richness for several islands. It assumes that the island
#' species-area relationship provides the expected value of richness, and the
#' error distribution is Poisson.
#' @inheritParams default_params_doc
#' @return The maximum likelhood parameters and corresponding likelihood and
#' AICc value.
#' @author Rampal S. Etienne
#' @export
isar_ML <- function(f_isar,
                    initparsopt,
                    parsfix,
                    idparsopt,
                    idparsfix,
                    area,
                    obs_richness,
                    optimmethod = 'simplex',
                    num_cycles = 5,
                    tol = c(1E-6, 1E-6, 1E-6),
                    maxiter = 10000 * round((1.25)^length(idparsopt)),
                    trial_settings = c(n = 10, sd = 0.1)) {
  optimpars <- c(tol, maxiter)
  ML <- -Inf
  for(i in 0:trial_settings[1]) {
    if(i > 0) {
      initparsopt2 <- initparsopt * exp(rnorm(n = length(initparsopt), mean = 0,sd = trial_settings[2]))
      if(length(parsfix) > 0) {
        parsfix2 <- parsfix * exp(rnorm(n = length(parsfix), mean = 0,sd = trial_settings[2]))
      }
    } else {
      initparsopt2 <- initparsopt
      parsfix2 <- parsfix
    }
    trparsopt <- initparsopt2/(1 + initparsopt2)
    trparsfix <- parsfix2/(1 + parsfix2)
    trparsfix[parsfix2 == Inf] = 1
    initloglik <- isar_loglik_choosepar(trparsopt = trparsopt,
                                        trparsfix = trparsfix,
                                        idparsopt = idparsopt,
                                        idparsfix = idparsfix,
                                        area = area,
                                        obs_richness = obs_richness,
                                        f_isar = f_isar)
    if(initloglik > -Inf) {
      out2 <- DDD::optimizer(optimmethod = optimmethod,
                             optimpars = optimpars,
                             fun = isar_loglik_choosepar,
                             trparsopt = trparsopt,
                             trparsfix = trparsfix,
                             idparsopt = idparsopt,
                             idparsfix = idparsfix,
                             area = area,
                             obs_richness = obs_richness,
                             num_cycles = num_cycles,
                             f_isar = f_isar)
      if(as.numeric(unlist(out2$fvalues)) > ML) {
        out <- out2
      }
    } else {
      cat("The initial parameter values lead to -Inf for the loglikelihood. Try again with different initial values.\n")
      out2 <- list(conv = 13)
    }
  }
  if(out$conv > 0)
  {
    cat("Optimization has not converged. Try again with different initial values.\n")
    out2 <- NULL
  } else {
    MLtrpars <- as.numeric(unlist(out$par))
    MLpars <- MLtrpars/(1 - MLtrpars)
    MLpars1 <- rep(0,length(idparsopt) + length(idparsfix))
    MLpars1[idparsopt] <- MLpars
    MLpars1[idparsfix] <- parsfix
    ML <- as.numeric(unlist(out$fvalues))
    AIC_c <- AICc(LogLik = ML, k = length(initparsopt), n = length(obs_richness))
    out2 <- c(MLpars1, ML, AIC_c)
    names(out2)[idparsopt] <- names(initparsopt)
    names(out2)[idparsfix] <- names(parsfix)
    names(out2)[(length(out2) - 1):length(out2)] <- c("ML","AICc")
  }
  return(out2)
}
