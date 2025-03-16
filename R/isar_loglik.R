#' Likelihood of isar model
#' @name isar_loglik
#' @title Computes the loglikelihood of a parameter set of a specified island
#' species-area relationships (isars) for a given data set of values of
#' area and richness.
#' @description Computes the logarithm of the likelihood of a model for the
#' island species-area with a given set of parameters for a given data set with
#' values of area and richness for several islands. It assumes that the island
#' species-area relationship provides the expected value of richness, and the
#' error distribution is Poisson.
#' @inheritParams default_params_doc
#' @return The loglikelihood.
#' @author Rampal S. Etienne
#' @export
isar_loglik <- function(f_isar,
                        pars,
                        area,
                        obs_richness) {
  exp_richness <- f_isar(area, pars)
  loglik <- sum(dpois(x = obs_richness, lambda = exp_richness, log = TRUE))
  return(loglik)
}
