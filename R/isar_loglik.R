isar_loglik <- function(f_isar,
                        pars,
                        area,
                        obs_richness) {
  exp_richness <- f_isar(area, pars)
  loglik <- sum(dpois(x = obs_richness, lambda = exp_richness, log = TRUE))
  return(loglik)
}
