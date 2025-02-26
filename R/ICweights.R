#' Information Criterion Weights
#' @inherit default_params_doc
#' @return The IC weights for a set of models for which loglikelihoods are provided.
#' @author Rampal S. Etienne

ICweights <- function(IC) {
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights)
}
