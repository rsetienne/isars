#' Default parameter documentation
#'
#' @param pars A numeric vector containing the model parameters:
#'   \itemize{
#'     \item{\code{pars[1]}: }
#'     \item{\code{pars[2]}: }
#'   }
#' @param initparsopt The initial values of the parameters that must be
#'   optimized; they are all positive.
#' @param idparsopt The ids of the parameters that must be optimized.
#' @param idparsfix The ids of the parameters that should not be optimized,
#' @param LogLik The loglikelhood or a vector of loglikelihood values
#' @param k The number of parameters
#' @param n The number of observations (sample size)
#' @param IC The value(s) of the Information Criterion
#' @param parsfix The values of the parameters that should not be optimized.
#' @param tol Sets the tolerances in the optimization. Consists of: \cr reltolx
#'   = relative tolerance of parameter values in optimization \cr reltolf =
#'   relative tolerance of function value in optimization \cr abstolx = absolute
#'   tolerance of parameter values in optimization.
#' @param maxiter Sets the maximum number of iterations in the optimization.
#' @param optimmethod Method used in likelihood optimization. For example
#'  \code{simplex} or \code{subplex}. See \code{\link[DDD]{optimizer}()} for
#'  further details.
#' @param jitter Numeric for \code{\link[DDD]{optimizer}()}. Jitters the
#'   parameters being optimized by the specified amount which should be very
#'   small, e.g. 1e-5. Jitter when \code{link{subplex}{subplex}()} produces
#'   incorrect output due to parameter transformation.
#' @param num_cycles The number of cycles the optimizer will go through.
#'   Default is 1.
#' @param f_isar Function for island species-area relationship
#' @param area Vector of values of areas for different islands
#' @param obs_richness vector of values of observed richness on islands. Should
#' be the same order as in area.
#' @param trial_settings Number of trials in optimization, seed of the pseudo-
#' random number generator, and standard deviation of the normal distribution
#' used to generate the trials (initial parameter sets).
#' @param working_directory Name of working directory (as a string) where results
#' result files will be put. Default is current working directory.
#' @param verbose Determines the level of verbosity in the output.
#' @param type_of_richness Type of richness that is used. Choice is one of
#' "D1.OC", "D2.OC_Int", "D3.OC_Ext", "D6.CC".
#' D1.OC: the diversity/richness from the Original Community
#' D2.OC_Int: Original community richness plus introduced species
#' D3.OC_Ext: Original community minus extinct species
#' D6.CC: Contemporary Community, that is, OC plus Introduced species minus
#' Extinct species.
#' @return Nothing
default_params_doc <- function(
    initparsopt,
    idparsopt,
    idparsfix,
    parsfix,
    tol,
    maxiter,
    optimmethod,
    jitter,
    num_cycles,
    f_isar,
    area,
    obs_richness,
    trial_settings,
    num_init,
    working_directory,
    verbose,
    type_of_richness
) {
  # Nothing
}
