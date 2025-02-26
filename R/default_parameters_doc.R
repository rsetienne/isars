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
) {
  # Nothing
}
