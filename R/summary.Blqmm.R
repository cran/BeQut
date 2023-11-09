#' @export
#'
summary.Blqmm <- function (object, ...)
{
  if (!inherits(object, "Blqmm"))
    stop("Use only with \"Blqmm\" objects.\n")

  #---- Global details
  cat("#-- Statistical model: Linear quantile mixed model \n")
  cat(paste("     - Quantile order: ", object$control$tau, "\n"))
  cat(paste("     - Number of observations: ", nrow(object$data), "\n"))
  cat(paste("     - Number of statistic units (e.g. subject): ", object$control$I, "\n"))
  cat("\n")

  #---- Parameter Estimations
  coefs <- object$mean
  CIs <- object$CIs
  Rhat <- object$Rhat
  # beta regression parameters
  beta_estim <- cbind(
    "Value" = coefs$beta,
    "2.5%" = CIs$beta[, 1],
    "97.5%" = CIs$beta[, 2],
    "Rhat" = Rhat$beta
  )
  sigma_estim <- cbind(
    "Value" = coefs$sigma,
    "2.5%" = CIs$sigma[1],
    "97.5%" = CIs$sigma[2],
    "Rhat" = Rhat$sigma
  )
  rownames(sigma_estim) <- "sigma"
  cat("#-- Estimation of longitudinal regression parameters and their credible interval bounds: \n")
  prmatrix(beta_estim, na.print = "")
  cat("\n")
  cat("#-- Estimation of 'sigma' parameter associated with asymmetric Laplace distribution: \n")
  prmatrix(sigma_estim, na.print = "")

  # Random effects for mixed regression model
  cat("\n")
  cat("#-- (Co)variance matrix of the random-effect(s): \n")
  if (object$control$RE_ind) {
    tmp <- diag(object$mean$covariance.b)
    colnames(tmp) <-
      rownames(tmp) <- names(object$mean$covariance.b)
    prmatrix(tmp, na.print = "")
  }
  if (!object$control$RE_ind)
    prmatrix(object$mean$covariance.b, na.print = "")

}
