#' @export
#'
summary.Blqm <- function (object, ...)
{
  if (!inherits(object, "Blqm"))
    stop("Use only with \"Blqm\" objects.\n")

  #---- Global details
  cat("#-- Statistical model: Linear quantile model \n")
  cat(paste("     - Quantile order: ", object$control$tau, "\n"))
  cat(paste("     - Number of observations: ", nrow(object$data), "\n"))
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
  cat("#-- Estimation of linear predictor parameters and their credible interval bounds: \n")
  prmatrix(beta_estim, na.print = "")
  cat("\n")
  cat(
    "#-- Estimation of 'sigma' parameter associated with asymmetric Laplace distribution: \n"
  )
  prmatrix(sigma_estim, na.print = "")

}
