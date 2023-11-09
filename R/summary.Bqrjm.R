#' @export
#'
summary.Bqrjm <- function (object, ...)
{
  if (!inherits(object, "Bqrjm"))
    stop("Use only with \"Bqrjm\" objects.\n")

  #---- Global details
  cat("#-- Statistical model: Quantile regression joint model \n")
  if(object$control$param=="sharedRE")
    cat(paste("     - Association structure: Shared random effect(s) \n"))
  if(object$control$param=="value")
    cat(paste("     - Association structure: Current quantile value of longitudinal process \n"))
  if(object$control$survMod=="weibull")
    cat(paste("     - Survival baseline risk function: Weibull distribution \n"))
  cat(paste("     - Quantile order(s): ", object$control$tau, "\n"))
  cat(paste("     - Number of observations: ", nrow(object$data), "\n"))
  cat(paste("     - Number of statistic units (e.g. subject): ", object$control$I, "\n"))
  cat(paste("     - Number of observed events: ", sum(object$control$event), "\n"))
  cat("\n")

  #---- Parameter Estimations
  coefs <- object$mean
  CIs <- object$CIs
  Rhat <- object$Rhat
  # beta regression parameters
  beta_estim <- cbind("Value" = coefs$beta,
                      "2.5%" = CIs$beta[, 1],
                      "97.5%" = CIs$beta[, 2],
                      "Rhat" = Rhat$beta)
  sigma_estim <- cbind("Value" = coefs$sigma,
                        "2.5%" = CIs$sigma[1],
                        "97.5%" = CIs$sigma[2],
                        "Rhat" = Rhat$sigma)
  rownames(sigma_estim) <- "sigma"

  cat("#-- Estimation of longitudinal regression parameters and their credible interval bounds: \n")
  prmatrix(beta_estim, na.print = "")
  cat("\n")
  cat("#-- Estimation of 'sigma' parameter associated with asymmetric Laplace distribution: \n")
  prmatrix(sigma_estim, na.print = "")

  # Random effects for "mixed regression model
  cat("\n")
  cat("#-- (Co)variance matrix of the random-effect(s): \n")
  if(object$control$RE_ind){
    tmp <- diag(object$mean$covariance.b)
    colnames(tmp) <- rownames(tmp) <- names(object$mean$covariance.b)
    prmatrix(tmp, na.print = "")
  }
  if(!object$control$RE_ind)
    prmatrix(object$mean$covariance.b, na.print = "")

  # survival parameters
  # alpha regression parameters# survival structure parameter
  if(object$control$survMod=="weibull"){
    param_estim1 <- cbind("Value" = object$mean$shape,
                          "2.5%"  = object$CIs$shape[1],
                          "97.5%" = object$CIs$shape[2],
                          "Rhat" = object$Rhat$shape)
    rownames(param_estim1) <- "shape"
  }
  # alpha parameters
  param_estim2 <- cbind("Value" = object$mean$alpha,
                        "2.5%" = object$CIs$alpha[, 1],
                        "97.5%" = object$CIs$alpha[, 2],
                        "Rhat" = object$Rhat$alpha)
  # association parameters
  if(length(object$mean$alpha.assoc)==1)
    param_estim3 <- cbind("Value" = object$mean$alpha.assoc,
                          "2.5%" = object$CIs$alpha.assoc[1],
                          "97.5%" = object$CIs$alpha.assoc[2],
                          "Rhat" = object$Rhat$alpha.assoc)
  else
    param_estim3 <- cbind("Value" = object$mean$alpha.assoc,
                          "2.5%" = object$CIs$alpha.assoc[, 1],
                          "97.5%" = object$CIs$alpha.assoc[, 2],
                          "Rhat" = object$Rhat$alpha.assoc)
  rownames(param_estim3) <- rep("alpha.assoc", length(object$mean$alpha.assoc))
  # Print
  cat("\n")
  cat("#-- Estimation of survival models: \n")
  prmatrix(rbind(param_estim1, param_estim2, param_estim3), na.print = "")

}
