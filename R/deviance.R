#' \code{deviance} returns the deviance based on the conditional likelihood associated with the survival part.
#'
#' @param object an object inheriting from class 'Bqrjm'.
#' @param M an integer indicating the number of draws used for the approximation of the integral with respect to random effects, M=1000 by default.
#' @param conditional is "survival" by default (only this one is implemented until now).
#' @param verbose A logical indicating if information about method's progress (included progress bars for each step) must be printed (default to TRUE). Adds a small extra overload.
#'
#' @return An object which is a list with the following elements:
#'    \describe{
#'   \item{\code{deviance}}{Numerical object returning the deviance}
#'   \item{\code{likelihood}}{(Conditional) likelihood}
#'   \item{\code{sims.list}}{list of individual quantities as likelihood, draws of random effects, hazard and survival functions}
#'   \item{\code{control}}{list of arguments giving details about the deviance}
#'  }
#'
#' @author Antoine Barbieri and Baptiste Courrèges
#'
#' @import MASS survival
#' @export
#'
#' @examples
#'
#' \donttest{
#' #---- load data
#' data(dataLong)
#'
#' #---- Fit quantile regression joint model for the median
#' qrjm_50 <- qrjm(formFixed = y ~ visit,
#'                formRandom = ~ visit,
#'                formGroup = ~ ID,
#'                formSurv = survival::Surv(time, event) ~ X1 + X2,
#'                survMod = "weibull",
#'                param = "value",
#'                timeVar= "visit",
#'                data = dataLong,
#'                save_va = TRUE,
#'                parallel = FALSE,
#'                tau = 0.5)
#'
#' deviance(qrjm_50, M=200)
#' }
#'
deviance <- function(object, M=1000, conditional="survival", verbose = TRUE){

  # # To do :
  # - Add the shared random effect case

  # condition on the object
  if(!inherits(object, "Bqrjm"))
    stop("Use only with 'Bqrjm' objects from 'qrjm()' function.\n")
  if(is.null(object$out_jagsUI$q50$va1))
    stop("Saving posterior of auxilary variable 'W' is needed.\n")

  # get and management of data

  # longitudinal
  data <- object$data[unique(c(all.vars(object$control$formGroup),all.vars(object$control$formFixed),all.vars(object$control$formRandom)))]
  y <- data[all.vars(object$control$formFixed)][, 1]
  mfX <- stats::model.frame(object$control$formFixed, data = data)
  X <- stats::model.matrix(object$control$formFixed, mfX)
  mfU <- stats::model.frame(object$control$formRandom, data = data)
  U <- stats::model.matrix(object$control$formRandom, mfU)
  id <- as.integer(data[all.vars(object$control$formGroup)][,1])
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))

  # survival
  tmp <- object$data[c(all.vars(object$control$formGroup),all.vars(object$control$formSurv))]
  tmp <- unique(tmp)
  Time <- tmp[all.vars(object$control$formSurv)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
  event <- tmp[all.vars(object$control$formSurv)][, 2]   # vector of event indicator (delta)
  mfZ <- stats::model.frame(object$control$formSurv, data = tmp)
  Z <- stats::model.matrix(object$control$formSurv, mfZ)

  # data management for shared latent structure
  if(object$control$param=="value"){
    data.id <- object$data[!duplicated(id), ]
    data.id[[object$control$timeVar]] <- Time
    mfX.id <- stats::model.frame(object$control$formFixed, data = data.id)
    Xtime <- stats::model.matrix(object$control$formFixed, mfX.id)
    mfU.id <- stats::model.frame(object$control$formRandom, data = data.id)
    Utime <- stats::model.matrix(object$control$formRandom, mfU.id)
    # approxitmation of the intergral via the Gaussian quadrature (Gauss Kronrod rule)
    sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
            0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
            -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
            0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
    wk <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
            0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
            0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
            0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
    K <- length(sk)
    P <- Time/2
    st <- outer(P, sk + 1)
    id.GK <- rep(seq_along(Time), each = K)
    data.id2 <- data.id[id.GK, ]
    data.id2[[object$control$timeVar]] <- c(t(st))
    mfX <- stats::model.frame(object$control$formFixed, data = data.id2)
    mfU <- stats::model.frame(object$control$formRandom, data = data.id2)
    Xs <- stats::model.matrix(object$control$formFixed, mfX)
    Us <- stats::model.matrix(object$control$formRandom, mfU)
  }

  # get parameter etimations
  beta <- object$mean$beta
  Cov_b <- object$median$covariance.b
  c1 <- (1-2*object$control$tau)/(object$control$tau*(1-object$control$tau))
  c2 <- 2/(object$control$tau*(1-object$control$tau))
  sigma <- object$mean$sigma
  if(object$control$survMod=="weibull"){
    alpha <- object$mean$alpha
    shape <- object$median$shape
    alpha.assoc <- object$mean$alpha.assoc
  }
  # get postrior median because distribution of this auxilary variable has skewed distribution
  w <- object$out_jagsUI$q50$va1
  Diag_Cov_y.b_w <- c2 * sigma * w

  # initialisation
  b <- array(NA, dim = c(I, M, ncol(U)))
  haz <- surv <- matrix(NA, nrow = I, ncol = M)
  if(object$control$survMod=="weibull")
    etaBaseline <- Z%*%alpha
  likelihood <- rep(NA, I)
  Xbeta <- X%*%beta
  Cov_b.tU <- Cov_b %*% t(U)
  if(object$control$param=="value"){
    SurvLong <- rep(NA, K)
    Xtimebeta <- Xtime%*%beta
  }
  # initialization of the progress bar
  if (verbose == TRUE) {
    pb <- utils::txtProgressBar(min = 0,
                                max = I,
                                initial = 0,
                                char = "*",
                                style = 3)
  }

  for (i in 1:I){
    # compute the individual contribution for the deviance

    # indice of observations for individual i
    indice_i <- offset[i]:(offset[i+1]-1)
    if(length(indice_i)==1){
      DD <- Diag_Cov_y.b_w[indice_i]
    }else{
      DD <- diag(Diag_Cov_y.b_w[indice_i])
    }
    ## draw sample of individual random effects b given y and w
    # Cov_y.w <- diag(Diag_Cov_y.b_w[indice_i]) + U[indice_i,] %*% Cov_b %*% t(U[indice_i, ])
    # mu_b.y_w <- Cov_b %*% t(U[indice_i, ]) %*% Inv_Cov_y.w %*% (y[indice_i] - Xbeta[indice_i] + c1 * w[indice_i])
    # Cov_b.y_w <- Cov_b - Cov_b %*% t(U[indice_i, ]) %*% Inv_Cov_y.w %*% U[indice_i, ] %*% Cov_b
    Cov_y.w <- DD + U[indice_i,] %*% Cov_b.tU[, indice_i]
    Inv_Cov_y.w <- solve(Cov_y.w)
    mu_b.y_w <- Cov_b.tU[, indice_i] %*% Inv_Cov_y.w %*% (y[indice_i] - Xbeta[indice_i] + c1 * w[indice_i])
    Cov_b.y_w <- Cov_b - Cov_b.tU[, indice_i] %*% Inv_Cov_y.w %*% U[indice_i, ] %*% Cov_b
    b[i, , ] <- MASS::mvrnorm(M, mu_b.y_w, Cov_b.y_w)

    # integration numeric for survival function
    if(object$control$survMod=="weibull")
      h0 <- shape * Time[i]^(shape-1) # without scale parameter
    for (m in 1:M) {
      if(object$control$survMod=="weibull"){
        if(object$control$param=="value"){
          # haz[i, m] <- h0 * exp( etaBaseline[i] + alpha.assoc * (Xtime[i, ]%*%beta +  Utime[i, ]%*%b[i, m, ]) )
          haz[i, m] <- h0 * exp( etaBaseline[i] + alpha.assoc * (Xtimebeta[i] +  Utime[i, ]%*%b[i, m, ]) )
          for (k in 1:K) {
            SurvLong[k] <- wk[k] * shape * st[i, k]^(shape-1) * exp(alpha.assoc * (Xs[K*(i-1)+k, ]%*%beta + Us[K*(i-1)+k, ]%*%b[i, m, ]))
          }
          surv[i, m] <- exp(-exp(etaBaseline[i]) * P[i] * sum(SurvLong[]))
        }
        if(object$control$param=="sharedRE"){
          tmp <- etaBaseline[i] + alpha.assoc%*%b[i, m, ]
          haz[i, m] <- h0 * exp( tmp )
          surv[i, m] <- exp( -exp(tmp) * Time[i]^shape )
        }
      }
    }
    likelihood[i] <- sum(surv[i, ]*haz[i, ]^event[i])/M

    #  Give method's progress
    if (verbose == TRUE) {
      utils::setTxtProgressBar(pb,i)
    }
  }

  if (verbose == TRUE)
    message("\n")

  # output management
  sims.list <- list(b = b,
                    hazard = haz,
                    survival = surv)

  deviance <- -2*sum(log(likelihood))
  deviance

  out <- list(deviance = deviance,
              likelihood = likelihood,
              sims.list = sims.list,
              control = list(M = M,
                             conditional = conditional,
                             verbose = verbose))
  out
}
