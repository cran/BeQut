#' \code{qrjm} fits quantile regression joint model
#'
#' Function using 'JAGS' software via \code{jagsUI} package to estimate the quantile regression joint model assuming asymmetric Laplace distribution for residual error.
#' Joint modeling concerns longitudinal data and time-to-event
#'
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without response variable
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param formSurv survival formula as formula in survival package for latency submodel
#' @param survMod specifying the baseline risk function for Cox proportional hazard model (only "weibull" is available until now)
#' @param param shared association including in joint modeling: the classical shared random effects or the current value denoting by "sharedRE" (default) or "value", respectively.
#' @param timeVar string specify the names of time variable (time of repeated measurements)
#' @param data dataset of observed variables
#' @param tau the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau}, otherwise the execution is stopped.
#' @param RE_ind Boolean denoting if the random effects are assumed independent ; default is \code{FALSE}
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is 5000
#' @param precision variance by default for vague prior distribution
#' @param C value used in the zero trick; default is 1000.
#' @param save_va If \code{TRUE} (is \code{FALSE} by default), the draws of auxiliary variable W is returned by the function
#' @param save_jagsUI If \code{TRUE} (by default), the output of \code{jagsUI} package is returned by the function
#' @param parallel see \code{jagsUI::jags()} function
#'
#'
#' @return A \code{Bqrjm} object is a list with the following elements:
#'  \describe{
#'   \item{\code{mean}}{list of posterior mean for each parameter}
#'   \item{\code{median}}{list of posterior median for each parameter}
#'   \item{\code{modes}}{list of posterior mode for each parameter}
#'   \item{\code{StErr}}{list of standard error for each parameter}
#'   \item{\code{StDev}}{list of standard deviation for each parameter}
#'   \item{\code{ICs}}{list of the credibility interval at 0.95 for each parameters excepted for covariance parameters in covariance matrix of random effects. Otherwise, use save_jagsUI=TRUE to have the associated quantiles.}
#'   \item{\code{data}}{data included in argument}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters and random effects}
#'   \item{\code{control}}{list of arguments giving details about the estimation}
#'   \item{\code{random_effect}}{list for each quantile including both posterior mean and posterior standard deviation of subject-specific random effects}
#'   \item{\code{out_jagsUI}}{only if \code{save_jagsUI=TRUE} in argument: list including posterior mean, median, quantiles (2.5%, 25%, 50%, 75%, 97.5%), standart deviation for each parameter and each random effect.
#'   Moreover, this list also returns the MCMC draws, the Gelman and Rubin diagnostics (see output of jagsUI objects)}
#'  }
#'
#' @author Antoine Barbieri
#'
#' @import jagsUI lqmm survival
#'
#' @export
#'
#' @references Ming Yang, Sheng Luo, and Stacia DeSantis (2019).
#' \emph{Bayesian quantile regression joint models: Inference and dynamic predictions}.
#' Statistical Methods in Medical Research, 28(8):2524-2537. doi: 10.1177/0962280218784757.
#'
#' @examples
#'
#' \donttest{
#' #---- load data
#' data(dataLong)
#'
#' #---- Fit quantile regression joint model for the first quartile
#' qrjm_75 <- qrjm(formFixed = y ~ visit,
#'                formRandom = ~ visit,
#'                formGroup = ~ ID,
#'                formSurv = Surv(time, event) ~ X1 + X2,
#'                survMod = "weibull",
#'                param = "value",
#'                timeVar= "visit",
#'                data = dataLong,
#'                tau = 0.75)
#'
#' #---- Visualize the trace for beta parameters
#' jagsUI::traceplot(qrjm_75$out_jagsUI, parameters = "beta")
#'
#' #---- Get the estimated coefficients: posterior means
#' qrjm_75$mean
#'
#' #---- Summary of output
#' summary(qrjm_75)
#' }
#'
qrjm <- function(formFixed,
                 formRandom,
                 formGroup,
                 formSurv,
                 survMod = "weibull",
                 param = "value",
                 timeVar,
                 data,
                 tau,
                 RE_ind = FALSE,
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 5000,
                 n.thin = 1,
                 n.adapt = 5000,
                 precision = 10,
                 C = 1000,
                 save_jagsUI = TRUE,
                 save_va = FALSE,
                 parallel = FALSE){

  # #
  # #   -- To do
  # #   verify with value.IG
  # #   add a stopping convergence criteria
  # #   initialize the values of parameter chains and to fix the intercept to test the convergence of beta's parameter.
  # #
  #

  #-- data management

  # control
  lag = 0

  #--- longitudinal part
  data_long <- data[unique(c(all.vars(formGroup),all.vars(formFixed),all.vars(formRandom)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- stats::model.frame(formFixed, data = data_long)
  X <- stats::model.matrix(formFixed, mfX)
  mfU <- stats::model.frame(formRandom, data = data_long)
  U <- stats::model.matrix(formRandom, mfU)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  if(!("id" %in% colnames(data_long)))
    data_long <- cbind(data_long, id = id)
  # use lqmm function to initiated values
  message("Initialisation of longitudinal parameter values using 'lqmm' package.")
  tmp_model <- lqmm::lqmm(fixed = formFixed,
                          random = formRandom,
                          group = id,
                          tau = tau,
                          data = data_long)
  # prior beta parameters
  priorMean.beta <- coef.lqmm(tmp_model)
  priorTau.beta <- diag(rep(1/10,length(priorMean.beta)))

  bis <- as.matrix(lqmm::ranef(tmp_model))
  bis[abs(bis)<.0001] <- 0
  initial.values <- list(b = bis,
                         beta = priorMean.beta,
                         sigma = tmp_model$scale)

  # list of data jags
  jags.data <- list(y = y,
                    X = X,
                    U = U,
                    tau = tau,
                    ncX = ncol(X),
                    ncU = ncol(U),
                    I = I,
                    offset = offset,
                    priorMean.beta = priorMean.beta,
                    priorTau.beta = priorTau.beta,
                    priorA.sigma = 1/precision,
                    priorB.sigma = 1/precision
                    )

  if(jags.data$ncU==1)
    RE_ind <- TRUE
  if(RE_ind){
    jags.data <- c(jags.data,
                   list(priorA.Sigma2 = 1/precision,
                        priorB.Sigma2 = 1/precision
                   )
    )
    initial.values$prec.Sigma2 <- 1/VarCorr(tmp_model)
  }else{
    jags.data <- c(jags.data,
                   list(priorR.Sigma2 = diag(rep(1/precision, ncol(U))),
                        priorK.Sigma2 = ncol(U),
                        mu0 = rep(0, ncol(U))
                   )
    )
    initial.values$prec.Sigma2 <- diag(1/lqmm::VarCorr(tmp_model))
  }

  #--- survival part
  tmp <- data[c(all.vars(formGroup),all.vars(formSurv))]
  tmp <- unique(tmp)
  Time <- tmp[all.vars(formSurv)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
  event <- tmp[all.vars(formSurv)][, 2]   # vector of event indicator (delta)
  nTime <- length(Time)                   # number of subject having Time
  zeros <- numeric(nTime)                 # for zero trick in Bayesian procedure
  # design matrice
  mfZ <- stats::model.frame(formSurv, data = tmp)
  Z <- stats::model.matrix(formSurv, mfZ)
  # use survival::coxph function to initiated values
  message("Initialisation of survival parameter values using 'survival' package.")
  tmp_model <- survival::coxph(formSurv,
                               data = tmp,
                               x = TRUE)
  # Complete the jags data
  priorMean.alpha <- c(0, tmp_model$coefficients)
  priorTau.alpha <- diag(c(1/precision, 1/(precision*diag(tmp_model$var))))
  jags.data <- c(jags.data,
                 list(C = C,
                      zeros = numeric(nTime),
                      Time = Time,
                      event = event,
                      Z = Z,
                      ncZ = ncol(Z),
                      priorMean.alpha = priorMean.alpha,
                      priorTau.alpha = priorTau.alpha,
                      priorTau.alphaA = 1/precision)
                 )
  # initialisation values of survival parameters
  initial.values$alpha <- c(0, tmp_model$coefficients)
  # if(survMod=="weibull")
  #   initial.values$shape <- 1
  if(param=="value")
    initial.values$alpha.assoc <- 0
  if(param=="sharedRE")
    initial.values$alpha.assoc <- rep(0, ncol(U))

  #--- shared current value case
  data.id <- data_long[!duplicated(id), ]
  if (!timeVar %in% names(data_long))
    stop("\n'timeVar' does not correspond to one of the columns in formulas")
  if (param %in% c("value")) {
    data.id[[timeVar]] <- pmax(Time - lag, 0)
    mfX.id <- stats::model.frame(formFixed, data = data.id)
    Xtime <- stats::model.matrix(formFixed, mfX.id)
    mfU.id <- stats::model.frame(formRandom, data = data.id)
    Utime <- stats::model.matrix(formRandom, mfU.id)
    # if (one.RE)
    #   Utime <- cbind(Utime, rep(0, nrow(Utime)))
    jags.data <- c(jags.data, list(Xtime = Xtime, Utime = Utime))

    #-- approxitmation of the intergral via the Gaussian quadrature (Gauss Kronrod rule)
    gaussKronrod <-
      function (k = 15) {
        sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
                0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
                -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
                0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
        wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                  0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                  0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                  0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
        wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
                 0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
        if (k == 7)
          list(sk = sk[1:7], wk = wk7)
        else
          list(sk = sk, wk = wk15)
      }

    wk <- gaussKronrod()$wk
    sk <- gaussKronrod()$sk
    K <- length(sk)
    P <- Time/2
    st <- outer(P, sk + 1)
    id.GK <- rep(seq_along(Time), each = K)
    data.id2 <- data.id[id.GK, ]
    data.id2[[timeVar]] <- c(t(st))
    mfX <- stats::model.frame(formFixed, data = data.id2)
    mfU <- stats::model.frame(formRandom, data = data.id2)
    Xs <- stats::model.matrix(formFixed, mfX)
    Us <- stats::model.matrix(formRandom, mfU)

    jags.data <- c(jags.data, list(K = K, P = P, st = st, wk = wk, Xs = Xs, Us = Us))
  }

  #---- Model for JAGS

  # for weibull baseline hazard function and current value as shared association
  if(param=="value"){
    if(RE_ind){
      # inverse gamma distribution when RE are considered as independent
      jags_code <- "model{
# constants
c1 <- (1-2*tau)/(tau*(1-tau))
c2 <- 2/(tau*(1-tau))
# likelihood
for (i in 1:I){
  # longitudinal part
  for(j in offset[i]:(offset[i+1]-1)){
    y[j] ~ dnorm(mu[j], prec[j])
    va1[j] ~ dexp(1/sigma)
    prec[j] <- 1/(sigma*c2*va1[j])
    mu[j] <- inprod(beta[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncU], U[j, 1:ncU]) + c1*va1[j]
  }#end of j loop
# random effects
  for(r in 1:ncU){
    b[i,r] ~ dnorm(0, prec.Sigma2[r])
  }
  # survival part
  etaBaseline[i] <- inprod(alpha[1: ncZ], Z[i, 1:ncZ])
  shareY[i] <- inprod(beta[1:ncX], Xtime[i, 1:ncX]) + inprod(b[i, 1:ncU], Utime[i, 1:ncU])
  log_h1[i] <- log(shape) + (shape - 1) * log(Time[i]) + etaBaseline[i] + alpha.assoc * shareY[i]
  for (k in 1:K) {
    shareY.s[i, k] <- inprod(beta[1:ncX], Xs[K * (i - 1) + k, 1:ncX]) + inprod(b[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])
    SurvLong[i, k] <- wk[k] * shape * pow(st[i, k], shape - 1) * exp(alpha.assoc * shareY.s[i, k])
  }
  log_S1[i] <- (-exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ]))
  logL[i] <- event[i]*log_h1[i] + log_S1[i]
  mlogL[i] <- -logL[i] + C
  zeros[i] ~ dpois(mlogL[i])
}#end of i loop
# priors for longitudinal parameters
for(rr in 1:ncU){
  prec.Sigma2[rr] ~ dgamma(priorA.Sigma2, priorB.Sigma2)
  covariance.b[rr] <- 1/prec.Sigma2[rr]
}
beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
sigma ~ dgamma(priorA.sigma, priorB.sigma)
# priors for survival parameters
alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
shape ~ dgamma(priorA.shape, priorB.shape)
alpha.assoc ~ dnorm(0, priorTau.alphaA)
}"
    }else{
      # wishart distribution when RE are considered as dependent
      jags_code <- "model{
# constants
c1 <- (1-2*tau)/(tau*(1-tau))
c2 <- 2/(tau*(1-tau))
# likelihood
for (i in 1:I){
  # longitudinal part
  for(j in offset[i]:(offset[i+1]-1)){
    y[j] ~ dnorm(mu[j], prec[j])
    va1[j] ~ dexp(1/sigma)
    prec[j] <- 1/(sigma*c2*va1[j])
    mu[j] <- inprod(beta[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncU], U[j, 1:ncU]) + c1*va1[j]
  }#end of j loop
  # random effects
  b[i, 1:ncU] ~ dmnorm(mu0[], prec.Sigma2[, ])
  # survival part
  etaBaseline[i] <- inprod(alpha[1: ncZ], Z[i, 1:ncZ])
  shareY[i] <- inprod(beta[1:ncX], Xtime[i, 1:ncX]) + inprod(b[i, 1:ncU], Utime[i, 1:ncU])
  log_h1[i] <- log(shape) + (shape - 1) * log(Time[i]) + etaBaseline[i] + alpha.assoc * shareY[i]
  for (k in 1:K) {
    shareY.s[i, k] <- inprod(beta[1:ncX], Xs[K * (i - 1) + k, 1:ncX]) + inprod(b[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])
    SurvLong[i, k] <- wk[k] * shape * pow(st[i, k], shape - 1) * exp(alpha.assoc * shareY.s[i, k])
  }
  log_S1[i] <- (-exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ]))
  logL[i] <- event[i]*log_h1[i] + log_S1[i]
  mlogL[i] <- -logL[i] + C
  zeros[i] ~ dpois(mlogL[i])
}#end of i loop
# priors for longitudinal parameters
prec.Sigma2[1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
covariance.b <- inverse(prec.Sigma2[, ])
beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
sigma ~ dgamma(priorA.sigma, priorB.sigma)
# priors for survival parameters
alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
shape ~ dgamma(priorA.shape, priorB.shape)
alpha.assoc ~ dnorm(0, priorTau.alphaA)
}"
    }
  # replace inprod
  # regression from X
  rplc <- paste(paste("beta[", 1:jags.data$ncX, "] * X[j, ", 1:jags.data$ncX, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX], X[j, 1:ncX])", rplc, jags_code, fixed = TRUE)
  # regression for random effects b
  rplc <- paste(paste("b[i, ", 1:jags.data$ncU, "] * U[j, ", 1:jags.data$ncU, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(b[i, 1:ncU], U[j, 1:ncU])", rplc, jags_code, fixed = TRUE)
  # regression on survival part for time-independant covariates
  rplc <- paste(paste("alpha[", 1:jags.data$ncZ, "] * Z[i, ", 1:jags.data$ncZ, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(alpha[1: ncZ], Z[i, 1:ncZ])", rplc, jags_code, fixed = TRUE)
  # regression on survival part for share association in hazard function
  rplc <- paste(paste("beta[", 1:jags.data$ncX, "] * Xtime[i, ", 1:jags.data$ncX, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX], Xtime[i, 1:ncX])", rplc, jags_code, fixed = TRUE)
  rplc <- paste(paste("b[i, ", 1:jags.data$ncU, "] * Utime[i, ", 1:jags.data$ncU, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(b[i, 1:ncU], Utime[i, 1:ncU])", rplc, jags_code, fixed = TRUE)
  # regression on survival part for share association in survival function with GK intergration
  rplc <- paste(paste("beta[", 1:jags.data$ncX, "] * Xs[K * (i - 1) + k, ", 1:jags.data$ncX, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX], Xs[K * (i - 1) + k, 1:ncX])", rplc, jags_code, fixed = TRUE)
  rplc <- paste(paste("b[i,", 1:jags.data$ncU, "] * Us[K * (i - 1) + k, ", 1:jags.data$ncU, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(b[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])", rplc, jags_code, fixed = TRUE)

  }

  # for weibull baseline hazard function and shared random effects as association
  if(param=="sharedRE"){
    if(RE_ind){
      # inverse gamma distribution when RE are considered as independent
      jags_code <- "model{
# constants
  c1 <- (1-2*tau)/(tau*(1-tau))
  c2 <- 2/(tau*(1-tau))
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      va1[j] ~ dexp(1/sigma)
      prec[j] <- 1/(sigma*c2*va1[j])
      mu[j] <- inprod(beta[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncU], U[j, 1:ncU]) + c1*va1[j]
    }#end of j loop
    # random effects
    for(r in 1:ncU){
      b[i,r] ~ dnorm(0, prec.Sigma2[r])
    }
    # survival part
    etaBaseline[i] <- inprod(alpha[1: ncZ], Z[i, 1:ncZ]) + inprod(alpha.assoc[1:ncU], b[i, 1:ncU])
    log_S1[i] <- -exp(etaBaseline[i]) * pow(Time[i], shape)
    log_h1[i] <- log(shape) + (shape-1)*log(Time[i]) + etaBaseline[i]
    logL[i] <- event[i]*log_h1[i] + log_S1[i]
    mlogL[i] <- -logL[i] + C
    zeros[i] ~ dpois(mlogL[i])
  }#end of i loop
  # priors for longitudinal parameters
  for(rr in 1:ncU){
    prec.Sigma2[rr] ~ dgamma(priorA.Sigma2, priorB.Sigma2)
    covariance.b[rr] <- 1/prec.Sigma2[rr]
  }
  beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  sigma ~ dgamma(priorA.sigma, priorB.sigma)
  # priors for survival parameters
  alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
  shape ~ dgamma(priorA.shape, priorB.shape)
  for(rr in 1:ncU){
    alpha.assoc[rr] ~ dnorm(0, priorTau.alphaA)
  }
}"
    }else{
    # wishart distribution when RE are considered as dependent
    jags_code <- "model{
# constants
  c1 <- (1-2*tau)/(tau*(1-tau))
  c2 <- 2/(tau*(1-tau))
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      va1[j] ~ dexp(1/sigma)
      prec[j] <- 1/(sigma*c2*va1[j])
      mu[j] <- inprod(beta[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncU], U[j, 1:ncU]) + c1*va1[j]
    }#end of j loop
    # random effects
    b[i, 1:ncU] ~ dmnorm(mu0[], prec.Sigma2[, ])
    # survival part
    etaBaseline[i] <- inprod(alpha[1: ncZ], Z[i, 1:ncZ]) + inprod(alpha.assoc[1:ncU], b[i, 1:ncU])
    log_S1[i] <- -exp(etaBaseline[i]) * pow(Time[i], shape)
    log_h1[i] <- log(shape) + (shape-1)*log(Time[i]) + etaBaseline[i]
    logL[i] <- event[i]*log_h1[i] + log_S1[i]
    mlogL[i] <- -logL[i] + C
    zeros[i] ~ dpois(mlogL[i])
  }#end of i loop
  # priors for longitudinal parameters
  prec.Sigma2[1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
  covariance.b <- inverse(prec.Sigma2[, ])
  beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  sigma ~ dgamma(priorA.sigma, priorB.sigma)
  # priors for survival parameters
  alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
  shape ~ dgamma(priorA.shape, priorB.shape)
  for(rr in 1:ncU){
    alpha.assoc[rr] ~ dnorm(0, priorTau.alphaA)
  }
}"
  }

  # replace inprod
  # regression from X
  rplc <- paste(paste("beta[", 1:jags.data$ncX, "] * X[j, ", 1:jags.data$ncX, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX], X[j, 1:ncX])", rplc, jags_code, fixed = TRUE)
  # regression for random effects b
  rplc <- paste(paste("b[i, ", 1:jags.data$ncU, "] * U[j, ", 1:jags.data$ncU, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(b[i, 1:ncU], U[j, 1:ncU])", rplc, jags_code, fixed = TRUE)
  # regression on survival part for time-independant covariates
  rplc <- paste(paste("alpha[", 1:jags.data$ncZ, "] * Z[i, ", 1:jags.data$ncZ, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(alpha[1: ncZ], Z[i, 1:ncZ])", rplc, jags_code, fixed = TRUE)
  # regression on survival part for shared association in both hazard survival function
  rplc <- paste(paste("alpha.assoc[", 1:jags.data$ncU, "] * b[i,  ", 1:jags.data$ncU, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(alpha.assoc[1:ncU], b[i, 1:ncU])", rplc, jags_code, fixed = TRUE)

  }

  # initialisation values
  if (n.chains == 3)
    inits <- list(initial.values,
                  initial.values,
                  initial.values)
  if (n.chains == 2)
    inits <- list(initial.values,
                  initial.values)
  if (n.chains == 1)
    inits <- initial.values

  # parameters to save in the sampling step
  parms_to_save <- c("alpha", "alpha.assoc", "beta", "sigma", "b", "covariance.b")
  if (save_va)
    parms_to_save <- c(parms_to_save, "va1")

  # complement given survMod
  if (survMod == "weibull") {
    jags.data <-
      c(jags.data,
        list(
          priorA.shape = 1 / precision,
          priorB.shape = 1 / precision
       ))
    parms_to_save <- c(parms_to_save, "shape")
  }

  # using jagsUI
  out_jags = jagsUI::jags(data = jags.data,
                          parameters.to.save = parms_to_save,
                          model.file = textConnection(jags_code),
                          inits = inits,
                          n.chains = n.chains,
                          parallel = parallel,
                          n.adapt = n.adapt,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          DIC = F)

  #---- output building

  #-- MCMClist management

  #- arguments
  out <- list(data = data)
  out$control <- list(formFixed = formFixed,
                      formRandom = formRandom,
                      formGroup = formGroup,
                      formSurv = formSurv,
                      timeVar = timeVar,
                      tau = tau,
                      call_function = "qrjm",
                      I = I,
                      C = C,
                      param = param,
                      survMod = survMod,
                      n.chains = n.chains,
                      parallel = parallel,
                      n.adapt = n.adapt,
                      n.iter = n.iter,
                      n.burnin = n.burnin,
                      n.thin = n.thin,
                      RE_ind = RE_ind,
                      event = event,
                      Time = Time)

  #- other outputs

  # sims.list output
  out$sims.list <- out_jags$sims.list
  out$sims.list$b <- NULL
  if(save_va)
    out$sims.list$va1 <- NULL

  # random effect output
  out$random_effect <- list(postMeans = out_jags$mean$b,
                            postSd = out_jags$sd$b)
  colnames(out$random_effect$postMeans) <- colnames(U)
  colnames(out$random_effect$postSd) <- colnames(U)

  # median : Posterior median of parameters
  out$median <- out_jags$q50
  out$median$b <- NULL

  # mean : Posterior mean of parameters
  out$mean <- out_jags$mean
  out$mean$b <- NULL

  # modes of parameters
  out$modes <- lapply(out$sims.list, function(x) {
    m <- function(x) {
      d <- stats::density(x, bw = "nrd", adjust = 3, n = 1000)
      d$x[which.max(d$y)]
    }
    if (is.matrix(x))
      as.array(apply(x, 2, m))
    else{
      if(is.array(x))
        apply(x, c(2,3), m)
      else m(x)
    }
  })

  # standard error of parameters
  out$StErr <- lapply(out$sims.list, function(x) {
    f <- function(x) {
      acf.x <- drop(stats::acf(x, lag.max = 0.4 * length(x), plot = FALSE)$acf)[-1]
      acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
      ess <- length(x)/(1 + 2 * sum(acf.x))
      sqrt(stats::var(x)/ess)
    }
    if (is.matrix(x))
      as.array(apply(x, 2, f))
    else{
      if(is.array(x))
        apply(x, c(2,3), f)
      else f(x)
    }
  })

  # standard deviation of parameters
  out$StDev <- out_jags$sd
  out$StDev$b <- NULL

  # Rhat : Gelman & Rubin diagnostic
  out$Rhat <- out_jags$Rhat
  out$Rhat$b <- NULL

  # clean regarding auxilary variable (information avalaible in jagsUI output)
  if(save_va)
    out$median$va1 <- out$mean$va1 <- out$StDev$va1 <- out$Rhat$va1 <- NULL

  # names
  names(out$mean$beta) <-
    names(out$median$beta) <-
    names(out$modes$beta) <-
    names(out$StErr$beta) <-
    names(out$Rhat$beta) <-
    names(out$StDev$beta) <- colnames(X)

  if(RE_ind){
    names(out$mean$covariance.b) <-
      names(out$median$covariance.b) <-
      names(out$modes$covariance.b) <-
      names(out$StErr$covariance.b) <-
      names(out$Rhat$covariance.b) <-
      names(out$StDev$covariance.b) <- colnames(U)
  }else{
    colnames(out$mean$covariance.b) <-
      rownames(out$mean$covariance.b) <-
      colnames(out$median$covariance.b) <-
      rownames(out$median$covariance.b) <-
      colnames(out$modes$covariance.b) <-
      rownames(out$modes$covariance.b) <-
      colnames(out$Rhat$covariance.b) <-
      rownames(out$Rhat$covariance.b) <-
      colnames(out$StErr$covariance.b) <-
      rownames(out$StErr$covariance.b) <-
      colnames(out$StDev$covariance.b) <-
      rownames(out$StDev$covariance.b) <- colnames(U)
  }

  names(out$mean$alpha) <-
    names(out$median$alpha) <-
    names(out$modes$alpha) <-
    names(out$Rhat$alpha) <-
    names(out$StErr$alpha) <-
    names(out$StDev$alpha) <- colnames(Z)

  # credible intervalles
  out$CIs$beta <- cbind(as.vector(t(out_jags$q2.5$beta)),
                        as.vector(t(out_jags$q97.5$beta)))
  rownames(out$CIs$beta) <- colnames(X)
  colnames(out$CIs$beta) <- c("2.5%", "97.5%")

  out$CIs$sigma <- c(out_jags$q2.5$sigma,
                     out_jags$q97.5$sigma)
  names(out$CIs$sigma) <- c("2.5%", "97.5%")

  out$CIs$shape <- c(out_jags$q2.5$shape,
                     out_jags$q97.5$shape)
  names(out$CIs$shape) <- c("2.5%", "97.5%")

  if(param %in% c("value")){
    out$CIs$alpha.assoc <- c(out_jags$q2.5$alpha.assoc,
                             out_jags$q97.5$alpha.assoc)
    names(out$CIs$alpha.assoc) <- c("2.5%", "97.5%")
  }
  if(param %in% c("sharedRE")){
    if(jags.data$ncU>1){
      out$CIs$alpha.assoc <- cbind(as.vector(t(out_jags$q2.5$alpha.assoc)),
                                   as.vector(t(out_jags$q97.5$alpha.assoc)))
      rownames(out$CIs$alpha.assoc) <- colnames(U)
      colnames(out$CIs$alpha.assoc) <- c("2.5%", "97.5%")
    }else{
      out$CIs$alpha.assoc <- c(out_jags$q2.5$alpha.assoc,
                               out_jags$q97.5$alpha.assoc)
      names(out$CIs$alpha.assoc) <- c("2.5%", "97.5%")
    }
  }

  out$CIs$alpha <- cbind(as.vector(t(out_jags$q2.5$alpha)),
                         as.vector(t(out_jags$q97.5$alpha)))
  rownames(out$CIs$alpha) <- colnames(Z)
  colnames(out$CIs$alpha) <- c("2.5%", "97.5%")

  # only for diagonal elements of covariance matrix of random effects
  if(RE_ind){
    out$CIs$variances.b <- cbind(as.vector(out_jags$q2.5$covariance.b),
                                 as.vector(out_jags$q97.5$covariance.b))
  }else{
    out$CIs$variances.b <- cbind(as.vector(diag(out_jags$q2.5$covariance.b)),
                                 as.vector(diag(out_jags$q97.5$covariance.b)))
  }
  rownames(out$CIs$variances.b) <- colnames(U)
  colnames(out$CIs$variances.b) <- c("2.5%", "97.5%")

  # save jags output if requires
  if(save_jagsUI)
    out$out_jagsUI <- out_jags

  #---- End of the function defining the class and returning the output
  class(out) <- "Bqrjm"
  out

}
