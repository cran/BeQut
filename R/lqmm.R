#' \code{lqmm} fits linear quantile mixed model
#'
#' Function using 'JAGS' software to estimate the linear quantile mixed model assuming asymmetric Laplace
#' distribution for residual error.
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without response variable
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param data dataset of observed variables
#' @param tau the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau}, otherwise the execution is stopped.
#' @param RE_ind Boolean denoting if the random effects are assumed independent ; default is \code{FALSE}
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param precision variance by default for vague prior distribution
#' @param save_jagsUI If \code{TRUE} (by default), the output of \code{jagsUI} package is returned by the function. Warning, if \code{TRUE}, the output can be large.
#' @param parallel see \code{jagsUI::jags()} function
#'
#'
#' @return A \code{Blqmm} object is a list with the following elements:
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
#' @import lqmm jagsUI
#'
#' @references Marco Geraci and Matteo Bottai (2014).
#' \emph{Linear quantile mixed models}.
#' Statistics and Computing, 24(3):461-479. doi: 10.1007/s11222-013-9381-9.
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' #---- Use dataLong dataset
#' data(dataLong)
#'
#' #---- Fit regression model for the first quartile
#' lqmm_075 <- lqmm(formFixed = y ~ visit,
#'                  formRandom = ~ visit,
#'                  formGroup = ~ ID,
#'                  data = dataLong,
#'                  tau = 0.75,
#'                  n.iter = 10000,
#'                  n.burnin = 1000)
#'
#' #---- Get the posterior means
#' lqmm_075$mean
#'
#' #---- Visualize the trace for beta parameters
#' jagsUI::traceplot(lqmm_075$out_jagsUI, parameters = "beta")
#'
#' #---- Summary of output
#' summary(lqmm_075)
#' }
#'
lqmm <- function(formFixed,
                 formRandom,
                 formGroup,
                 data,
                 tau,
                 RE_ind = FALSE,
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 5000,
                 n.thin = 1,
                 n.adapt = NULL,
                 precision = 10,
                 save_jagsUI = TRUE,
                 parallel = FALSE){


  #-- data management
  data_long <- data[unique(c(all.vars(formGroup),all.vars(formFixed),all.vars(formRandom)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- stats::model.frame(formFixed, data = data_long)
  X <- stats::model.matrix(formFixed, mfX)
  mfU <- stats::model.frame(formRandom, data = data_long)
  U <- stats::model.matrix(formRandom, mfU)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data_long)))
    data_long <- cbind(data_long, id = id)
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  # use lqmm function to initiated values
  message("Initiation of parameter values using lqmm package.")
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
                        priorB.Sigma2 = 1/precision)
                   )
    initial.values$prec.Sigma2 <- 1/VarCorr(tmp_model)
  }else{
    jags.data <- c(jags.data,
                   list(priorR.Sigma2 = diag(rep(1/precision, ncol(U))),
                        priorK.Sigma2 = ncol(U),
                        mu0 = rep(0, ncol(U))
                        )
                   )
    initial.values$prec.Sigma2 <- diag(1/VarCorr(tmp_model))
  }

  #---- write jags model in txt from R function
  if(RE_ind){
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
  }#end of i loop
  # priors for parameters
  for(rr in 1:ncU){
    prec.Sigma2[rr] ~ dgamma(priorA.Sigma2, priorB.Sigma2)
    covariance.b[rr] <- 1/prec.Sigma2[rr]
  }
  beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  sigma ~ dgamma(priorA.sigma, priorB.sigma)
}"
  }else{
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
  }#end of i loop
  # priors for parameters
  prec.Sigma2[1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
  covariance.b <- inverse(prec.Sigma2[, ])
  beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  sigma ~ dgamma(priorA.sigma, priorB.sigma)
}"
  }

  # regression from X
  rplc <- paste(paste("beta[", 1:jags.data$ncX, "] * X[j, ", 1:jags.data$ncX, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX], X[j, 1:ncX])", rplc, jags_code, fixed = TRUE)
  # regression on random effects b
  rplc <- paste(paste("b[i, ", 1:jags.data$ncU, "] * U[j, ", 1:jags.data$ncU, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(b[i, 1:ncU], U[j, 1:ncU])", rplc, jags_code, fixed = TRUE)

  # initialisation values
  if(n.chains==3)
    inits <- list(initial.values,
                  initial.values,
                  initial.values)
  if(n.chains==2)
    inits <- list(initial.values,
                  initial.values)
  if(n.chains==1)
    inits <- initial.values


  # parameters to save in the sampling step
  parms_to_save <- c("beta", "sigma", "b", "covariance.b")

  #---- use JAGS sampler via jagsUI
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

  #---- output
  out <- list(data = data)
  out$control <- list(formFixed = formFixed,
                      formRandom = formRandom,
                      formGroup = formGroup,
                      tau = tau,
                      call_function = "lqmm",
                      n.chains = n.chains,
                      parallel = parallel,
                      n.adapt = n.adapt,
                      n.iter = n.iter,
                      n.burnin = n.burnin,
                      n.thin = n.thin,
                      RE_ind = RE_ind,
                      I = I)

  #- other outputs

  # sims.list output
  out$sims.list <- out_jags$sims.list
  out$sims.list$b <- NULL

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

  # credible intervalles
  out$CIs$beta <- cbind(as.vector(t(out_jags$q2.5$beta)),
                        as.vector(t(out_jags$q97.5$beta)))
  rownames(out$CIs$beta) <- colnames(X)
  colnames(out$CIs$beta) <- c("2.5%", "97.5%")

  out$CIs$sigma <- c(out_jags$q2.5$sigma,
                     out_jags$q97.5$sigma)
  names(out$CIs$sigma) <- c("2.5%", "97.5%")

  # only for diagonal elements of covariance matrix of random effects
  if(RE_ind){
    out$CIs$covariance.b <- cbind(as.vector(t(out_jags$q2.5$covariance.b)),
                                  as.vector(t(out_jags$q97.5$covariance.b)))
    rownames(out$CIs$covariance.b) <- colnames(X)
    colnames(out$CIs$covariance.b) <- c("2.5%", "97.5%")
  }else{
    out$CIs$covariance.b <- cbind(as.vector(diag(out_jags$q2.5$covariance.b)),
                                 as.vector(diag(out_jags$q97.5$covariance.b)))
    rownames(out$CIs$covariance.b) <- colnames(U)
    colnames(out$CIs$covariance.b) <- c("2.5%", "97.5%")
  }

  # save jags output if requires
  if(save_jagsUI)
    out$out_jagsUI <- out_jags

  #---- End of the function defining the class and retruning the output
  class(out) <- "Blqmm"
  out

  }
