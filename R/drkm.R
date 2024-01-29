#' Estimate the survival curves in two treatment groups
#'
#' Doubly robust estimation of the survival curves in two exposure groups
#' and their difference using the method of Sjolander and Vansteelandt (2017) .
#'
#' Following Bai et al. 2013 and Sjölander and Vansteelandt 2017, consider the
#' following estimating equation for \eqn{S_x(t) = p\{T(x) > t\}}:
#' \deqn{\sum_{i = 1}^n \left[S_x(t) - \frac{I_{X_i = x}
#' I_{\tilde{T}_i > t}}{\bar{g}(\boldsymbol{Z}_i, \boldsymbol{\alpha})G_c(t,
#' X_i, \boldsymbol{Z}_i)} - \frac{I_{X_i = x} - \bar{g}(\boldsymbol{Z}_i,
#' \boldsymbol{\alpha})}{\bar{g}(\boldsymbol{Z}_i, \boldsymbol{\alpha})} H(t,
#' \boldsymbol{Z}_i, X = x) - \right. \\ \left. \frac{I_{X_i = x} -
#' \bar{g}(\boldsymbol{Z}_i, \boldsymbol{\alpha})}{\bar{g}(\boldsymbol{Z}_i,
#' \boldsymbol{\alpha})} H(t, \boldsymbol{Z}_i, X = x) \int_0^t\frac{d\, M_c(u,
#' \boldsymbol{Z}_i, X_i, \tilde{T}_i, \Delta_i)}{G_c(u, X_i,
#' \boldsymbol{Z}_i)H(u, \boldsymbol{Z}_i, X_i)}\right] = 0, }
#' where \eqn{\bar{g}(\boldsymbol{Z}_i, \boldsymbol{\alpha}) = g(\boldsymbol{Z}_i,
#' \boldsymbol{\alpha})^x (1 - g(\boldsymbol{Z}_i,
#' \boldsymbol{\alpha}))^{(1-x)}}, \eqn{H(t, \boldsymbol{Z}, X)} is a model for
#' \eqn{p\{T > t | \boldsymbol{Z}, X\}}, and \eqn{M_c(t, \boldsymbol{Z}, X, \tilde{T},
#' \Delta)} is the martingale increment for the censoring distribution. The
#' above is an unbiased estimating equation for \eqn{S_x(t)} if either \eqn{H(t,
#' \boldsymbol{Z}, X)} is correctly specified for both censoring and
#' confounding, i.e., is a correctly specified model for \eqn{p\{T(x) > t |
#' \boldsymbol{Z}, X\}} or both \eqn{G_c(u, X, \boldsymbol{Z})} and
#' \eqn{g(\boldsymbol{Z}, \boldsymbol{\alpha})} are correctly specified for
#' censoring and confounding, respectively. To obtain estimates of the
#' difference in survival probabilities, one must specify models for the unknown
#' functions \eqn{g}, \eqn{G_c}, and \eqn{H}, get estimates of those, plug them
#' into the estimating equations, and solve for \eqn{S_x(t)} under \eqn{x \in
#' \{0, 1\}}. In this package, one can use semiparametric Cox models or
#' parametric survival models for the outcome and the censoring distributions,
#' and logistic regression for the propensity score model. Bai et al. 2013
#' provide an expression for a variance estimator that accounts for the
#' uncertainty due to the estimation of the propensity score \eqn{g} and the
#' censoring distribution \eqn{G_c}.
#'
#' @references Xiaofei Bai, Anastasios A Tsiatis, and Sean~M O'Brien.
#'   "Doubly-robust estimators of treatment-specific survival distributions in
#'   observational studies with stratified sampling." Biometrics,
#'   69(4):830--839, 2013.
#'
#'   Sjölander, Arvid, and Stijn Vansteelandt. "Doubly
#'   robust estimation of attributable fractions in survival analysis."
#'   Statistical methods in medical research 26(2):948-969, 2017.
#'
#' @param oformula The outcome formula, the left side should be a call to \link[survival]{Surv}
#' @param ofunc The model type for the outcome, one of "coxph" or "survreg"
#' @param oarg Arguments passed to \code{ofunc}
#' @param cformula The censoring model formula, the left side should be a call to \link[survival]{Surv}
#' @param cfunc The model type for the censoring model, either "coxph" or "survreg"
#' @param carg Arguments passed to \code{cfunc}
#' @param eformula The exposure model formula, which will be fit by logistic regression in a call to \link[stats]{glm}
#' @param earg Arguments passed to the glm for the exposure model
#' @param method Estimation method, either "IPW" or "DR"
#' @param times Vector of times at which to estimate the survival curves
#' @param rel.tol Convergence tolerance
#' @param jacobian.method Method of computing the jacobian, passed to \link[numDeriv]{jacobian}
#' @param se.type Method of calculating standard errors, either "none" for no standard errors, "sandwich", or "boot"
#' @param data Data frame in which to find the variables in the model formulas
#' @param weights Vector of case weights
#' @param subset Logical vector
#' @param R Number of bootstrap replicates, used if \code{se.type = "boot"}
#' @param parallel Parallel processing for bootstrap, see \link[boot]{boot}
#' @param cl Optional cluster if parallel processing, see \link[boot]{boot}
#' @param ncpus Number of cpus to use for parallel processing, see \link[boot]{boot}
#' @returns A list with the estimates survival probabilities in each group, their difference, and standard errors, if requested
#'
#' @examples
#' df <- rotterdam
#' df$time <- pmin(df$rtime, df$dtime) / 365.25
#' df$status <- ifelse(df$recur == 1 | df$death == 1, 1, 0)
#' df$censor <- 1 - df$status
#' drFit <-
#'   sjovan_survdiff(
#'     oformula = Surv(time, status) ~ chemo + year + age + meno +
#'       size + factor(grade) + nodes + pgr + er + hormon,
#'     ofunc = "survreg",
#'     cformula = Surv(time, censor) ~ chemo + year + age,
#'     cfunc = "survreg",
#'     eformula = chemo ~ year + age + meno + size +
#'       factor(grade) + nodes + pgr + er + hormon,
#'     method = "DR",
#'     times = c(2.5, 5, 7.5),
#'     se.type = "sandwich",
#'     data = df
#'   )
#' drFit
#'
#' @export


sjovan_survdiff <- function(oformula = NULL, ofunc = "coxph", oarg = list(),
                       cformula = NULL, cfunc = "coxph", carg = list(),
                       eformula = NULL, earg = list(),
                       method = "DR", times = NULL,
                       rel.tol = .Machine$double.eps^0.1, jacobian.method = "simple",
                       se.type = "none", data = NULL, weights = NULL, subset = NULL,
                       R = 50, parallel = "no", cl = NULL, ncpus = NULL){

  if(!is.null(subset)) data <- data[subset, ]
  data <- na.omit(data)
  n <- nrow(data)
  if(is.null(weights)) weights <- rep(1, n)

  #---check that no left truncation, except if method == "ML"
  if((method == "IPW" | method == "DR") &
     length(as.character(oformula[[2]])) == 4)
    stop("Left truncation only allowed if method == ML", call. = FALSE)

  if(se.type == "boot"){

    bb <- boot::boot(data = data, statistic = bootfun, R = R, parallel = parallel,
               cl = cl, ncpus = ncpus,
               oformula = oformula, ofunc = ofunc, oarg = oarg,
               cformula = cformula, cfunc = cfunc, carg = carg,
               eformula = eformula, earg = earg,
               method = method, times = times, rel.tol = rel.tol,
               ww = weights)

    est <- bb$t0
    se <-  apply(X = bb$t, MARGIN = 2, FUN = sd)
    out <- list(est.S1 = est[seq(1,3*length(times)-2,3)],
                se.S1 = se[seq(1,3*length(times)-2,3)],
                est.S0 = est[seq(2,3*length(times)-1,3)],
                se.S0 = se[seq(2,3*length(times)-1,3)],
                est.diff = est[seq(3,3*length(times),3)],
                se.diff = se[seq(3,3*length(times),3)])

  }

  if(se.type == "sandwich"){

    if((ofunc == "coxph" & (method == "ML" | method == "DR")) |
       (cfunc == "coxph" & (method == "IPW" | method == "DR")))
      stop("Sandwich estimator of variance not allowed for Cox models",
           call. = FALSE)

    fit <- afSurvival.fit(oformula = oformula, ofunc = ofunc, oarg = oarg,
                          cformula = cformula, cfunc = cfunc, carg = carg,
                          eformula = eformula, earg = earg,
                          method = method, times = times,
                          rel.tol = rel.tol, jacobian.method = jacobian.method,
                          se.fit = TRUE, data = data, weights = weights)
    est <- fit$est
    se <- fit$se
    out <- list(est.S1 = est[seq(1,3*length(times)-2,3)],
                se.S1 = se[seq(1,3*length(times)-2,3)],
                est.S0 = est[seq(2,3*length(times)-1,3)],
                se.S0 = se[seq(2,3*length(times)-1,3)],
                est.diff = est[seq(3,3*length(times),3)],
                se.diff = se[seq(3,3*length(times),3)])

  }

  if(se.type == "none"){

    fit <- afSurvival.fit(oformula = oformula, ofunc = ofunc, oarg = oarg,
                          cformula = cformula, cfunc = cfunc, carg = carg,
                          eformula = eformula, earg = earg,
                          method = method, times = times,
                          rel.tol = rel.tol, jacobian.method = jacobian.method,
                          se.fit = FALSE, data = data, weights = weights)
    est <- fit$est
    se <- fit$se
    out <- list(est.S1 = est[seq(1,3*length(times)-2,3)],
                est.S0 = est[seq(2,3*length(times)-1,3)],
                est.diff = est[seq(3,3*length(times),3)])

  }

  return(out)

}

bootfun <- function(data, indicies,
                    oformula, ofunc, oarg,
                    cformula, cfunc, carg,
                    eformula, earg,
                    method, times, rel.tol,
                    ww){

  weights <- ww
  n <- nrow(data)
  N <- sum(weights)
  #bootstrapping by weights: if(any(indicies != 1:n)) then the current call to
  #bootfun is not for the original data but for a bootstrap replicate. Then,
  #create a bootstrap replicate by assigning new weights according to the
  #distribution of the original weights. Otherwise, just pass on the
  #original weights.
  if(any(indicies != 1:n))
    weights <- rmultinom(n = 1, size = N, prob = weights/N)
  zeros <- (weights == 0)
  if(any(zeros)){
    weights <- weights[-which(zeros)]
    data <- data[-which(zeros), ]
  }
  out <- afSurvival.fit(oformula = oformula, ofunc = ofunc, oarg = oarg,
                        cformula = cformula, cfunc = cfunc, carg = carg,
                        eformula = eformula, earg = earg,
                        method = method, times = times, rel.tol = rel.tol,
                        se.fit = FALSE, data = data, weights = weights)
  return(out$est)

}

afSurvival.fit <- function(oformula, ofunc, oarg,
                           cformula, cfunc, carg,
                           eformula, earg,
                           method, times, rel.tol, jacobian.method,
                           se.fit, data, weights){

  #---SUBROUTINES---

  ML.fun <- function(theta){
    ofit.temp <- ofit
    if(!missing(theta)){
      if(ofunc == "survreg"){
        ofit.temp$coefficients <- theta[1:(nopar-1)]
        ofit.temp$scale <- exp(theta[nopar])
      }
      if(ofunc == "coxph") ofit.temp$coefficients <- theta
    }
    opred <- predict(object = ofit.temp, newdata = data, type = "lp")
    opred0 <- predict(object = ofit.temp, newdata = data0, type = "lp")
    if(ofunc == "survreg"){
      Stres <- 1-survival::psurvreg(q = t, mean = opred, scale = ofit.temp$scale,
                          distribution = ofit.temp$dist)
      St0res <- 1-survival::psurvreg(q = t, mean = opred0, scale = ofit.temp$scale,
                           distribution = ofit.temp$dist)
    }
    if(ofunc == "coxph"){
      #the standardized survival function could be estimated with survexp,
      #but that is MUCH slower
      sst <- survival::survfit(formula = ofit.temp, se.fit = FALSE, censor = FALSE)
      Ststep <- stats::stepfun(sst$time, c(1, sst$surv))
      Stres <- Ststep(t)^exp(opred)
      St0res <- Ststep(t)^exp(opred0)
    }
    out <- matrix(c(Stres, St0res), nrow = n, ncol = 2)
    if(!missing(theta))
      return(apply(X = out, MARGIN = 2, FUN = weighted.mean, w = weights,
                   na.rm=TRUE))
    else return(out)
  }

  IPW.fun <- function(theta){
    cfit.temp <- cfit
    efit.temp <- efit
    if(!missing(theta)){
      if(cfunc == "survreg"){
        cfit.temp$coefficients <- theta[1:(ncpar-1)]
        cfit.temp$scale <- exp(theta[ncpar])
      }
      if(cfunc == "coxph") cfit.temp$coefficients <- theta[1:ncpar]
      efit.temp$coefficients <- theta[(ncpar+1):(ncpar+nepar)]
    }
    cpred <- predict(object = cfit.temp, newdata = data, type = "lp")
    if(cfunc == "survreg")
      Scres <- 1-survival::psurvreg(q = t, mean = cpred, scale = cfit.temp$scale,
                          distribution = cfit.temp$dist)
    if(cfunc == "coxph"){
      #the standardized survival function could be estimated with survexp,
      #but that is MUCH slower
      ssc <- survival::survfit(formula = cfit.temp, se.fit = FALSE, censor = FALSE)
      Scstep <- stats::stepfun(ssc$time, c(1, ssc$surv))
      Scres <- Scstep(t)^exp(cpred)
    }
    p <- predict(object = efit.temp, newdata = data, type = "response")
    Stres <- (data[, U] > t)/Scres
    St0res <- (1-data[, A])*(data[, U] > t)/((1-p)*Scres)
    out <- matrix(c(Stres, St0res), nrow = n, ncol = 2)
    if(!missing(theta))
      return(apply(X = out, MARGIN = 2, FUN = weighted.mean, w = weights,
                   na.rm=TRUE))
    else return(out)
  }

  DR.fun <- function(theta){
    ofit.temp <- ofit
    cfit.temp <- cfit
    efit.temp <- efit
    if(!missing(theta)){
      if(ofunc == "survreg"){
        ofit.temp$coefficients <- theta[1:(nopar-1)]
        ofit.temp$scale <- exp(theta[nopar])
      }
      if(ofunc == "coxph") ofit.temp$coefficients <- theta[1:nopar]
      if(cfunc == "survreg"){
        cfit.temp$coefficients <- theta[(nopar+1):(nopar+ncpar-1)]
        cfit.temp$scale <- exp(theta[nopar+ncpar])
      }
      if(cfunc == "coxph")
        cfit.temp$coefficients <- theta[(nopar+1):(nopar+ncpar)]
      efit.temp$coefficients <- theta[(nopar+ncpar+1):(nopar+ncpar+nepar)]
    }
    opred <- predict(object = ofit.temp, newdata = data, type = "lp")
    opred1 <- predict(object = ofit.temp, newdata = data1, type = "lp")
    opred0 <- predict(object = ofit.temp, newdata = data0, type = "lp")
    cpred <- predict(object = cfit.temp, newdata = data, type = "lp")
    if(ofunc == "survreg")
      St.func <- function(r, opred)
        1-survival::psurvreg(q = r, mean = opred, scale = ofit.temp$scale,
                   distribution = ofit.temp$dist)
    if(ofunc == "coxph"){
      #the standardized survival function could be estimated with survexp,
      #but that is MUCH slower
      sst <- survival::survfit(formula = ofit.temp, se.fit = FALSE, censor = FALSE)
      Ststep <- stats::stepfun(sst$time, c(1, sst$surv))
      St.func <- function(t, pred) Ststep(t)^exp(pred)
    }
    int.val <- vector(length = n)
    if(cfunc  == "survreg"){
      Sc.func <- function(r, cpred)
        1-survival::psurvreg(q = r, mean = cpred, scale = cfit.temp$scale,
                   distribution = cfit.temp$dist)
      int.func <- function(r, opred, cpred){
        fc.r <- survival::dsurvreg(x = r, mean = cpred, scale = cfit.temp$scale,
                         distribution = cfit.temp$dist)
        Sc.r <- Sc.func(r, cpred)
        St.r <- St.func(r, opred)
        return(fc.r/(Sc.r^2*St.r))
      }
      int.val <- vector(length = n)
      for(i in seq(n)){
        int.try <- try(integrate(f = int.func, lower = 0,
                                 upper = min(data[i, U], t), opred = opred[i], cpred = cpred[i],
                                 rel.tol = rel.tol), silent = TRUE)
        if(inherits(int.try, 'try-error')) int.val[i] <- NA
        else int.val[i] <- int.try$value
      }
    }
    if(cfunc == "coxph"){
      #the standardized survival function could be estimated with survexp,
      #but that is MUCH slower
      ssc <- survival::survfit(formula = cfit.temp, se.fit = FALSE, censor = FALSE)
      Scstep <- stats::stepfun(ssc$time, c(1, ssc$surv))
      Sc.func <- function(t, pred) Scstep(t)^exp(pred)
      hc <- diff(-log(c(1, ssc$surv)))
      times.c <- ssc$time
      for(i in seq(n)){
        temp <- times.c < min(data[i, U], t)
        if(!any(temp)) int.val[i] <- 0
        else{
          upper <- max(which(times.c<min(data[i, U], t)))
          hc.i <- hc[1:upper]*exp(cpred[i])
          Sc.i <- Sc.func(times.c[1:upper], cpred[i])
          St.i <-  St.func(times.c[1:upper], opred[i])
          int.val[i] <- sum(hc.i/(Sc.i*St.i))
        }
      }
    }
    p <- predict(object = efit.temp, newdata = data, type = "response")
    H.t <- St.func(t, opred)
    H0.t <- St.func(t, opred0)
    H1.t <- St.func(t, opred1)
    K.t <- Sc.func(t, cpred)
    H.U <- St.func(data[, U], opred)
    K.U <- Sc.func(data[, U], cpred)
    int.t <- H.t*((1-data[, D])*(data[, U]<t)/(K.U*H.U)-int.val)
    int0.t <- (1-data[, A])/(1-p)*int.t
    int1.t <- (data[,A])/p * int.t
    ipw.term.t <- (data[, U] >= t)/K.t
    Stres <- ipw.term.t+int.t
    St0res <- (1-data[, A])/(1-p)*ipw.term.t-
      ((1-data[, A])-(1-p))/(1-p)*H0.t+int0.t
    St1res <- (data[, A])/(p)*ipw.term.t-
      ((data[, A])-(p))/(p)*H1.t+int1.t

    out <- matrix(c(St1res, St0res), nrow = n, ncol = 2)
    if(!missing(theta))
      return(apply(X = out, MARGIN = 2, FUN = weighted.mean, w = weights,
                   na.rm=TRUE))
    else return(out)
  }

  #---PREPARATION---

  n <- nrow(data)
  N <- sum(weights)
  A <- as.character(eformula[[2]])
  data0 <- data1 <- data
  data0[, A] <- 0
  data1[, A] <- 1
  est <- vector(length = length(3*times))
  se <- vector(length = length(3*times))

  if(method == "ML" | method == "DR"){
    D <- as.character(oformula[[2]])[length(as.character(oformula[[2]]))]
    oarg$formula <- oformula
    oarg$data <- data
    oarg$weights <- weights
    if(ofunc == "coxph") oarg$ties <- "breslow" #According to the manual,
    #the default is "efron", but it seems like the coxph function overrides
    #this and always uses "breslow" if the ties arise because of weights.
    #So for consistency, we set ties to "breslow".
    ofit <- do.call(ofunc, oarg)
    if(ofunc == "survreg") nopar <- length(ofit$coefficients)+1
    if(ofunc == "coxph") nopar <- length(ofit$coefficients)
    ofit <- do.call(ofunc, oarg)
    if(se.fit){
      if(ofunc == "coxph") ores <- residuals(object = ofit, type = "score")
      if(ofunc == "survreg"){
        rr <- residuals(object = ofit, type = "matrix")
        dldg <- rr[, "dg"]
        dgdx <- model.matrix(object = ofit)
        ores <- dldg*dgdx
        if(nrow(ofit$var) > length(ofit$coef)) ores <- cbind(ores, rr[, "ds"])
      }
      ovar <- -solve(ofit$var)/N
    }
  }

  if(method == "IPW" | method == "DR"){
    U <- as.character(cformula[[2]])[2]
    carg$formula <- cformula
    carg$data <- data
    carg$weights <- weights
    if(cfunc == "coxph") carg$ties <- "breslow" #According to the manual,
    #the default is "efron", but it seems like the coxph function overrides
    #this and always uses "breslow" if the ties arise because of weights.
    #So for consistency, we set ties to "breslow".
    cfit <- do.call(cfunc, carg)
    if(cfunc == "survreg") ncpar <- length(cfit$coefficients)+1
    if(cfunc == "coxph") ncpar <- length(cfit$coefficients)
    earg$formula <- eformula
    earg$data <- data
    earg$weights <- weights
    earg$family <- "binomial"
    efit <- do.call("glm", earg)
    nepar <- length(efit$coefficients)
    if(se.fit){
      if(cfunc == "coxph") cres <- residuals(object = cfit, type = "score")
      if(cfunc == "survreg"){
        rr <- residuals(object = cfit, type = "matrix")
        dldg <- rr[, "dg"]
        dgdx <- model.matrix(object = cfit)
        cres <- dldg*dgdx
        if(nrow(cfit$var) > length(cfit$coef)) cres <- cbind(cres, rr[, "ds"])
      }
      dldg <- residuals(object = efit, type = "response")
      dgdx <- model.matrix(object = efit)
      eres <- dldg*dgdx

      cvar <- -solve(cfit$var)/N
      evar <- -solve(vcov(object = efit))/N
    }
  }

  #---LOOP OVER TIMES---

  for(k in 1:length(times)){

    t <- times[k]

    if(method == "ML"){
      S <- ML.fun()
      Stres <- S[, 1]
      St0res <- S[, 2]
      if(se.fit){
        res <- cbind(Stres, St0res, ores)
        #this subsetting is needed because cov.wt does not have a na.rm argument
        complete <- complete.cases(res)
        res <- res[complete, ]
        wt <-  weights[complete]
        #weighted covariance is needed because res is created with the residual
        #function, which gives unweighted residuals
        J <- cov.wt(x = res, wt = wt)$cov
        if(ofunc == "survreg")
          theta <- c(ofit$coefficients, log(ofit$scale))
        if(ofunc == "coxph")
          theta <- ofit$coefficients
        SI <- cbind(c(-1, 0), c(0, -1),
                    numDeriv::jacobian(func = ML.fun, x = theta, method = jacobian.method))
        oI <- cbind(rep(0, nopar), rep(0, nopar), ovar)
        I <- rbind(SI, oI)
      }
    }

    if(method == "IPW"){
      S <- IPW.fun()
      Stres <- S[, 1]
      St0res <- S[, 2]
      if(se.fit){
        res <- cbind(Stres, St0res, cres, eres)
        #this subsetting is needed because cov.wt does not have a na.rm argument
        complete <- complete.cases(res)
        res <- res[complete, ]
        wt <-  weights[complete]
        #weighted covariance is needed because res is created with the residual
        #function, which gives unweighted residuals
        J <- cov.wt(x = res, wt = wt)$cov
        if(cfunc == "survreg")
          theta <- c(cfit$coefficients, log(cfit$scale), efit$coefficients)
        if(cfunc == "coxph")
          theta <- c(cfit$coefficients, efit$coefficients)
        SI <- cbind(c(-1, 0), c(0, -1),
                    numDeriv::jacobian(func = IPW.fun, x = theta, method = jacobian.method))
        cI <- cbind(rep(0, ncpar), rep(0, ncpar), cvar,
                    matrix(rep(0, ncpar*nepar), nrow = ncpar, ncol = nepar))
        eI <- cbind(rep(0, nepar), rep(0, nepar),
                    matrix(rep(0, nepar*ncpar), nrow = nepar, ncol = ncpar), evar)
        I <- rbind(SI, cI, eI)
      }
    }

    if(method == "DR"){
      S <- DR.fun()
      St1res <- S[, 1]
      St0res <- S[, 2]
      if(se.fit){
        res <- cbind(St1res, St0res, ores, cres, eres)
        #this subsetting is needed because cov.wt does not have a na.rm argument
        complete <- complete.cases(res)
        res <- res[complete, ]
        wt <-  weights[complete]
        #weighted covariance is needed because res is created with the residual
        #function, which gives unweighted residuals
        J <- cov.wt(x = res, wt = wt)$cov
        theta <- NULL
        if(ofunc == "survreg")
          theta <- c(theta, ofit$coefficients, log(ofit$scale))
        if(ofunc == "coxph")
          theta <- c(theta, ofit$coefficients)
        if(cfunc == "survreg")
          theta <- c(theta, cfit$coefficients, log(cfit$scale))
        if(cfunc == "coxph")
          theta <- c(theta, cfit$coefficients)
        theta <- c(theta, efit$coefficients)
        SI <- cbind(c(-1, 0), c(0, -1),
                    numDeriv::jacobian(func = DR.fun, x = theta, method = jacobian.method))
        oI <- cbind(rep(0, nopar), rep(0, nopar), ovar,
                    matrix(rep(0, nopar*ncpar), nrow = nopar, ncol = ncpar),
                    matrix(rep(0, nopar*nepar), nrow = nopar, ncol = nepar))
        cI <- cbind(rep(0, ncpar), rep(0, ncpar),
                    matrix(rep(0, ncpar*nopar), nrow = ncpar, ncol = nopar), cvar,
                    matrix(rep(0, ncpar*nepar), nrow = ncpar, ncol = nepar))
        eI <- cbind(rep(0, nepar), rep(0, nepar),
                    matrix(rep(0, nepar*nopar), nrow = nepar, ncol = nopar),
                    matrix(rep(0, nepar*ncpar), nrow = nepar, ncol = ncpar), evar)
        I <- rbind(SI, oI, cI, eI)
      }
    }

    #---ESTIMATES---

    St1 <- weighted.mean(x = St1res, w = weights, na.rm = TRUE)
    St0 <- weighted.mean(x = St0res, w = weights, na.rm = TRUE)
    est[((k-1)*3+1):((k-1)*3+3)] <- c(St1, St0, St1 - St0)

    #---STANDARD ERRORS---
    if(se.fit){
      V <- (solve(I)%*%J%*%t(solve(I))/N)[1:2, 1:2]
      #dAF.dSt.dSt0 <- matrix(c(-(1-St0)/(1-St)^2, 1/(1-St)), nrow = 2, ncol = 1)
      se[((k-1)*3+1):((k-1)*3+3)] <- c(sqrt(V[1,1]), sqrt(V[2,2]),
                                       sqrt(V[1,1] + V[2,2] - 2 * V[1,2]))
    }
  }

  if(se.fit) return(list(est = est, se = se))
  else return(list(est = est))

}


