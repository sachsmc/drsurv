#' Estimate the survival probabilities in two treatment groups at several times
#'
#' Doubly robust estimation of the survival probabilities in two exposure groups
#' and their difference using the method of Blanche et al 2023. By default it is
#' assumed that censoring is completely independent, but covariate-dependence
#' can be modeled by either specifying the right side of \code{cformula} which
#' will be passed to a Cox model for censoring, or estimating censoring weights
#' externally and passing them to the \code{cens.weights} argument.
#'
#' As presented in Blanche et al. 2023, we model the causal survival
#' probabilities \eqn{p\{T(x)>t\}}. We assume there is a way to correctly
#' specify for confounding a model for \eqn{p\{T(x)\leq t|Z\}} as in a logistic
#' regression model
#' \deqn{p\{T \leq t | X, \boldsymbol{Z}\} = Q_t(X, \boldsymbol{Z}; \beta^\dagger,
#'  \boldsymbol{\gamma}^\dagger) = \frac{\exp(\beta^\dagger
#'  X+h(X,\boldsymbol{Z};\boldsymbol{\gamma}^\dagger))}{1 + \exp(\beta^\dagger
#'  X+h(X,\boldsymbol{Z};\boldsymbol{\gamma}^\dagger))}.}
#' The working propensity model for exposure is a logistic regression model. In
#' addition to these models, we
#' model the censoring distribution as \eqn{p\{C > u | X, \boldsymbol{Z}\} = G_c(u,
#' X, \boldsymbol{Z}).} If \code{cens.weights} is not NULL, then we estimate \eqn{\hat{G}_c(u, X, \boldsymbol{Z})}
#' by a Cox model. This can be done externally by Aalen's additive hazard model,
#' nonparametric estimators stratified on covariates, or if censoring is assumed
#' independent of covariates, \eqn{G_c} may be estimated with the Kaplan-Meier
#' estimator or any other nonparametric estimator and passed in the
#' \code{cens.weights} argument. Let the estimated inverse
#' probability of censoring weights be \deqn{\hat{U}_i(t) =
#'  \frac{I_{\tilde{T_i} \leq t} \Delta_i}{\hat{G}_c(\tilde{T}_i, X_i,
#'  \boldsymbol{Z}_i)} + \frac{I_{\tilde{T_i}> t}}{\hat{G}_c(t, X_i,
#'  \boldsymbol{Z}_i)}.} Then let \eqn{\hat{\beta}^\dagger,
#' \hat{\boldsymbol{\gamma}}^\dagger} denote the parameter estimates resulting
#' from the inverse probability of censoring weighted, but propensity score
#' unweighted, score equations for the logistic model estimated via maximum
#' likelihood among the subset of individuals who either had the event before
#' time t or remained under observation until at least time t. Then the doubly
#' robust estimator of the causal event probability (one minus survival) under
#' exposure level \eqn{x = 1} is \deqn{\frac{1}{n}\sum_{i = 1}^n
#'  W(1, \boldsymbol{Z}_i; \hat{\boldsymbol{\alpha}})\hat{U}_i(t)\left(I_{T_i \leq
#'  t} - Q_k(1, \boldsymbol{Z}_i, \hat{\beta}^\dagger,
#' \hat{\boldsymbol{\gamma}}^\dagger)\right) + \frac{1}{n}\sum_{i = 1}^n Q_k(1,
#' \boldsymbol{Z}_i, \hat{\beta}^\dagger, \hat{\boldsymbol{\gamma}}^\dagger),
#' } where the term \eqn{I_{T_i \leq t}} is defined to be 0
#' for individuals censored before time \eqn{t}, but it is otherwise fully
#' observed. Similarly, we can substitute in \eqn{x = 0} and take the difference
#' between them to obtain estimates of differences in survival probabilities. A
#' standard error estimator is suggested in Blanche et al. 2023 and implemented
#' in the R package \link[mets]{logitIPCWATE}. This function is basically a
#' wrapper to that function to be consistent with the other estimators.
#'
#'
#' @param oformula The outcome formula, the left side should be a call to \link[timereg]{Event}
#' @param cformula The censoring model formula, a one sided formula with no outcome
#' @param cens.weights A vector of inverse probability of censoring weights. If not NULL, then \code{cformula} will be ignored.
#' @param eformula The exposure model formula, which will be fit by logistic regression in a call to \link[stats]{glm}
#' @param times Vector of times at which to estimate the survival probabilities
#' @param data Data frame in which to find the variables in the model formulas
#' @param weights Vector of case weights
#' @param subset Logical vector
#'
#' @returns A list with estimated survival probabilities in each group, their difference, and estimated standard errors.
#'
#' @references Blanche, Paul Frédéric, Anders Holt, and Thomas Scheike. "On logistic regression with right censored data, with or without competing risks, and its use for estimating treatment effects." Lifetime data analysis 29(2):441-482, 2023.
#'
#' @examples
#' df <- rotterdam
#' df$time <- pmin(df$rtime, df$dtime) / 365.25
#' df$status <- ifelse(df$recur == 1 | df$death == 1, 1, 0)
#' df$censor <- 1 - df$status
#' df$chemo <- factor(df$chemo)
#' drFit <-
#'   blanche_survdiff(
#'     oformula = Event(time, status) ~ chemo + year + age + meno +
#'       size + factor(grade) + nodes + pgr + er + hormon,
#'     cformula = ~ chemo + meno + size,
#'     eformula = chemo ~ year + age + meno + size +
#'       factor(grade) + nodes + pgr + er + hormon,
#'     times = c(2.5, 5, 7.5),
#'     data = df
#'   )
#' drFit
#'
#' @export
blanche_survdiff <- function(oformula,
                            cformula = ~ 1, cens.weights = NULL,
                            eformula,
                            times,
                            data, weights = NULL, subset = NULL){


  fits <- lapply(times, \(t) {

    resin <- logitIPCWATE(
      formula = oformula,
      cens.model = cformula,
      cens.weights = cens.weights,
      treat.model = eformula,
      data = data,
      time = t,
      weights = weights,
      subset = subset
    )

    resin2 <- summary(resin)$ateDR[,1:2]
    resin2[1:2, 1] <- 1-resin2[1:2, 1]
    resin2[3,1] <- -resin2[3,1]
    resin2

  })

  out <- list(est.S1 = rep(NA, length(times)),
              se.S1 = rep(NA, length(times)),
              est.S0 = rep(NA, length(times)),
              se.S0 = rep(NA, length(times)),
              est.diff = rep(NA, length(times)),
              se.diff = rep(NA, length(times)))

  for(i in 1:length(times)) {

    out$est.S1[i] <- fits[[i]][1,1]
    out$est.S0[i] <- fits[[i]][2,1]
    out$est.diff[i] <- fits[[i]][3,1]
    out$se.S1[i] <- fits[[i]][1,2]
    out$se.S0[i] <- fits[[i]][2,2]
    out$se.diff[i] <- fits[[i]][3,2]

  }

  out

}
