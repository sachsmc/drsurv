#' Estimate the survival probabilities in two treatment groups at several times
#'
#' Doubly robust estimation of the survival probabilities in two exposure groups
#' and their difference using the method of Wang (2018). By default it is
#' assumed that censoring is completely independent, but covariate-dependence
#' can be modeled by either specifying the right side of \code{cformula} and \code{cens.method} which
#' will be passed to \link[eventglm]{cumincglm}.
#'
#' As presented in Wang 2018, let the jackknife pseudo observation for the $i$th
#' subject be \deqn{\hat{S}^i(t) = n \hat{S}(t) - (n - 1)\hat{S}^{-i}(t),} where
#' \eqn{\hat{S}(t)} is an estimator of the survival probability at time \eqn{t}
#' using all \eqn{n} subjects, and \eqn{\hat{S}^{-i}(t)} is an estimator of the
#' survival probability at time \eqn{t} leaving the \eqn{i}th subject out of the
#' estimation procedure. The assumption about censoring is implied by the way
#' these pseudo observations are computed, i.e., these could be inverse
#' probability of censoring weighted estimators where the censoring distribution
#' may depend on covariates, nonparametric estimators stratified on covariates,
#' or simply the Kaplan-Meier estimators if censoring is assumed to be
#' independent of covariates. We then use a working outcome model for the
#' conditional survival probability:
#' \deqn{\log[-\log\{v_i(X_i)\}]=\log[-\log\{p(T > t | X_i, \boldsymbol{Z}_i)\}] =
#' \beta^\star X_i+h(X_i,\boldsymbol{Z}_i;\boldsymbol{\gamma}^\star),} where
#' \eqn{\log(-\log(x))} is the link function that coincides with how a
#' proportional hazards model fits survival at this time point. It is of note
#' that, unlike the Cox model, this estimator does not assume anything about the
#' survival function prior to time point \eqn{t}. The parameters in the above
#' model are estimated by solving the generalized estimating equations with
#' \eqn{\hat{S}^i(t)} used as the
#' outcome variable. Let \eqn{\hat{v}_i(x) = \exp[-\exp\{\hat{\beta}^\star
#' x+h(x,\boldsymbol{Z}_i;\hat{\boldsymbol{\gamma}}^\star)\}]}. A logistic regression model for the propensity of the exposure is assumed and estimated. The doubly robust
#' estimator of the difference in survival probability at time \eqn{t} is
#' \deqn{\frac{1}{n}\sum_{i = 1}^n\left[ \frac{X_i \hat{S}^i(t) - (X_i -
#' g(\boldsymbol{Z}_i,
#' \hat{\boldsymbol{\alpha}}))\hat{v}_i(1)}{g(\boldsymbol{Z}_i,
#' \hat{\boldsymbol{\alpha}})} - \frac{(1 - X_i) \hat{S}^i(t) + (X_i -
#' g(\boldsymbol{Z}_i, \hat{\boldsymbol{\alpha}}))\hat{v}_i(0)}{1 -
#' g(\boldsymbol{Z}_i, \hat{\boldsymbol{\alpha}})}\right].} Although Wang
#' proposed a variance estimator, it was found to be biased upward, and they
#' suggest instead using bootstrap variance estimation. Thus, we implement the
#' bootstrap variance estimator.
#'
#' @param oformula The outcome formula, the left side should be a call to \link[survival]{Surv}
#' @param cformula The censoring model formula, a one sided formula with no outcome
#' @param cens.method Method to use for the censoring model, passed to \link[eventglm]{cumincglm}
#' @param eformula The exposure model formula, which will be fit by logistic regression in a call to \link[stats]{glm}
#' @param times Vector of times at which to estimate the survival probabilities
#' @param data Data frame in which to find the variables in the model formulas
#' @param weights Vector of case weights
#' @param subset Logical vector
#' @param R Number of bootstrap replicates for standard error
#' @param parallel Parallel processing for bootstrap, see \link[boot]{boot}
#' @param cl Optional cluster if parallel processing, see \link[boot]{boot}
#' @param ncpus Number of cpus to use for parallel processing, see \link[boot]{boot}
#'
#' @returns A list with estimated survival probabilities in each group, their difference, and estimated standard errors.
#'
#' @references Wang, J. A simple, doubly robust, efficient estimator for survival functions using pseudo observations. Pharmaceutical Statistics. 2018; 17: 38–48. https://doi.org/10.1002/pst.1834
#'
#' @examples
#' df <- rotterdam
#' df$time <- pmin(df$rtime, df$dtime) / 365.25
#' df$status <- ifelse(df$recur == 1 | df$death == 1, 1, 0)
#' df$censor <- 1 - df$status
#' drFit <-
#'   pseudo_survdiff(
#'     oformula = Surv(time, status) ~ chemo + year + age + meno +
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
pseudo_survdiff <- function(oformula,
                             cformula = ~ 1, cens.method = NULL,
                             eformula,
                             times,
                             data, weights = NULL, subset = NULL,
                            R = 50, parallel = "no", cl = NULL, ncpus = NULL){


  pseudo_estimate <- function(data, indices) {

    data <- data[indices,]
    wmod <-
      glm(
        eformula,
        data = data,
        family = "binomial"
      )

    exposure <- wmod$terms[[2]]
    phatr <- predict(wmod, type = "response")

    df1 <- df0 <- data
    df1[[exposure]] <- 1
    df0[[exposure]] <- 0

    # data$Wprop <- data[[exposure]] / phatr +
    #   (1 - data[[exposure]]) / (1 - phatr)

    out <- matrix(NA, nrow = length(times), ncol = 3)

    for(tt in 1:length(times)) {
      ofitr <-
        cumincglm(
          oformula,
          data = data,
          time = times[tt],
          link = "cloglog"
        )


      dr1rw <-
        ofitr$y * (data[[exposure]] == 1) / phatr -
        predict(ofitr, newdata = df1, type = "response") *
        ((data[[exposure]] == 1) - phatr) / phatr
      dr0rw <-
        ofitr$y * (data[[exposure]] == 0) / (1 - phatr) +
        predict(ofitr, newdata = df0, type = "response") *
        ((data[[exposure]] == 1) - phatr) / (1 - phatr)

      s1 <- 1 - mean(dr1rw)
      s0 <- 1 - mean(dr0rw)
      sdiffs.po <- c(s1 - s0)

      out[tt, ] <- c(s1, s0, sdiffs.po)

    }

    out
  }

  ests <- boot(data, pseudo_estimate,
               R = R, parallel = parallel, ncpus = ncpus, cl = cl)

  basest <- ests$t0
  ses <- matrix(apply(ests$t, 2, sd), nrow = length(times), ncol = 3)

  out <- list(est.S1 = basest[, 1],
              se.S1 = ses[, 1],
              est.S0 = basest[, 2],
              se.S0 = ses[, 2],
              est.diff = basest[, 3],
              se.diff = ses[, 3]
              )
  out

}
