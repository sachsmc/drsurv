
<!-- README.md is generated from README.Rmd. Please edit that file -->

# drsurv

<!-- badges: start -->
<!-- badges: end -->

The goal of `drsurv` is to implement several doubly robust estimators
for the survival difference at a given time point and one more complex
doubly robust estimator for the survival curve process. The estimators
are doubly robust in the sense that they are consistent if the censoring
model is correctly specified for censoring and either the outcome model
is correctly specified for confounding or the exposure model is
correctly specified for confounding. See
<https://arxiv.org/abs/2310.16207> for more details and examples.

## Installation

You can install `drsurv` like so:

``` r
remotes::install_github("sachsmc/drsurv")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(drsurv)
#> Loading required package: survival
## basic example code

df <- rotterdam
df$time <- pmin(df$rtime, df$dtime) / 365.25
df$status <- ifelse(df$recur == 1 | df$death == 1, 1, 0)
df$censor <- 1 - df$status
drFit.sjovan <-
 sjovan_survdiff(
   oformula = Surv(time, status) ~ chemo + year + age + meno +
     size + factor(grade) + nodes + pgr + er + hormon,
   ofunc = "survreg",
   cformula = Surv(time, censor) ~ chemo + year + age,
   cfunc = "survreg",
   eformula = chemo ~ year + age + meno + size +
     factor(grade) + nodes + pgr + er + hormon,
   method = "DR",
   times = c(2.5, 5, 7.5),
   se.type = "sandwich",
   data = df
 )
drFit.sjovan
#> $est.S1
#> [1] 0.7661741 0.6032550 0.5224925
#> 
#> $se.S1
#> [1] 0.03440273 0.02906199 0.14520126
#> 
#> $est.S0
#> [1] 0.7190130 0.5515392 0.4670844
#> 
#> $se.S0
#> [1] 0.07008758 0.07427899 0.10500988
#> 
#> $est.diff
#> [1] 0.04716109 0.05171581 0.05540801
#> 
#> $se.diff
#> [1] 0.06863352 0.07611256 0.14153032
```

The other methods are `blanche_survdiff` and `pseudo_survdiff`.
