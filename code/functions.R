# =============================================================================
# Checks if x is a whole number
# =============================================================================
is.wholenumber <- function(x) {abs(x - round(x)) < .Machine$double.eps^0.5}

# =============================================================================
# Computes a Wald confidence interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the confidence level
# overshoot is truncated to [0,1]
# degenerate zero width interval for x=0 and x=n
# =============================================================================
waldCI <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  q <- qnorm(p = (1 + conf.level)/2)
  p <- x/n
  lower <- max(p - q*sqrt(p*(1 - p)/n), 0)
  upper <- min(p + q*sqrt(p*(1 - p)/n), 1)
  return(cbind(lower, upper))
}

# =============================================================================
# Computes a (modified) logit Wald confidence interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the confidence level
# modification: add 0.5 successes and 0.5 failures
# (otherwise for x=0 and x=n no interval computable)
# =============================================================================
logitwaldCI <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  q <- qnorm(p = (1 + conf.level)/2)
  x <- x + 0.5
  n <- n + 1
  lower <- plogis(log(x/(n-x)) - q*sqrt(1/x + 1/(n - x)))
  upper <- plogis(log(x/(n-x)) + q*sqrt(1/x + 1/(n - x)))
  return(cbind(lower, upper))
}

# =============================================================================
# Computes a variance-stabilized Wald confidence interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the confidence level
# overshoot is truncated to [0,pi/2] on var.-stab. scale
# =============================================================================
varstabwaldCI <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  q <- qnorm(p = (1 + conf.level)/2)
  p <- x/n
  lower <- sin(max(asin(sqrt(p)) - q*sqrt(1/(4*n)), 0))^2
  upper <- sin(min(asin(sqrt(p)) + q*sqrt(1/(4*n)), pi/2))^2
  return(cbind(lower, upper))
}

# =============================================================================
# Computes a Wilson confidence interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the confidence level
# =============================================================================
wilsonCI <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  q <- qnorm(p = (1 + conf.level)/2)
  p <- x/n
  pseudo.est <- (x + q^2/2)/(n + q^2)
  pseudo.se <- sqrt(n)/(n + q^2) * sqrt(p*(1 - p) + q^2/(4*n))
  lower <- pseudo.est - q*pseudo.se
  upper <- pseudo.est + q*pseudo.se
  return(cbind(lower, upper))
}

# =============================================================================
# Computes an Agresti-Coull confidence interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the confidence level
# overshoot is truncated to [0,1]
# =============================================================================
agrestiCI <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  q <- qnorm(p = (1 + conf.level)/2)
  x <- x + 2
  n <- n + 4
  p <- x/n
  lower <- max(p - q*sqrt(p*(1 - p)/n), 0)
  upper <- min(p + q*sqrt(p*(1 - p)/n), 1)
  return(cbind(lower, upper))
}

# =============================================================================
# Computes a Clopper-Pearson confidence interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the confidence level
# =============================================================================
clopperpearsonCI <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  alpha <- 1 - conf.level
  if(x == 0) {
    lower <- 0
    upper <- 1 - (alpha/2)^(1/n)
  }
  else if(x == n) {
    lower <- (alpha/2)^(1/n)
    upper <- 1
  }
  else {
    lower <- qbeta(p = alpha/2, shape1 = x, shape2 = n - x + 1)
    upper <- qbeta(p = 1 - alpha/2, shape1 = x + 1, shape2 = n - x)
  }
  return(cbind(lower, upper))
}

# =============================================================================
# Computes a likelihood confidence interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the confidence level
# numerical solutions for 0<x<n
# =============================================================================
likelihoodCI <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  p <- x/n
  eps <- 1e-12
  loglik <- function(p, x, n) {x*log(p) + (n - x)*log(1 - p)}
  f <- function(theta) {
    loglik(p = theta, x = x, n = n) - loglik(p = p, x = x, n = n) + 
      1/2*qchisq(p = conf.level, df = 1)
  }
  if(x == 0) {
    lower <- 0
    upper <- 1 - exp(-1/2*qchisq(p = conf.level, df = 1)/n)
  }
  else if(x == n) {
    lower <- exp(-1/2*qchisq(p = conf.level, df = 1)/n)
    upper <- 1
  }
  else {
    lower <- uniroot(f, interval = c(eps, p))$root
    upper <- uniroot(f, interval = c(p, 1 - eps))$root
  }
  return(cbind(lower, upper))
}

# =============================================================================
# Computes a Jeffreys equal-tailed credible interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the credibility level
# =============================================================================
jeffreysET <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  alpha <- 1 - conf.level
  lower <- qbeta(p = alpha/2, shape1 = x + 0.5, shape2 = n - x + 0.5)
  upper <- qbeta(p = 1 - alpha/2, shape1 = x + 0.5, shape2 = n - x + 0.5)
  return(cbind(lower, upper))
}

# =============================================================================
# Computes a uniform equal-tailed credible interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the credibility level
# =============================================================================
uniformET <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  alpha <- 1 - conf.level
  lower <- qbeta(p = alpha/2, shape1 = x + 1, shape2 = n - x + 1)
  upper <- qbeta(p = 1 - alpha/2, shape1 = x + 1, shape2 = n - x + 1)
  return(cbind(lower, upper))
}

# =============================================================================
# Helper function for HPD intervals with beta posteriors,
# returns the probability of all points for which the beta density is smaller 
# than h (two tails) as well as the two boundaries of these regions
# - p1 and p2 are the parameters of the beta density
# - h is the function value of the beta density (height)
# =============================================================================
outerdens <- function(h, p1, p2){
  modus <- (p1 - 1)/(p1 + p2 - 2)
  schnitt.l <- uniroot(function(x){dbeta(x, p1, p2) - h}, 
                       interval = c(0, modus))$root
  schnitt.u <- uniroot(function(x){dbeta(x, p1, p2) - h}, 
                       interval = c(modus, 1))$root
  tails <- pbeta(schnitt.l, p1, p2) + pbeta(schnitt.u, p1, p2, 
                                            lower.tail = FALSE)
  return(c(tails, schnitt.l, schnitt.u))
}

# =============================================================================
# Computes a Jeffreys HPD interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the credibility level
# extreme cases x=0, x=n handled differently (no mode in these cases)
# =============================================================================
jeffreysHPD <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  alpha <- 1 - conf.level
  p1 <- x + 0.5
  p2 <- n - x + 0.5
  modus <- (p1 - 1)/(p1 + p2 - 2)
  eps <- 1e-15
  if(x == 0) {
    lower <- 0
    upper <- qbeta(p = conf.level, shape1 = 0 + 0.5, shape2 = n - 0 + 0.5)
  }
  else if(x == n) {
    lower <- qbeta(p = 1 - conf.level, shape1 = n + 0.5, shape2 = n - n + 0.5)
    upper <- 1
  }
  else {
    height <- uniroot(function(h) {
      outerdens(h = h, p1 = p1, p2 = p2)[1] - alpha}, 
      interval = c(eps, dbeta(modus, p1, p2) - 10*eps))[["root"]]
    lower <- outerdens(h = height, p1 = p1, p2 = p2)[2]
    upper <- outerdens(h = height, p1 = p1, p2 = p2)[3]
  }
  return(cbind(lower, upper))
}

# =============================================================================
# Computes a uniform HPD interval
# - x is the number of successes
# - n is the sample size
# - conf.level is the credibility level
# extreme cases x=0, x=n handled differently (no mode in these cases)
# =============================================================================
uniformHPD <- function(x, n, conf.level) {
  stopifnot(is.wholenumber(x), is.wholenumber(n), x <= n, n >= 1, x >= 0,
            conf.level > 0, conf.level < 1)
  alpha <- 1 - conf.level
  p1 <- x + 1
  p2 <- n - x + 1
  modus <- (p1 - 1)/(p1 + p2 - 2)
  eps <- 1e-15
  if(x == 0) {
    lower <- 0
    upper <- qbeta(p = conf.level, shape1 = 0 + 1, shape2 = n - 0 + 1)
  }
  else if(x == n) {
    lower <- qbeta(p = 1 - conf.level, shape1 = n + 1, shape2 = n - n + 1)
    upper <- 1
  }
  else {
    height <- uniroot(function(h) {
      outerdens(h = h, p1 = p1, p2 = p2)[1] - alpha}, 
      interval = c(eps, dbeta(modus, p1, p2) - 10*eps))[["root"]]
    lower <- outerdens(h = height, p1 = p1, p2 = p2)[2]
    upper <- outerdens(h = height, p1 = p1, p2 = p2)[3]
  }
  return(cbind(lower, upper))
}

# =============================================================================
# Returns a list with lower and upper limits for all possible x (nb of success)
# for the CIs: Wald, Logit Wald, Var stab Wald, Wilson, Agresti, Likelihood,
#              Clopper-Pearson, , Jeffreys HPD, Jeffreys equal-tailed,
#              Uniform HPD, Uniform equal-tailed
# - n is the sample size
# - conf.level is the confidence level
# =============================================================================
CIprop <- function(n, conf.level) {
  x <- 0:n
  nwald <- t(sapply(x, waldCI, n = n, conf.level = conf.level))
  lwald <- t(sapply(x, logitwaldCI, n = n, conf.level = conf.level))
  vwald <- t(sapply(x, varstabwaldCI, n = n, conf.level = conf.level))
  wil <- t(sapply(x, wilsonCI, n = n, conf.level = conf.level))
  agr <- t(sapply(x, agrestiCI, n = n, conf.level = conf.level))
  lik <- t(sapply(x, likelihoodCI, n = n, conf.level = conf.level))
  CP <- t(sapply(x, clopperpearsonCI, n = n, conf.level = conf.level))
  jeffET <- t(sapply(x, jeffreysET, n = n, conf.level = conf.level))
  jeffHPD <- t(sapply(x, jeffreysHPD, n = n, conf.level = conf.level))
  unifET <- t(sapply(x, uniformET, n = n, conf.level = conf.level))
  unifHPD <- t(sapply(x, uniformHPD, n = n, conf.level = conf.level))
  return(list(nwald, lwald, vwald, wil, agr, lik, CP, 
              jeffHPD, jeffET, unifHPD, unifET))
}

# =============================================================================
# Computes the width of a confidence interval
# - ci is a confidence interval (vector with lower and upper limit)
# =============================================================================
width <- function(ci) {
  lower <- ci[1]
  upper <- ci[2]
  return(upper - lower)
}

# =============================================================================
# Computes the mean width
# - true.pi is the true proportion
# - n is the sample size
# - ci contains all possible confidence intervals (for all possible x = 0:n), 
#   a matrix where each row contains lower and upper limit for one x
# conf.level is not needed here but in CIpropmeasures() because of the score
# =============================================================================
meanwidth <- function(true.pi, n, ci, conf.level) {
  sum(dbinom(x = 0:n, size = n, prob = true.pi)*apply(ci, 1, width))
}

# =============================================================================
# Computes the coverage of a confidence interval
# - ci is a confidence interval (vector with lower and upper limit)
# - true.pi is the true proportion
# =============================================================================
coverage <- function(ci, true.pi) {
  lower <- ci[1]
  upper <- ci[2]
  return((true.pi >= lower) & (true.pi <= upper))
}

# =============================================================================
# Computes the coverage probability
# - true.pi is the true proportion
# - n is the sample size
# - ci contains all possible confidence intervals (for all possible x = 0:n), 
#   a matrix where each row contains lower and upper limit for one x
# conf.level is not needed here but in CIpropmeasures() because of the score
# =============================================================================
meancoverage <- function(true.pi, n, ci, conf.level) {
  sum(dbinom(x = 0:n, size = n, prob = true.pi)*
        apply(ci, 1, coverage, true.pi = true.pi))
}

# =============================================================================
# Computes the interval score of a confidence interval
# - ci is a confidence interval (vector with lower and upper limit)
# - true.pi is the true proportion
# - conf.level is the corresponding (!) confidence level of ci
# =============================================================================
score <- function(ci, true.pi, conf.level) {
  alpha <- 1 - conf.level
  return(width(ci) + 
           2/alpha*min(abs(true.pi - ci))*(1 - coverage(ci, true.pi)))
}

# =============================================================================
# Computes the mean interval score
# - true.pi is the true proportion
# - n is the sample size
# - ci contains all possible confidence intervals (for all possible x = 0:n), 
#   a matrix where each row contains lower and upper limit for one x
# - conf.level is the corresponding (!) confidence level of ci
# - varstab.rescale decides if the var.-stab. transformation should be used
# =============================================================================
meanscore <- function(true.pi, n, ci, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  sum(dbinom(x = 0:n, size = n, prob = true.pi)*
        apply(ci, 1, score, conf.level = conf.level, true.pi = true.pi))
}
# needs to be vectorized for the numerical integration
meanscore <- Vectorize(meanscore, vectorize.args = c("true.pi"))

# =============================================================================
# Computes the weighted interval score of confidence intervals for different 
# confidence levels
# - ci is a vector with confidence intervals (lower and upper limit), 
#   concatenated for the different levels
# - true.pi is the true proportion
# - conf.level is a vector with the corresponding (!) confidence levels of ci
# =============================================================================
wscore <- function(ci, true.pi, conf.level) {
  ci <- matrix(ci, ncol = 2, byrow = TRUE)
  res <- apply(ci, 1, score, true.pi = true.pi, conf.level = conf.level)
  return(sum(diag(res)))
}

# =============================================================================
# Computes the mean weighted interval score
# - true.pi is the true proportion
# - n is the sample size
# - ci contains all possible confidence intervals (for all possible x = 0:n), 
#   a matrix where each row contains lower and upper limits for one x, 
#   concatenated for all considered levels (per row)
# - conf.level is a vector with the corresponding (!) confidence levels of ci
# - varstab.rescale decides if the var.-stab. transformation should be used
# =============================================================================
meanwscore <- function(true.pi, n, ci, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  sum(dbinom(x = 0:n, size = n, prob = true.pi)*
        apply(ci, 1, wscore, conf.level = conf.level, true.pi = true.pi))
}
# needs to be vectorized for the numerical integration
meanwscore <- Vectorize(meanwscore, vectorize.args = c("true.pi"))

# =============================================================================
# Returns mean width, coverage probability and mean interval score for a finite 
# grid of true proportions and the following confidence intervals: 
# Wald, Logit Wald, Var stab Wald, Wilson, Agresti, Likelihood, Clopper-Pearson, 
# Jeffreys HPD, Jeffreys equal-tailed, Uniform HPD, Uniform equal-tailed 
# (as a concatenated vector)
# - n is the sample size
# - conf.level is the confidence level
# - nbpoints is the number of points in the equidist. grid of true proportions
# measures symmetric around 0.5 --> grid of (0,0.5]
# =============================================================================
CIpropmeasures <- function(n, conf.level, nbpoints) {
  CIs <- CIprop(n, conf.level)
  pvector <- seq(0, 0.5, length = nbpoints + 1)[-1]
  methodwisemeanmeasure <- function(ci, measure.fct) {
    res <- sapply(pvector, measure.fct, n = n, ci = ci, conf.level = conf.level)
    return(res)
  }
  w <- lapply(CIs, methodwisemeanmeasure, measure.fct = meanwidth)
  w <- unlist(w)
  c <- lapply(CIs, methodwisemeanmeasure, measure.fct = meancoverage)
  c <- unlist(c)
  s <- lapply(CIs, methodwisemeanmeasure, measure.fct = meanscore)
  s <- unlist(s)
  return(data.frame("width" = w, "coverage" = c, "score" = s))
}

# =============================================================================
# Returns the integral of the mean interval score over the unit interval of
# underlying true proportions for the following confidence intervals: 
# Wald, Logit Wald, Var stab Wald, Wilson, Agresti, Likelihood, Clopper-Pearson, 
# Jeffreys HPD, Jeffreys equal-tailed, Uniform HPD, Uniform equal-tailed 
# (in a vector with the same order)
# - n is the sample size
# - conf.level is the confidence level
# - nbpoints is the number of subdivisions used for the integration
# - varstab.rescale decides if the var.-stab. transformation should be used
# =============================================================================
CIpropscoreintegral <- function(n, conf.level, nbpoints, 
                                varstab.rescale = FALSE) {
  CIs <- CIprop(n, conf.level)
  if(varstab.rescale == FALSE) {u <- 1}
  if(varstab.rescale == TRUE) {u <- pi/2}
  res <- lapply(CIs, function(ci) integrate(meanscore, n = n, ci = ci, 
                                            conf.level = conf.level, 
                                            varstab.rescale = varstab.rescale,
                                            lower = 0, upper = u, 
                                            subdivisions = nbpoints)$value)
  return(unlist(res))
}

# =============================================================================
# Returns the mean weighted interval scores for a finite grid of true 
# proportions and confidence intervals: Wald, Logit Wald, Var stab Wald, Wilson,
#                                       Agresti, Likelihood, Clopper-Pearson, 
#                                       Jeffreys HPD, Jeffreys equal-tailed,
#                                       Uniform HPD, Uniform equal-tailed
# (as a concatenated vector)
# - n is the sample size
# - conf.level is a vector of the used confidence levels
# - nbpoints is the number of points in the equidist. grid of true proportions
# uses the linearity of the expectation
# =============================================================================
CIpropwscore <- function(n, conf.level, nbpoints) {
  score.matrix <- sapply(conf.level, 
                         function(x) CIpropmeasures(conf.level = x, 
                                                    nbpoints = nbpoints, 
                                                    n = n)[, "score"])
  return(apply(score.matrix, 1, sum))
}

# =============================================================================
# Returns the integral of the mean weighted interval score over the unit 
# interval of underlying true proportions for the confidence intervals: 
# Wald, Logit Wald, Var stab Wald, Wilson, Agresti, Likelihood, Clopper-Pearson, 
# Jeffreys HPD, Jeffreys equal-tailed, Uniform HPD, Uniform equal-tailed
# (in a vector with the same order)
# - n is the sample size
# - conf.level is a vector of the used confidence levels
# - nbpoints is the number of subdivisions used for the integration
# - tol is the relative tolerance (of errors) in the integration
# - varstab.rescale decides if the var.-stab. transformation should be used
# =============================================================================
CIpropwscoreintegral <- function(n, conf.level, nbpoints, tol, 
                                 varstab.rescale = FALSE) {
  res <- lapply(conf.level, CIprop, n = n)
  CIs <- res[[1]]
  for(i in 2:length(res)) {
    CIs <- mapply(cbind, CIs, res[[i]], SIMPLIFY = FALSE)
  }
  if(varstab.rescale == FALSE) {u <- 1}
  if(varstab.rescale == TRUE) {u <- pi/2}
  res <- lapply(CIs, function(ci) integrate(meanwscore, n = n, ci = ci, 
                                            conf.level = conf.level, 
                                            varstab.rescale = varstab.rescale,
                                            lower = 0, upper = u, 
                                            subdivisions = nbpoints,
                                            rel.tol = tol)$value)
  return(unlist(res))
}

# =============================================================================
# Asymptotical reference for expected (weighted) interval score/ expected width
# and integrals thereof
# - true.pi is the true proportion
# - n is the sample size
# - conf.level is the confidence level
# - varstab.rescale decides if the var.-stab. transformation should be used
# - nbpoints is the number of subdivisions used for the integration
# =============================================================================
refmeanscore <- function(true.pi, n, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  q <- qnorm(p = (1 + conf.level)/2)
  s <- sqrt(true.pi*(1-true.pi)/n)
  2*s*dnorm(q)/(1-pnorm(q))
}

refmeanwscore <- function(true.pi, n, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  res <- sapply(conf.level, refmeanscore, true.pi = true.pi, n = n, 
                varstab.rescale = varstab.rescale)
  return(sum(res))
}
# needs to be vectorized for the numerical integration
refmeanwscore <- Vectorize(refmeanwscore, vectorize.args = c("true.pi"))

refmeanwidth <- function(true.pi, n, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  q <- qnorm(p = (1 + conf.level)/2)
  s <- sqrt(true.pi*(1-true.pi)/n)
  2*q*s
}

refmeanscoreintegral <- function(n, conf.level, nbpoints, 
                                 varstab.rescale = FALSE) {
  if(varstab.rescale == FALSE) {u <- 1}
  if(varstab.rescale == TRUE) {u <- pi/2}
  res <- integrate(refmeanscore, n = n, conf.level = conf.level, 
                   varstab.rescale = varstab.rescale,
                   lower = 0, upper = u, 
                   subdivisions = nbpoints)$value
  return(res)
}

refmeanwscoreintegral <- function(n, conf.level, nbpoints, 
                                  varstab.rescale = FALSE) {
  if(varstab.rescale == FALSE) {u <- 1}
  if(varstab.rescale == TRUE) {u <- pi/2}
  res <- integrate(refmeanwscore, n = n, conf.level = conf.level, 
                   varstab.rescale = varstab.rescale,
                   lower = 0, upper = u, 
                   subdivisions = nbpoints)$value
  return(res)
}

# =============================================================================
# Asymptotical reference for variance of interval score and integrals thereof
# - true.pi is the true proportion
# - n is the sample size
# - conf.level is the confidence level
# - varstab.rescale decides if the var.-stab. transformation should be used
# - nbpoints is the number of subdivisions used for the integration
# =============================================================================
refvarscore <- function(true.pi, n, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  q <- qnorm(p = (1 + conf.level)/2)
  s <- sqrt(true.pi*(1-true.pi)/n)
  a <- 1-conf.level
  res <- 4*s^2*(1/a + q^2*(1/a-1) + q*(2-1/a)*dnorm(q)/(1-pnorm(q)) - 
                  (dnorm(q))^2/(1-pnorm(q))^2)
  return(res)
}

normscore <- function(true.pi, n, ci, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  res <- sum(dbinom(x = 0:n, size = n, prob = true.pi)*
               apply(ci, 1, score, conf.level = conf.level, true.pi = true.pi))
  refExp <- refmeanscore(true.pi = true.pi, n = n, conf.level = conf.level, 
                          varstab.rescale = varstab.rescale)
  refVar <- refvarscore(true.pi = true.pi, n = n, conf.level = conf.level, 
                         varstab.rescale = varstab.rescale)
  return((res-refExp)/sqrt(refVar))
}

CIpropnormscores <- function(n, conf.level, nbpoints) {
  CIs <- CIprop(n, conf.level)
  pvector <- seq(0, 0.5, length = nbpoints + 1)[-1]
  methodwisemeanmeasure <- function(ci, measure.fct) {
    res <- sapply(pvector, measure.fct, n = n, ci = ci, conf.level = conf.level)
    return(res)
  }
  res <- lapply(CIs, methodwisemeanmeasure, measure.fct = normscore)
  res <- unlist(res)
  return(res)
}

scaledscore <- function(true.pi, n, ci, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  res <- sum(dbinom(x = 0:n, size = n, prob = true.pi)*
               apply(ci, 1, score, conf.level = conf.level, true.pi = true.pi))
  refVar <- refvarscore(true.pi = true.pi, n = n, conf.level = conf.level, 
                        varstab.rescale = varstab.rescale)
  return(res/sqrt(refVar))
}
# needs to be vectorized for the numerical integration
scaledscore <- Vectorize(scaledscore, vectorize.args = c("true.pi"))

scaledref <- function(true.pi, n, conf.level, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  refExp <- refmeanscore(true.pi = true.pi, n = n, conf.level = conf.level, 
                         varstab.rescale = varstab.rescale)
  refVar <- refvarscore(true.pi = true.pi, n = n, conf.level = conf.level, 
                        varstab.rescale = varstab.rescale)
  return(refExp/sqrt(refVar))
}
# needs to be vectorized for the numerical integration
scaledref <- Vectorize(scaledref, vectorize.args = c("true.pi"))

CIpropscaledscoreintegral <- function(n, conf.level, nbpoints, 
                                      varstab.rescale = FALSE) {
  CIs <- CIprop(n, conf.level)
  if(varstab.rescale == FALSE) {u <- 1}
  if(varstab.rescale == TRUE) {u <- pi/2}
  eps <- 1e-12
  res <- lapply(CIs, function(ci) integrate(scaledscore, n = n, ci = ci, 
                                            conf.level = conf.level, 
                                            varstab.rescale = varstab.rescale,
                                            lower = 0 + eps, upper = u - eps, 
                                            subdivisions = nbpoints)$value)
  return(unlist(res))
}

CIpropscaledrefintegral <- function(n, conf.level, nbpoints, 
                                    varstab.rescale = FALSE) {
  if(varstab.rescale == FALSE) {u <- 1}
  if(varstab.rescale == TRUE) {u <- pi/2}
  eps <- 1e-12
  res <- integrate(scaledref, n = n, conf.level = conf.level, 
                   varstab.rescale = varstab.rescale,
                   lower = 0 + eps, upper = u - eps, 
                   subdivisions = nbpoints)$value
  return(res)
}

# =============================================================================
# Returns local average of coverage probabilities (as in BayarriBerger2004)
# for a finite grid of true proportions and the following confidence intervals: 
# Wald, Logit Wald, Var stab Wald, Wilson, Agresti, Likelihood, Clopper-Pearson, 
# Jeffreys HPD, Jeffreys equal-tailed, Uniform HPD, Uniform equal-tailed 
# (as a concatenated vector)
# - n is the sample size
# - conf.level is the confidence level
# - nbpoints is the number of points in the equidist. grid of true proportions
# measures symmetric around 0.5 --> grid of (0,0.5]
# =============================================================================
CIproplocalcoverage <- function(n, conf.level, nbpoints) {
  CIs <- CIprop(n, conf.level)
  x <- 0:n
  pvector <- seq(0, 0.5, length = nbpoints + 1)[-1]
  a <- function(p) {
    if(p <= 0.025) {
      NA # (p*(1-p)*p^(-2) - 1)*p #1 - 2*0.025
    }
    else if(p >= (1 - 0.025)) {
      NA # (p*(1-p)*(1-p)^(-2) - 1)*p #1/0.025 - 3 + 2*0.025
    }
    else {
      (p*(1-p)*0.025^(-2) - 1)*p
    }
  }
  local.meancoverage <- function(p, ci) {
    ap <- a(p)
    a1mp <- a(1-p)
    alpha <- ap + x
    beta <- a1mp + n - x
    values.gamma <- (lchoose(n, x) 
                     + lgamma(ap + a1mp) - lgamma(ap) - lgamma(a1mp)
                     + lgamma(ap + x) + lgamma(a1mp + n - x) 
                     - lgamma(ap + a1mp + n))
    values.integral <- log(pbeta(ci[, 2], alpha, beta) 
                           - pbeta(ci[, 1], alpha, beta))
    return(sum(exp(values.gamma + values.integral)))
  }
  methodwisemeanmeasure <- function(ci) {
    res <- sapply(pvector, local.meancoverage, ci = ci)
    return(res)
  }
  c <- lapply(CIs, methodwisemeanmeasure)
  return(unlist(c))
}

# =============================================================================
# GIS for HPD intervals (computation of tail probabilities, generalized interval
# score and expected generalized interval score)
# - n is the sample size
# - true.pi is the true proportion
# - prior is the used prior (uniform or Jeffreys)
# - tailprobs: ci is a vector with lower limit, upper limit, nb of success x
# - gscore: ci is a vector with lower limit, upper limit, tail prob. 1 and 2
# - meangscore: ci is a matrix with rows for all x=0:n, vector like in gscore
#   for each row
# - varstab.rescale decides if the var.-stab. transformation should be used
# =============================================================================
tailprobs <- function(ci, n, prior) {
  lower <- ci[1]
  upper <- ci[2]
  x <- ci[3]
  if (prior == "uniform") s <- 1
  if (prior == "jeffreys") s <- 0.5
  if (x == 0) {
    alpha1 <- 0
    alpha2 <- pbeta(upper, x + s, n - x + s, lower.tail = FALSE)
  }
  else if (x == n) {
    alpha1 <- pbeta(lower, x + s, n - x + s)
    alpha2 <- 0
  }
  else {
    alpha1 <- pbeta(lower, x + s, n - x + s)
    alpha2 <- pbeta(upper, x + s, n - x + s, lower.tail = FALSE)
  }
  return(c(alpha1, alpha2))
}

gscore <- function(ci, true.pi) {
  lower <- ci[1]
  upper <- ci[2]
  alpha1 <- ci[3]
  alpha2 <- ci[4]
  if (alpha1 == 0) {
    k <- 1/alpha2
  }
  else if (alpha2 == 0) {
    k <- 1/alpha1
  }
  else {
    k <- ifelse(abs(true.pi - lower) < abs(true.pi - upper), 1/alpha1, 1/alpha2) 
  }
  k <- k*min(abs(c(true.pi - lower, true.pi - upper)))
  return(width(ci) + k*(1 - coverage(ci, true.pi)))
}

meangscore <- function(true.pi, n, ci, varstab.rescale = FALSE) {
  if(varstab.rescale == TRUE) {true.pi <- sin(true.pi)^2}
  sum(dbinom(x = 0:n, size = n, prob = true.pi)*
        apply(ci, 1, gscore, true.pi = true.pi))
}
# needs to be vectorized for the numerical integration
meangscore <- Vectorize(meangscore, vectorize.args = c("true.pi"))

# =============================================================================
# Plot function for binomial proportion
# - df is a dataframe with columns x (x values), y (y values), CI (method of
#   confidence interval), type (sample size or confidence level setting)
# - xlab and ylab are the labels of the x- and y-axis
# - layers decides if we have several plots (types)
# - ylimfree decides if a common y range should be used
# =============================================================================
library(ggplot2)
library(scales)
ggprop <- function(df, xlab, ylab, layers = FALSE, ylimfree = FALSE) {
  colours <- hue_pal()(6)
  colours <- c(colours[c(1,4,3)], "red", colours[c(6,2)], "orange", "grey", 
               "black", colours[5], "blue")
  if(layers == FALSE) {
    relevelnb <- order(df$y[1:11])
    ranking10 <- levels(df$CI)[relevelnb]
    df$CI <- factor(df$CI, levels = ranking10)
    colours <- colours[relevelnb]
  }
  names(colours) <- levels(df$CI)
  p <- ggplot(df, aes(x = x, y = y, group = CI, colour = CI)) +
    geom_line(alpha = 0.5) +
    scale_colour_manual(name = "CI", values = colours) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(legend.justification = "top", aspect.ratio = 1)
  if(layers == TRUE) p <- p + facet_wrap(~ type)
  if(ylimfree == TRUE) p <- p + facet_wrap(~ type, scales = "free")
  return(p)
}
