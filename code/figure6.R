# =============================================================================
# Parameter setup
# =============================================================================
conf.level <- c(0.9, 0.95, 0.99)
n_values <- c(seq(10, 20, by = 2), seq(25, 50, by = 5), seq(60, 100, by = 10))
nbparas <- length(n_values)

# =============================================================================
# Plot
# =============================================================================
CInames <- factor(c("Wald", "Rindskopf", "Arcsine Wald", "Wilson", 
                    "Agresti-Coull", "Likelihood ratio", "Clopper-Pearson", 
                    "Jeffreys HPD", "Jeffreys equal-tailed", "Uniform HPD", 
                    "Uniform equal-tailed"), 
                  levels = c("Wald", "Rindskopf", "Arcsine Wald", "Wilson", 
                             "Agresti-Coull", "Likelihood ratio", 
                             "Clopper-Pearson", 
                             "Jeffreys HPD", "Jeffreys equal-tailed", 
                             "Uniform HPD", "Uniform equal-tailed"))
nbmethods <- length(CInames)


y <- foreach(n = n_values, .combine = 'c') %dopar% {
  CIpropwscoreintegral(n = n, conf.level = conf.level, nbpoints = 10000, 
                       tol = 0.01)
}

# difference to the asymptotical reference
ref <- foreach(n = n_values, .combine = 'c') %dopar% {
  r <- refmeanwscoreintegral(n = n, conf.level = conf.level, nbpoints = 10000)
  return(rep(r, times = nbmethods))
}

y <- y - ref

df <- data.frame("x" = rep(n_values, each = nbmethods),
                 "y" = y,
                 "CI" = rep(CInames, times = nbparas))

print(
ggprop(df = df, xlab = "n", 
       ylab = "Integral of (EWIS - asympt. EWIS) on uniform scale") +
  scale_x_continuous(breaks = seq(10, 100, 10)) +
  coord_cartesian(ylim = c(-0.25, 0.4))
)


y <- foreach(n = n_values, .combine = 'c') %dopar% {
  CIpropwscoreintegral(n = n, conf.level = conf.level, nbpoints = 10000, 
                       tol = 0.01, varstab.rescale = TRUE)
}

# difference to the asymptotical reference
ref <- foreach(n = n_values, .combine = 'c') %dopar% {
  r <- refmeanwscoreintegral(n = n, conf.level = conf.level, nbpoints = 10000, 
                             varstab.rescale = TRUE)
  return(rep(r, times = nbmethods))
}

y <- y - ref

df <- data.frame("x" = rep(n_values, each = nbmethods),
                 "y" = y,
                 "CI" = rep(CInames, times = nbparas))

print(
ggprop(df = df, xlab = "n", 
       ylab = "Integral of (EWIS - asympt. EWIS) on var.-stab. scale") +
  scale_x_continuous(breaks = seq(10, 100, 10)) +
  coord_cartesian(ylim = c(-0.2, 0.45))
)
