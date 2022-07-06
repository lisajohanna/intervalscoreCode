# =============================================================================
# Parameter setup
# =============================================================================
nbpoints <- 100
p <- seq(0, 0.5, length = nbpoints + 1)[-1]
conf.level <- 0.95
n_values <- c(10, 25, 50, 100)
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
typenames <- factor(c("n = 10", "n = 25", "n = 50", "n = 100"), 
                    levels = c("n = 10", "n = 25", "n = 50", "n = 100"))

y <- foreach(n = n_values, .combine = 'c') %dopar% {
  CIpropmeasures(n = n, conf.level = conf.level, nbpoints = nbpoints)[, "score"]
}

# difference to the asymptotical reference
ref <- foreach(n = n_values, .combine = 'c') %dopar% {
  r <- sapply(p, refmeanscore, n = n, conf.level = conf.level)
  return(rep(r, times = nbmethods))
}

y <- y - ref

# change y-scales manually
l <- nbmethods*nbpoints
y_modif <- y[1:l]
too_large <- which(y_modif > 0.5)
y_modif[too_large] <- rep(NA, times = length(too_large))
y[1:l] <- y_modif

y_modif <- y[(l+1):(2*l)]
too_large <- which(y_modif > 0.25)
y_modif[too_large] <- rep(NA, times = length(too_large))
y[(l+1):(2*l)] <- y_modif

y_modif <- y[(2*l+1):(3*l)]
too_large <- which(y_modif > 0.125)
y_modif[too_large] <- rep(NA, times = length(too_large))
y[(2*l+1):(3*l)] <- y_modif

y_modif <- y[(3*l+1):(4*l)]
too_large <- which(y_modif > 0.0625)
y_modif[too_large] <- rep(NA, times = length(too_large))
y[(3*l+1):(4*l)] <- y_modif

df <- data.frame("x" = rep(p, nbmethods*nbparas),
                 "y" = y,
                 "CI" = rep(rep(CInames, each = nbpoints), times = nbparas),
                 "type" = rep(typenames, each = nbmethods*nbpoints))

print(
ggprop(df = df, xlab = expression(paste("True ", pi)), 
       ylab = "EIS - asympt. EIS", 
       layers = TRUE, ylimfree = TRUE) +
  xlim(0, 0.5)
)
