# =============================================================================
# Parameter setup
# =============================================================================
nbpoints <- 999
p <- seq(0, 1, length = nbpoints + 2)[-c(1, nbpoints + 2)]
conf.level <- 0.95
n_values <- c(10, 25, 50, 100)
nbparas <- length(n_values)

# =============================================================================
# Plot
# =============================================================================
asympt_refs <- function(n, conf.level, nbpoints) {
  p <- seq(0, 1, length = nbpoints + 2)[-c(1, nbpoints + 2)]
  y_eis <- sapply(p, refmeanscore, n = n, conf.level = conf.level)
  y_ew <- sapply(p, refmeanwidth, n = n, conf.level = conf.level)
  return(c(y_eis, y_ew))
}

CInames <- factor(c("Expected interval score", "Expected width"), 
                  levels = c("Expected interval score", "Expected width"))
nbmethods <- length(CInames)
typenames <- factor(c("n = 10", "n = 25", "n = 50", "n = 100"), 
                    levels = c("n = 10", "n = 25", "n = 50", "n = 100"))

y <- foreach(n = n_values, .combine = 'c') %dopar% {
  asympt_refs(n = n, conf.level = conf.level, nbpoints = nbpoints)
}

df <- data.frame("x" = rep(p, nbmethods*nbparas),
                 "y" = y,
                 "CI" = rep(rep(CInames, each = nbpoints), times = nbparas),
                 "type" = rep(typenames, each = nbmethods*nbpoints))

print(
ggplot(df, aes(x = x, y = y, group = CI, colour = CI)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "Asymptotical reference value", 
       color = "Evaluation method") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  theme_bw() +
  theme(legend.justification = "top", aspect.ratio = 1) + 
  facet_wrap(~ type)
)
