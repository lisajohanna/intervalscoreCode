# =============================================================================
# Parameter setup
# =============================================================================
nbpoints <- 1000
p <- seq(0, 0.5, length = nbpoints + 1)[-1]
conf.level <- 0.95
n <- 50

# =============================================================================
# Plot
# =============================================================================
y <- CIpropmeasures(n = n, conf.level = conf.level, nbpoints = nbpoints)
y.local <- CIproplocalcoverage(n = n, conf.level = conf.level, 
                               nbpoints = nbpoints)

CInames <- c("Wald", "Rindskopf", "Arcsine Wald", "Wilson", "Agresti", 
             "Likelihood", "Clopper-Pearson", "JeffreysHPD", 
             "JeffreysET", "uniformHPD", "uniformET")
CInames <- rep(CInames, each = nbpoints)

y_cov <- y[, "coverage"][CInames == "Wald"]
y_width <- y[, "width"][CInames == "Wald"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "Wald"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p1 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Wald", x = "", y = "CP") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p2 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "EW") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "Rindskopf"]
y_width <- y[, "width"][CInames == "Rindskopf"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "Rindskopf"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p3 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Rindskopf", x = "", y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p4 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "Arcsine Wald"]
y_width <- y[, "width"][CInames == "Arcsine Wald"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "Arcsine Wald"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p5 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Arcsine Wald", x = "", y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p6 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "Wilson"]
y_width <- y[, "width"][CInames == "Wilson"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "Wilson"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p7 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Wilson", x = "", y = "CP") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p8 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "EW") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "Likelihood"]
y_width <- y[, "width"][CInames == "Likelihood"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "Likelihood"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p9 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Likelihood ratio", x = "", y = "CP") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p10 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "EW") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "Clopper-Pearson"]
y_width <- y[, "width"][CInames == "Clopper-Pearson"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "Clopper-Pearson"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p11 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Clopper-Pearson", x = "", y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p12 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "Agresti"]
y_width <- y[, "width"][CInames == "Agresti"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "Agresti"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p13 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Agresti-Coull", x = "", y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p14 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "JeffreysET"]
y_width <- y[, "width"][CInames == "JeffreysET"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "JeffreysET"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p15 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Jeffreys equal-tailed", x = "", y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p16 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "JeffreysHPD"]
y_width <- y[, "width"][CInames == "JeffreysHPD"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "JeffreysHPD"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p17 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Jeffreys HPD", x = "", y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p18 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "uniformET"]
y_width <- y[, "width"][CInames == "uniformET"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "uniformET"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p19 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Uniform equal-tailed", x = "", y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p20 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

y_cov <- y[, "coverage"][CInames == "uniformHPD"]
y_width <- y[, "width"][CInames == "uniformHPD"]
df <- data.frame("x" = p, "y" = y_cov)
y.local_cov <- y.local[CInames == "uniformHPD"]
df.local <- data.frame("x" = p, "y" = y.local_cov)
p21 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y), alpha = 0.5) +
  geom_line(data = df.local, aes(x = x, y = y), alpha = 1, color = "black") +
  labs(title = "Uniform HPD", x = "", y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0.8, 1)) +
  geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
  theme_bw()
df$y <- y_width
p22 <- ggplot(df, aes(x = x, y = y)) +
  geom_line(alpha = 0.5) +
  labs(x = expression(paste("True ", pi)), y = "") +
  xlim(0, 0.5) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()

p23 <- ggplot() + 
  theme_void()
p24 <- ggplot() + 
  theme_void()

library(ggpubr)
# should be rescaled to fig.height=12.5, fig.width=9
print(
ggarrange(p1, p3, p5, p13, p2, p4, p6, p14,
          p7, p19, p21, p11, p8, p20, p22, p12, 
          p9, p15, p17, p23, p10, p16, p18, p24, nrow = 6, ncol = 4)
)
