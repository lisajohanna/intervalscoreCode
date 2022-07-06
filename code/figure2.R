# =============================================================================
# Parameter setup
# =============================================================================
conf.level <- 0.95
n <- 50
x <- 0:n
nbpoints <- 999
p <- seq(0, 1, length = nbpoints + 2)[-c(1, nbpoints + 2)]
p_trafo <- asin(sqrt(p)) - (asin(sqrt(0.5))-0.5)

# =============================================================================
# Plot
# =============================================================================
ci <- t(sapply(x, waldCI, n = n, conf.level = conf.level))
y <- sapply(p, meanwidth, n = n, ci = ci, conf.level = conf.level)

df <- data.frame("x" = p, "y" = y)
df$group <- "Uniform scale"

df_trafo <- data.frame("x" = p_trafo, "y" = y)
df_trafo$group <- "Variance-stabilized scale"

df <- rbind(df, df_trafo)

print(
ggplot(df, aes(x = x, y = y, group = group, fill = group)) +
  geom_line(size = 0.5) + 
  geom_ribbon(data = df, aes(x = x, ymax = y), ymin = 0, alpha = 0.3) +
  scale_fill_manual(name = '', 
                    values = c("Uniform scale" = "blue", 
                               "Variance-stabilized scale" = "red")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "#0000FF70"),
        axis.text.x.top = element_text(colour = "#FF000070")) +
  theme(legend.position = "top") +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ . + (asin(sqrt(0.5))-0.5), 
                                         breaks = c(0, pi/4, pi/2),
                                         labels = c("0.0", 
                                                    expression(pi/4), 
                                                    expression(pi/2)),
                                         name = expression(
                                           paste("arcsin(", sqrt(pi), ")", 
                                                 " - arcsin(", sqrt(0.5), 
                                                 ") + 0.5")))) +
  labs(x = expression(pi), y = "EW")
)
