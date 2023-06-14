# fig:binomial

D <- 0.9
s <- 0.1
N <- 1e8
rep <- 1e6
m <- 1 + rgeom(rep, D ^ (1 + s))
b <- rbinom(rep, m, D)
b2 <- rbinom(rep, N*D, D ^ -(1 + s) / N)

# Estimate the PMF of 'b'
pmf_b <- table(b) / rep
pmf_b2 <- table(b2) / rep

# Create a dataframe for plotting
df <- data.frame(value = as.numeric(names(pmf_b)), pmf_b = as.numeric(pmf_b),
                 pmf_b2 = as.numeric(pmf_b2[names(pmf_b)]))

# Plot the PMFs
ggplot(df, aes(x = value)) +
  geom_point(aes(y = pmf_b, colour = 'Present Analysis'), size = 3, shape = 1) +
  geom_point(aes(y = pmf_b2, colour = 'Wahl et al., 2002'), size = 3) +
  scale_color_manual(values = c("Present Analysis" = "red", "Wahl et al., 2002" = "blue")) +
  theme_light() +
  scale_y_log10() +
  labs(
    x = "Value",
    y = "Probability mass",
    color = NULL
  ) +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "bottom"
  )
