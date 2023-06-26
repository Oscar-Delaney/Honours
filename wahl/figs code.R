source("wahl/wahl_code.R")
source("wahl/data_analysis.R")
library(hypergeo)
library(scales)
library(scico)
library(gridExtra)


base_plot <- function(summary) {
    plot <- ggplot(summary) +
        scale_x_continuous(trans = log10_trans(),
            breaks = 10^seq(-7, 0),
            labels = trans_format("log10", math_format(10^.x))) +
        scale_y_continuous(trans = log10_trans(),
            breaks = 10^seq(-7, 0),
            labels = trans_format("log10", math_format(10^.x))) +
        labs(
            x = expression(italic("D")),
            y = expression(paste("fixation rate (loci hour"^"-1" ~ italic(mu)^-1 ~ italic(N)^-1*")"))
        ) +
        theme_light() +
        theme(
            axis.title = element_text(size = 25),
            axis.text = element_text(size = 25)
        )
    return(plot)
}

### fig:binomial
# define parameters
D <- 0.9
s <- 0.1
N <- 1e8
rep <- 1e6
set.seed(0)

# define a tibble with the source and pmf values
pmfs <- tibble(
  source = rep(c("det_poi", "det_bin", "sto_poi", "sto_bin"), each = rep),
  value = c(
    rpois(rep, D^-s),
    {e <- D^-(1 + s); f <- floor(e);
        rbinom(rep, f, D) + rbinom(rep, 1, (e - f) * D)},
    {m <- 1 + rgeom(rep, D ^ (1 + s)); rpois(rep, m)},
    {m <- 1 + rgeom(rep, D ^ (1 + s)); rbinom(rep, m, D)}
  )
) %>%
  group_by(source, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pmf = n / sum(n))

label <- c("det_poi", "det_bin", "sto_poi", "sto_bin")
shapes <- 15:18  # example shapes for 4 sources

# Plot the PMFs
pdf("wahl/figs/binomial.pdf", width = 10, height = 10)
ggplot(pmfs, aes(x = value, y = pmf, colour = source, shape = source)) +
  geom_point(size = 5) +
  scale_color_manual(values = setNames(scico(4, palette = "roma"), label)) +
  scale_shape_manual(values = setNames(shapes, label)) +
  theme_light() +
  scale_y_continuous(trans = scales::log10_trans(),
        breaks = 10^seq(-7, 0),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(
    x = expression(paste("mutants remaining after bottlenecking, ", italic(M(tau^"+")))),
    y = "probability mass",
    colour = NULL,
    shape = NULL
  ) +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "bottom"
  )

dev.off()

### fig dynamics
data1 <- simulate(
    seed = 1,
    time = 30,
    dt = 1e-3,
    max_step = 1e-3,
    tau = log(10),
    D = 0.1,
    N0 = 1e9,
    m1 = 3e-9,
    init_W = 1e8,
    init_M = 0,
    # num_mutants = 10,
    loci = 3,
    r = 1,
    s = 0.1,
    k = 0,
    alpha = 0
)[[1]]

data <- list()
for (i in 1:4) {
    k_ratio <- 10 ^ (i - 3)
    data[[i]] <- simulate(
    seed = 1,
    time = 30,
    dt = 1e-3,
    max_step = 1e-3,
    tau = log(10),
    D = 0.1,
    N0 = 1e9,
    m1 = 3e-9,
    init_W = 1e8,
    init_M = 0,
    # num_mutants = 10,
    loci = 3,
    r = 2 * (k_ratio + 1),
    s = 0.1,
    k = k_ratio * 1e9,
    alpha = 1
)[[1]]
}

dynamics <- function(data, part) {
    # Filter data to include only non-zero mutants and remove "N"
    non_zero_data <- data %>%
    group_by(variable) %>%
    filter(sum(value) > 0, variable != "N")

    # Define colors using scico
    unique_vars <- unique(non_zero_data$variable)
    colors <- setNames(scico(length(unique_vars), palette = "roma"), unique_vars)

    # Create plot
    p <- ggplot(non_zero_data, aes(x = time, y = value, color = variable)) +
    geom_line(linewidth = 1.5) +
    scale_color_manual(values = colors) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
        breaks = 10^seq(0, 9),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(
        x = "time (hours)",
        y = "population size",
        color = "Strain",
        fill = "Strain"
    ) +
    annotate("text", x = 0, y = 3e9, size = 15,
           label = part, fontface = "bold", hjust = 0.1, vjust = 0.3) +
    theme_light() +
    theme(
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )
    return(p)
}
# save as pdf
pdf("wahl/figs/dynamics.pdf", width = 20, height = 10)
grid.arrange(dynamics(data1, "A"), dynamics(data2, "B"), ncol = 2)
dev.off()

grid.arrange(
    dynamics(data[[1]], "A"),
    dynamics(data[[2]], "B"),
    dynamics(data[[3]], "C"),
    dynamics(data[[4]], "D"),
    ncol = 2
)


### fig optimality

# resource unconstrained
summary <- expand.grid(D = 10 ^ - seq(0.1, 1, by = 0.1))
summary$tau <- - log(summary$D)
summary <- run_sims(summary, rep = 1e3, res = FALSE)
optimal(summary, 1, 0.1, "")

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_optimality_unconstrained.rdata")

# resource constrained
summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1))
summary$tau <- - log(summary$D)
summary <- run_sims(summary, rep = 1e3, res = TRUE, r = 2)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_optimality_constrained.rdata")

optimal <- function(summary, r, s, part) {
    # Define linetypes for the theory lines
    linetype_palette <- c("theory" = "solid", "approx1" = "dotted", "wahl1" = "solid", "wahl2" = "dotted")

    # Define a palette for the theory lines
    color_palette <- c("theory" = scico(2, palette = "vik")[1], 
                    "approx1" = scico(2, palette = "vik")[1], 
                    "wahl1" = scico(2, palette = "vik")[2], 
                    "wahl2" = scico(2, palette = "vik")[2])

    # Define theory data
    theory_long <- pivot_longer(data.frame(
    D = summary$D,
    theory = theory(summary$D, r, s),
    approx1 = approx1_theory(summary$D, r, s),
    wahl1 = 2 * s * summary$D * (log(summary$D) ^ 2),
    wahl2 = 2 * s * summary$D * -log(summary$D)
    ), cols = -D, names_to = "variable", values_to = "value")

    # Adjust the levels of variable to match the order in color_palette and linetype_palette
    theory_long$variable <- factor(theory_long$variable, levels = names(color_palette))


    # Plotting
    p <- ggplot() +
    geom_point(data = summary, aes(x = D, y = rate), color = "black", size = 5) +
    geom_errorbar(data = summary, aes(x = D, ymin = ci_lower, ymax = ci_upper), color = "black", size = 1) +
    geom_line(data = theory_long, aes(x = D, y = value, color = variable, linetype = variable), linewidth = 2) +
    scale_color_manual(values = color_palette, labels = c("Theory", "Approximation", "Wahl original", "Wahl updated")) +
    scale_linetype_manual(values = linetype_palette, labels = c("Theory", "Approximation", "Wahl original", "Wahl updated")) +
    scale_x_continuous(trans = scales::log10_trans(),
        breaks = 10^seq(-7, 0),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(trans = scales::log10_trans(),
        breaks = 10^seq(-7, 0),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(
        x = expression(italic("D")),
        y = expression(paste("fixation rate (loci hour"^"-1" ~ italic(mu)^-1 ~ italic(N)^-1*")")),
        color = "Model", # This will be the title of the unified legend
        linetype = "Model"
    ) +
    guides(
        color = guide_legend(override.aes = list(linetype = rep(c("solid", "dotted"), 2)), ncol = 2),
        linetype = "none"
    ) +
    annotate("text", x = 1e-4, y = 1e-1, size = 15,
           label = part, fontface = "bold", hjust = 0, vjust = 1) +
    theme_light() +
    theme(
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.key.size = unit(2, "cm"),
        legend.spacing.x = unit(0.5, "cm"),
        legend.text.align = 0
    )
    return(p)
}


# load the saved file
load("C:/Users/s4528540/Downloads/results/fig_optimality_unconstrained.rdata")
summary1 <- summary
load("C:/Users/s4528540/Downloads/results/fig_optimality_constrained.rdata")
summary2 <- summary
# save as pdf
pdf("wahl/figs/optimal.pdf", width = 20, height = 10)
grid.arrange(optimal(summary1, 1, 0.1, "A"),
    optimal(summary2, 1, 0.1, "B"), ncol = 2)
dev.off()


### fig constrained

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.4), tau = 24 * 2 ^ - seq(0, 6, by = 1))
summary <- run_sims(summary, rep = 1e2, r = 2)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_constrained.rdata")
load("C:/Users/s4528540/Downloads/results/fig_constrained.rdata")

p <- ggplot(summary, aes(x = D, y = tau / 24)) +
    geom_tile(aes(fill = log10(rate))) +
    scale_x_continuous(trans = scales::log10_trans(),
        breaks = 10^seq(-4, 0),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(trans = scales::log2_trans(),
        breaks = 2^seq(-6, 0),
        labels = scales::trans_format("log2", scales::math_format(2^.x))) +
    scale_fill_gradient(low = "white", high = "blue",
        breaks = pretty_breaks(n = 2), labels = scales::math_format(10^.x)) +
    labs(x = expression(italic("D")), 
        y = expression(paste(italic(tau), " (days)")), 
        fill = "fixation rate") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 30, hjust = 0.5),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.position = "bottom",
        legend.key.size = unit(2, "cm"),
        legend.spacing.x = unit(0.5, "cm")
    )

print(p)

# QQQ make panel B later with more zoomed in view

### fig:tau24

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1), tau = 24)
summary <- run_sims(summary, rep = 1e2, r = 2)

plot <- base_plot(summary) +
    geom_point(aes(x = D, y = rate), size = 5) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper), size = 1)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_tau24.rdata")
load("C:/Users/s4528540/Downloads/results/fig_tau24.rdata")

### fig:k_variation_optimal

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1), k_ratio = 10 ^ seq(-2, 1, by = 1))
summary$tau <- - log(summary$D)
summary$k_ratio <- as.factor(round(as.numeric(summary$k_ratio), 3))
summary <- run_sims(summary, rep = 1e2, r = 2)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_k_variation_optimal.rdata")
load("C:/Users/s4528540/Downloads/results/fig_k_variation_optimal.rdata")

plot <- base_plot(summary) +
    geom_point(aes(x = D, y = rate, color = k_ratio), size = 5) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper, color = k_ratio), size = 1) +
    scale_color_discrete() +
    theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )


### fig:s_distribution

summary1 <- expand.grid(D = 5 * 10 ^ - (1:4))
summary1$tau <- - log(summary1$D)
summary2 <- summary1
summary1 <- run_sims(summary1, rep = 1e1, res = FALSE, func = s_dist)
summary1$res <- FALSE
summary2 <- run_sims(summary2, rep = 1e1, res = TRUE, r = 2, func = s_dist)
summary2$res <- TRUE
summary <- rbind(summary1, summary2)

plot <- base_plot(summary)


### fig:t_distribution

### fig:fixes

r <- 1
s <- 0.1
summary <- expand.grid(D = 10 ^ - seq(0.01, 4, by = 0.01), div_tau = c(FALSE, TRUE))
summary$tau <- -log(summary$D) / r
summary <- with(summary, {
    summary$det_poi <- 2 * s * D * (log(D)^2) / (tau ^ div_tau)
    summary$det_bin <- D * (1 - phi(D, s)) * (tau ^ !div_tau)
    summary$sto_poi <- D * (1 - D ^ s) * (tau ^ !div_tau)
    summary$sto_bin <- theory(D, r, s) * (tau ^ !div_tau)
    return(summary)
})

# pivot longer
summary <- summary %>%
    pivot_longer(cols = !c("D", "div_tau", "tau"), names_to = "model", values_to = "rate")

ggplot(summary) +
    geom_line(aes(x = D, y = rate, color = model, linetype = div_tau), linewidth = 2) +
    scale_x_continuous(trans = log10_trans(),
        breaks = 10^seq(-7, 0),
        labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(trans = log10_trans(),
        breaks = 10^seq(-7, 0),
        labels = trans_format("log10", math_format(10^.x))) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    labs(
        x = expression(italic("D")),
        y = expression(paste("fixation rate (loci hour"^"-1" ~ italic(mu)^-1 ~ italic(N)^-1*")")),
        color = "model",
        linetype = expression(paste("optimizing for ", italic(tau)))
    ) +
    theme_light() +
    theme(
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )
