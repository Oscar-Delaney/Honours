source("wahl/wahl_code.R")
source("wahl/data_analysis.R")
library(hypergeo)
library(scales)
library(scico)
library(gridExtra)


### fig:binomial
D <- 0.9
s <- 0.1
N <- 1e8
rep <- 1e6
set.seed(0)
m <- 1 + rgeom(rep, D ^ (1 + s))
b <- rbinom(rep, m, D)
b2 <- rbinom(rep, N*D, D ^ -(1 + s) / N)

# Estimate the PMF of 'b'
pmf_b <- table(b) / rep
pmf_b2 <- table(b2) / rep

# Create a dataframe for plotting
df <- data.frame(value = as.numeric(names(pmf_b)), pmf_b = as.numeric(pmf_b),
                 pmf_b2 = as.numeric(pmf_b2[names(pmf_b)]))
label <- c("Present Analysis", "Wahl et al., 2002")
# Plot the PMFs
pdf("wahl/figs/binomial.pdf", width = 10, height = 10)
ggplot(df, aes(x = value)) +
  geom_point(aes(y = pmf_b, colour = label[1], shape = label[1]), size = 5) +
  geom_point(aes(y = pmf_b2, colour = label[2], shape = label[2]), size = 5) +
  scale_color_manual(values = setNames(scico(2, palette = "roma"), label)) +
  scale_shape_manual(values = setNames(c(19, 17), label)) +
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
  guides(colour = guide_legend(override.aes = list(shape = c(19, 17)))) +
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
    data[i] <- simulate(
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
    r = 2.1,
    s = 0.1,
    k = 1e9,
    alpha = 1
)[[1]]
}
data2 <- 

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

### fig optimality

s <- 0.1

# resource unconstrained
summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1))
summary$m1 <- summary$D ^ - 0.5 * 1e-9
data <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    m1 <- summary$m1[i]
    N0 <- 1e9
    data <- simulate(
        seed = i,
        rep = 1e3,
        time = 100,
        dt = 1e-1,
        tau = - log(D),
        max_step = Inf,
        D = D,
        N0 = N0,
        k = 0,
        alpha = 0,
        r = 1.023,
        s = s,
        init_W = round(N0 * D),
        num_mutants = 1e2
    )
    summary[i, c("rate", "ci_lower", "ci_upper")] <- metric(data)
    print(i / nrow(summary))
}

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_optimality_unconstrained.rdata")

# resource constrained
summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1))
summary$tau <- - log(summary$D)
summary$m1 <- summary$D ^ - 0.5 * 1e-9
data <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    m1 <- summary$m1[i]
    N0 <- 1e9
    data <- simulate(
        seed = i,
        rep = 1e3,
        time = 100,
        dt = 1e-1,
        tau = - log(D),
        max_step = Inf,
        D = D,
        N0 = N0,
        k = N0,
        alpha = 1,
        r = 1.023 * 4,
        s = s,
        init_W = round(N0 * D),
        num_mutants = 1e2
    )
    summary[i, c("rate", "ci_lower", "ci_upper")] <- metric(data)
    print(i / nrow(summary))
}

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
optimal(summary2, 0.5, 0.1, 1)

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

results_plot(summary)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_tau24.rdata")

### fig:k_variation_optimal

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1), k_ratio = 10 ^ seq(-2, 2, by = 0.5))
summary$tau <- - log(summary$D)
summary$k_ratio <- as.factor(round(as.numeric(summary$k_ratio), 3))
summary <- run_sims(summary, rep = 1e2, r = 2)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_k_variation_optimal.rdata")

results_plot(summary)
results_plot <- function(summary) {
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
    
    if("k_ratio" %in% names(summary)){
        plot <- plot +
            geom_point(aes(x = D, y = rate, color = k_ratio), size = 5) +
            geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper, color = k_ratio), size = 1) +
            scale_color_discrete() +
            theme(legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                legend.position = "bottom"
            )
    } else {
        plot <- plot +
            geom_point(aes(x = D, y = rate), size = 5) +
            geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper), size = 1)
    }

    return(plot)
}
