source("wahl/wahl_code.R")
source("wahl/data_analysis.R")
library(hypergeo)
library(scales)
library(scico)
library(gridExtra)

custom_theme <- theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(1, "cm"),
  )

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
        custom_theme
    return(plot)
}

dynamics <- function(data, part) {
    # Filter data to include only non-zero mutants and remove "N"
    non_zero_data <- data %>%
    group_by(variable) %>%
    filter(sum(value) > 0, variable != "N")

    # Define colors using scico
    unique_vars <- unique(non_zero_data$variable)
    colors <- setNames(scico(length(unique_vars), palette = "roma"), unique_vars)
    colors["W"] <- "black"

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
        custom_theme
    return(p)
}

constrained <- function(summary) {
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
        custom_theme
    return(p)
}

### fig dynamics
data_unconstrained <- simulate(
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
    num_mutants = 9,
    r = 1,
    s = 0.1,
    k = 0,
    alpha = 0
)[[1]]

# save the plot
pdf("wahl/figs/dynamics.pdf", width = 10, height = 10)
dynamics(data_unconstrained, "")
dev.off()

### fig dynamics constrained
data <- list()
for (i in 1:4) {
    k_ratio <- 10 ^ (i - 3)
    data[[i]] <- simulate(
    seed = 2,
    time = 30,
    dt = 1e-3,
    max_step = 1e-3,
    tau = log(10),
    D = 0.1,
    N0 = 1e9,
    m1 = 3e-9,
    init_W = 1e8,
    init_M = 0,
    num_mutants = 9,
    r = 2 * (k_ratio + 1),
    s = 0.1,
    k = k_ratio * 1e9,
    alpha = 1
)[[1]]
}

# save the plot
pdf("wahl/figs/dynamics-constrained.pdf", width = 20, height = 20)
grid.arrange(
    dynamics(data[[1]], "A"),
    dynamics(data[[2]], "B"),
    dynamics(data[[3]], "C"),
    dynamics(data[[4]], "D"),
    ncol = 2
)
dev.off()

### fig:optimal
summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1))
summary$tau <- - log(summary$D)
summary <- run_sims(summary, rep = 1e3, r = 1, s = 0.1, res = FALSE)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_optimality_unconstrained.rdata")

# Define linetypes and palettes for theory lines
labels <- c("theory", "approximation")
linetype_palette <- setNames(c("solid", "dashed"), labels)
color_palette <- setNames(scico(2, palette = "vik"), labels)

# Define theory data
theory_long <- pivot_longer(data.frame(
    D = summary$D,
    theory = theory(summary$D, r, s),
    approximation = approx1_theory(summary$D, r, s)
    ),
cols = -D, names_to = "variable", values_to = "value")

# save the plot
pdf("wahl/figs/optimal.pdf", width = 10, height = 10)
base_plot(summary) +
    geom_point(aes(x = D, y = rate), size = 3) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper), linewidth = 0.8) +
    geom_line(data = theory_long, aes(x = D, y = value, color = variable, linetype = variable), linewidth = 1) +
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = linetype_palette) +
    labs(color = "Model", linetype = "Model")
dev.off()

### fig constrained

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1), tau = 24 * 2 ^ - seq(0, 6.5, by = 0.5))
summary <- run_sims(summary, rep = 1e2, r = 2)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_constrained.rdata")

# save the plot
pdf("wahl/figs/constrained.pdf", width = 10, height = 10)
constrained(summary)
dev.off()

### fig:tau24

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1), tau = 24)
summary <- run_sims(summary, rep = 1e3, r = 2, res = TRUE)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_tau24.rdata")

# save the plot
pdf("wahl/figs/tau24.pdf", width = 10, height = 10)
base_plot(summary) +
    geom_point(aes(x = D, y = rate), size = 3) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper), linewidth = 0.8)
dev.off()

### fig:k_variation_optimal

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1),
    k_ratio = 10 ^ seq(-2, 1, by = 1))
summary$tau <- - log(summary$D)
summary$k_ratio <- as.factor(round(as.numeric(summary$k_ratio), 3))
summary <- run_sims(summary, rep = 1e3, r = 1.2)
theory_data <- data.frame(
    D = unique(summary$D),
    rate = theory(unique(summary$D), r = 1, s = 0.1),
    theory = "unconstrained"
)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_k_variation_optimal.rdata")

# save the plot
pdf("wahl/figs/k_variation_optimal.pdf", width = 10, height = 10)
base_plot(summary) +
    geom_point(aes(x = D, y = rate, color = k_ratio), size = 3) +
    geom_line(aes(x = D, y = rate, color = k_ratio), size = 1) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper, color = k_ratio),
        linewidth = 0.8) +
    geom_line(data = theory_data, aes(x = D, y = rate, linetype = theory), size = 1) +
    scale_color_scico_d(palette = "roma") +
    scale_linetype_manual(values = c("unconstrained" = "dashed"))
dev.off()


### fig:t_distribution
D <- 10^-0.1
r <- 1
s <- 0.1
data <- simulate(seed = 1, time = 50, rep = 1e3, dt = 1e-2, max_step = 1e-2,
    D = D,  s = s, tau = -log(D), N0 = 1e9, init_W = round(1e9 * D),
    m1 = 1e-8, k = 0, alpha = 0, r = r, num_mutants = 1e2)
final_counts <- fixed(data)[[1]]

mutation_times <- data[[1]] %>%
    filter(value == 0) %>%
    group_by(rep, variable) %>%
    summarise(last_zero = max(time)  %% data[[2]]$tau, .groups = "keep")
t_data <- final_counts %>%
    inner_join(mutation_times, by = c("rep", "variable"))

# save the data
save(t_data, file = "C:/Users/s4528540/Downloads/results/fig_t_distribution.rdata")

t_theory <- data.frame(
  t = seq(0, data[[2]]$tau, by = 0.001)
)
t_theory$rate <- rate_at_t(D, r = r, s = s, t = t_theory$t)
t_theory$rate <- t_theory$rate / sum(t_theory$rate * 0.001)

# save the plot
pdf("wahl/figs/t_distribution.pdf", width = 10, height = 10)
ggplot(t_data, aes(x = last_zero)) +
  geom_density(aes(y = ..density.., weight = p_fix), adjust = 1/2, fill = "grey", alpha = 1, bw = 0.01) +
  geom_line(data = t_theory, aes(x = t, y = rate, color = "theory")) +
  scale_color_manual(name = NULL, values = c("theory" = "blue")) +
  labs(
    x = expression(italic(t)),
    y = "probability density"
  ) +
  theme_minimal() +
  custom_theme
dev.off()

### fig:s_distribution
s_data <- final_counts %>%
    inner_join(data[[2]]$s_all, by = c("rep", "variable"))

# save the data
save(s_data, file = "C:/Users/s4528540/Downloads/results/fig_s_distribution.rdata")

s_theory <- data.frame(
  s = seq(0, 1, by = 0.001)
)
s_theory$fix <- s_theory$s * exp(-s_theory$s / s) / s^2
s_theory$arise <- exp(-s_theory$s / s) / s

# save the plot
pdf("wahl/figs/s_distribution.pdf", width = 10, height = 10)
ggplot(s_data, aes(x = value)) +
  # raise p_fix to a high power to ensure mutations actually very likely to fix
  # are the ones that dominate the distribution
  geom_density(aes(y = ..density.., weight = p_fix ^ 20),
    adjust = 1/2, fill = "grey", alpha = 1, bw = 0.01) +
  geom_line(data = s_theory, aes(x = s, y = fix, color = "fixing")) +
  geom_line(data = s_theory, aes(x = s, y = arise, color = "arising")) +
  scale_color_manual(name = "mutations", values = c("fixing" = "blue", "arising" = "red")) +
  labs(
    x = expression(italic(t)),
    y = "probability density"
  ) +
  theme_minimal() +
  custom_theme
dev.off()

# weighted mean of s for mutations that go onto fix
weighted.mean(s_data$value, s_data$p_fix ^ 20)

### fig:binomial
# define parameters
D <- 0.9
N <- 1e8
rep <- 1e6
label <- c("det_poi", "det_bin", "sto_poi", "sto_bin")
shapes <- 15:18
set.seed(0)

# define a tibble with the source and pmf values
pmfs <- tibble(
  source = rep(c( "det_bin", "det_poi", "sto_bin", "sto_poi"), each = rep),
  value = c(
    {e <- D^-(1 + s); f <- floor(e);
        rbinom(rep, f, D) + rbinom(rep, 1, (e - f) * D)},
    rpois(rep, D^-s),
    {m <- 1 + rgeom(rep, D ^ (1 + s)); rbinom(rep, m, D)},
    {m <- 1 + rgeom(rep, D ^ (1 + s)); rpois(rep, m)}
  )
) %>%
  group_by(source, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pmf = n / sum(n))

# save the plot
pdf("wahl/figs/binomial.pdf", width = 10, height = 10)
ggplot(pmfs, aes(x = value, y = pmf, colour = source, shape = source)) +
  geom_jitter(size = 5, width = 0.2, height = 0) +
  scale_color_manual(values = setNames(scico(4, palette = "roma"), label)) +
  scale_shape_manual(values = setNames(shapes, label)) +
  theme_light() +
  scale_x_continuous(breaks = 0:7, limits = c(-0.2, 7.2)) +
  scale_y_continuous(trans = scales::log10_trans(),
        breaks = 10^seq(-7, 0),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(
    x = expression(paste("mutants remaining after bottlenecking, ", italic(M(tau^"+")))),
    y = "probability mass",
    colour = NULL,
    shape = NULL
  ) +
  custom_theme
dev.off()


### fig:methodology

summary <- expand.grid(D = 10 ^ - seq(0.01, 4, by = 0.01), div_tau = c(FALSE, TRUE))
summary$tau <- -log(summary$D) / r
summary <- with(summary, {
    summary$det_bin <- D * (1 - phi(D, s)) * (tau ^ !div_tau)
    summary$det_poi <- 2 * s * D * (log(D)^2) / (tau ^ div_tau)
    summary$sto_bin <- theory(D, r, s) * (tau ^ !div_tau)
    summary$sto_poi <- D * (1 - D ^ s) * (tau ^ !div_tau)
    return(summary)
})

# pivot longer
summary <- summary %>%
    pivot_longer(cols = !c("D", "div_tau", "tau"), names_to = "model", values_to = "rate")

# save the plot
pdf("wahl/figs/methodology.pdf", width = 10, height = 10)
base_plot(summary) +
    geom_line(aes(x = D, y = rate, color = model, linetype = div_tau), linewidth = 1) +
    scale_color_manual(values = setNames(scico(4, palette = "roma"), label)) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    labs(color = "model", linetype = "time-optimized") +
    guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
dev.off()

### fig:ci

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1))
summary$tau <- - log(summary$D)
summary <- run_sims(summary, rep = 2e2, time = 1e3, r = 1.2,
    loci = 4, num_mutants = NULL)

# Save the summary to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/fig_ci.rdata")

# save the plot
pdf("wahl/figs/ci.pdf", width = 10, height = 10)
base_plot(summary) +
    labs(y = "Average mutation frequency",) +
    geom_point(aes(x = D, y = rate), size = 3) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper), linewidth = 0.8)
dev.off()
