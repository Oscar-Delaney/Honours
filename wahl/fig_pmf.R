source("wahl/wahl_code.R")
library(hypergeo)
library(scales)
library(scico)
library(gridExtra)


### fig:binomial
c("Present Analysis" = scico(2, palette = "vik")[1], 
                                "Wahl et al., 2002" = scico(2, palette = "vik")[2])

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
pdf("wahl/figs/binomial.pdf", width = 10, height = 10)
ggplot(df, aes(x = value)) +
  geom_point(aes(y = pmf_b, colour = 'Present Analysis', shape = 'Present Analysis'), size = 5) +
  geom_point(aes(y = pmf_b2, colour = 'Wahl et al., 2002', shape = 'Wahl et al., 2002'), size = 5) +
  scale_color_manual(values = setNames(scico(2, palette = "roma"),
    c("Present Analysis", "Wahl et al., 2002"))) +
  scale_shape_manual(values = c("Present Analysis" = 19, "Wahl et al., 2002" = 17)) +
  theme_light() +
  scale_y_log10() +
  labs(
    x = expression(italic(M(tau^"+"))),
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
    rep = 1,
    seed = 1,
    deterministic = FALSE,
    time = 30,
    dt = 1e-3,
    max_step = 1e-3,
    tau = log(10),
    D = 0.1,
    N0 = 1e9,
    m1 = 3e-9,
    init_W = 1e8,
    init_M = 0,
    num_mutants = 10,
    r = 1,
    s = 0.1,
    k = 0,
    alpha = 0
)[[1]]

data2 <- simulate(
    rep = 1,
    seed = 1,
    deterministic = FALSE,
    time = 30,
    dt = 1e-3,
    max_step = 1e-3,
    tau = log(10),
    D = 0.1,
    N0 = 1e9,
    m1 = 3e-9,
    init_W = 1e8,
    init_M = 0,
    num_mutants = 10,
    r = 2.1,
    s = 0.1,
    k = 1e9,
    alpha = 1
)[[1]]

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
    geom_line(size = 1.5) +  # Set line thickness here
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

# Probability a new mutant at the beginning of a growth phase will go extinct
phi <- function(D, s) {
    return(ifelse(D == 1, 1 / (1 + s), (1 - D) / (D ^ -s - D)))
}

# Hypergeometric function
hyper <- function(s, z) {
    return(Re(hypergeo(1, 1 / (1 + s), 1 + 1 / (1 + s), z)))
}

# theoretical rate function
theory <- function(D, r, s) {
    return(-r / log(D) * (hyper(s, (1 - D) / (D * (D ^ s - 1))) -
        D * hyper(s, (1 - D) * (D ^ s) /  (D ^ s - 1))))
}

# Approximate rate function
approx1_theory <- function(D, r, s) {
    return(r * s * log(D^-1) / (D^-1 -1))
}

# Even more approximate rate function
approx2_theory <- function(D, r, s) {
    return(r * s * sqrt(D))
}

# evaluate a set of simulation results
metric <- function(data) {
    # extract the mutation rate and population size
    m1 <- data[[2]]$m1
    N0 <- data[[2]]$N0
    # find the time just after the last bottleneck
    endpoint <- data[[1]] %>%
        filter(variable == "W", rep == 1) %>%
        filter(value - lag(value) < 0) %>%
        tail(1) %>%
        pull(time)
    # find the likelihood of a new mutation at t=0 going extinct 
    current_phi <- phi(data[[2]]$D, s)
    # count the number of each mutant at the endpoint
    final_counts <- data[[1]] %>%
        group_by(rep, variable) %>%
        filter(time == endpoint, !(variable %in% c("W", "N"))) %>%
        summarise(final_value = value, .groups = "keep") %>%
        mutate(p_fix = 1 - current_phi ^ final_value)
    # estimate the number of these mutations that will go on to fix
    fixed <- final_counts %>%
        group_by(rep) %>%
        summarise(n = sum(final_value > 1e1), n_hat = sum(p_fix))
    # estimate the fixation rate and store this
    fixation_rate <- fixed$n_hat / endpoint / (N0 * m1)
    se <- sd(fixation_rate) / sqrt(length(fixation_rate))
    ci <- mean(fixation_rate) + se * qnorm(c(0.5, 0.025, 0.975))
    return(ci)
}

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
        seed = NULL,
        rep = 2e2,
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
# load the saved file
load("C:/Users/s4528540/Downloads/results/fig_optimality_unconstrained.rdata")


# Define theory data
theory_data <- data.frame(
  D = summary$D,
  theory = theory(summary$D, 1, s),
  approx1 = approx1_theory(summary$D, 1, s),
  wahl1 = s * summary$D * (log(summary$D) ^ 2),
  wahl2 = s * summary$D * -log(summary$D)
)

# Pivot the theory data to long format
theory_long <- theory_data %>%
    pivot_longer(cols = -D, names_to = "variable", values_to = "value")

# Define color palette
my_colors <- setNames(scico(5, palette = "roma"), 
                      c("full", "approx1", "approx2", "wahl1", "wahl2"))

ggplot() +
  geom_point(data = summary, aes(x = D, y = rate), size = 3) +
  geom_errorbar(data = summary, aes(x = D, ymin = ci_lower, ymax = ci_upper)) +
  geom_line(data = theory_long, aes(x = D, y = value, color = variable), size = 1) +
  scale_color_manual(values = my_colors) +
#   scale_shape_manual(values = my_shapes) +
  theme_light() +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(italic("D")),
    y = expression("fixation rate (loci hour"^-1*"mu"^-1*"N"^-1*")"),
    color = "theory"
  ) +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "bottom"
  )
