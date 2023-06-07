# Critique of Wahl papers
source("wahl/wahl_code.R")
library(hypergeo)

# Hypergeometric function
hyper <- function(s, z) {
    return(Re(hypergeo(1, 1 / (1 + s), 1 + 1 / (1 + s), z)))
}

# theoretical rate function
theory <- function(D, r, s) {
    return(-r / log(D) * (hyper(s, (1 - D) / (D * (D ^ s - 1))) -
        D * hyper(s, (1 - D) * (D ^ s) /  (D ^ s - 1))))
}

phi <- function(D, s) {
    return(ifelse(D == 1, 1 / (1 + s), (1 - D) / (D ^ -s - D)))
}

# evaluate a set of simulation results
metric <- function(summary, data) {
    for (i in seq_len(nrow(summary))) {
        # extract the mutation rate and population size
        m1 <- data[[i]][[2]]$m1
        N0 <- data[[i]][[2]]$N0
        # find the time just after the last bottleneck
        endpoint <- data[[i]][[1]] %>%
            filter(variable == "W", rep == 1) %>%
            filter(value - lag(value) < 0) %>%
            tail(1) %>%
            pull(time)
        # find the likelihood of a new mutation at t=0 going extinct 
        current_phi <- phi(data[[i]][[2]]$D, s)
        # count the number of each mutant at the endpoint
        final_counts <- data[[i]][[1]] %>%
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
        summary$rate[i] <- mean(fixation_rate)
        se <- sd(fixation_rate) / sqrt(length(fixation_rate))
        ci <- summary$rate[i] + se * qnorm(c(0.025, 0.975))
        summary$ci_lower[i] <- ci[1]
        summary$ci_upper[i] <- ci[2]
        # Note how far through the analysis we are
        print(i / nrow(summary))
    }
    return(summary)
}

summary <- metric(summary, data2)
summary$theory <- theory(summary$D, r, s)
summary

log_plot(data2[[20]][[1]][data2[[1]][[1]]$rep == 1 & data2[[1]][[1]]$variable == "M100", ])

data[[1]][[1]] %>%
    filter(rep == 1, variable == "W") %>%
    summarise(max(value)) %>%
    pull() / data[[1]][[2]]$N0

s <- 0.1
time <- 50
m1 <- 1e-9
r <- 1.023

# Wahl 1 No resource constraints
summary <- expand.grid(D = 10 ^ - seq(0.1, 4.0, by = 0.1))
summary$N0 <- summary$D ^ - 0.5 * 10 ^ 9
data <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    N0 <- summary$N0[i]
    # m1 <- summary$m1[i]
    data[[i]] <- simulate(
        seed = NULL,
        rep = 2e2,
        time = time,
        dt = 1e-1,
        tau = - log(D),
        D = D,
        N0 = N0,
        k = 0,
        alpha = 0,
        mu = c(W = r, M = r * (1 + s)),
        init_W = round(N0 * D),
        num_mutants = 1e2
    )
    print(i / nrow(summary))
}

# Wahl 2 Constant resource concentration in dilution media
r <- 1.1
N0 <- 1e9
k_ratio <- 1e-2
summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.6), tau = 24 * 2 ^ - seq(0, 3, by = 1))
summary$m1 <- summary$D ^ - 0.5 * 1e-9
data2 <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    tau <- summary$tau[i]
    m1 <- summary$m1[i]
    k <- k_ratio * N0
    data2[[i]] <- simulate(
        seed = NULL,
        rep = 1e1,
        time = time,
        dt = 1e-3,
        tau = tau,
        D = D,
        N0 = N0,
        m1 = m1,
        k = k,
        alpha = 1,
        mu = c(W = r, M = r * (1 + s)),
        init_W = round(N0 * D),
        num_mutants = 1e2
    )
    print(i/nrow(summary))
}

# Wahl 3 Constant total resource supply per time
summary <- data.frame(D = 10 ^ - seq(0.05, 4, by = 0.2))
# summary$D <- exp(-summary$tau)
data3 <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    data3[[i]] <- simulate(
        seed = NULL,
        rep = 1e3,
        time = time,
        tau = - log(D),
        D = D,
        N0 = N0 * log(D) / (D - 1),
        k = 1e-2 * N0,
        alpha = 1,
        mu = c(W = r, M = r * (1 + s)),
        init = c(W = round(N0 * D), M = 0)
    )
    print(i/nrow(summary))
}

summary <- metric(summary, data2)
summary$k_ratio <- factor(round(summary$k_ratio, 2))
# png("images/test.png", width = 12, height = 10, units = "in", res = 300)
ggplot(summary, aes(x = D, y = rate, group = k_ratio, color = k_ratio)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)) +
    # geom_line(aes(y = theory, color = "new")) +
    # scale_color_manual(values = c("new" = "red")) +
    theme_light() +
    scale_x_log10() +
    scale_y_log10() +
    # scale_y_continuous(limits = c(0,1)) +
    labs(
        # title = "Optimal Dilution Ratio (resource unconstrained)",
        x = "mutation rate",
        y = "fixation rate (loci per hour per (mutations per generation))",
        color = "k_ratio"
    ) +
    theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )

ggplot(summary, aes(x = D, y = wins)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)) +
    geom_line(aes(y = 1 - phi, color = "theory")) +
    theme_light() +
    scale_x_log10() +
    # scale_y_log10() +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values = c("theory" = "red")) +
    labs(
        title = "Optimal Dilution Ratio (resource unconstrained)",
        x = "D",
        y = "1 - phi: probability of a single mutant at t = 0 eventually fixing",
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
# dev.off()

data5 <- list()
r <- 1.2
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    data5[[i]] <- simulate(
        seed = NULL,
        rep = 1e0,
        dt = 0.01,
        time = time,
        tau = - log(D),
        D = D,
        influx = c(A1 = 0, A2 = 0),
        N0 = N0,
        k = 1e-1 * N0,
        alpha = 1,
        m1 = 1e-9,
        m2 = 0,
        mu = c(W = r, M = r * (1 + s)),
        init = c(W = round(N0 * D), M = 0, R2 = 0, R12 = 0)
    )
    print(i/nrow(summary))
}

summary$N_plus <- 0
summary$N_minus <- 0
for (i in seq_len(nrow(summary))) {
    # Find the time at which the final bottleneck occurred
    sols <- data2[[i]][[1]]
    sols <- sols[sols$variable == "N", ]
    last_bneck <- sols %>%
        filter(rep == 1) %>%
        filter(value > lag(value)) %>%
        slice_tail(n = 1) %>%
        pull(time)
    dt <- 0.1
    summary$N_plus[i] <- mean(sols[near(sols$time, last_bneck), ]$value)
    summary$N_minus[i] <- mean(sols[near(sols$time, last_bneck - dt), ]$value)
}

ggplot(summary[summary$k_ratio == 0.1,]) +
    geom_point(aes(x = D, y = N_minus, color = "Before", shape = k_ratio), size = 3) +
    geom_point(aes(x = D, y = N_plus, color = "After", shape = k_ratio), size = 3) +
    theme_light() +
    scale_x_log10() +
    # scale_y_log10() +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      breaks = 10^seq(0, 20, by = 1),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = c("Before" = "blue", "After" = "red"), 
                       name = "Nutrient pulse: ") +
    # scale_y_continuous(limits = c(0,1)) +
    labs(
        title = "Nutrient use (constant media resource conc)",
        # title = "Nutrient use (constant resource flux)",
        x = "D",
        y = "Nutrient concentration"
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

ggplot(summary, aes(x = D, y = tau)) +
    geom_tile(aes(fill = wins)) +
    scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
    labs(x = "D", y = "tau", fill = "proportion",
            title = "Optimal Dilution Ratio (constant resource flux)",
            subtitle = "proportion of runs where M becomes established") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )

# # Save the simulation results to a file
save(data, file = "C:/Users/s4528540/Downloads/results/data_05_06.RData")
# read the data file we just wrote back as data
load("C:/Users/s4528540/Downloads/results/data4.RData")
