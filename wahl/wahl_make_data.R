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
s <- 0.1
time <- 100
r <- 1.1
N0 <- 1e9
k_ratio <- 1e-2
summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1), tau = 24 * 2 ^ - seq(0, 7, by = 1))
summary$m1 <- summary$D ^ - 0.5 * 1e-9
data2 <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    tau <- summary$tau[i]
    m1 <- summary$m1[i]
    k <- k_ratio * N0
    data2[[i]] <- simulate(
        seed = NULL,
        rep = 1e2,
        time = time,
        dt = 1e-1,
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

# png("images/test.png", width = 12, height = 10, units = "in", res = 300)
ggplot(summary, aes(x = D, y = rate)) +
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

ggplot(summary, aes(x = D, y = tau / 24)) +
    geom_tile(aes(fill = log10(rate))) +
    scale_x_log10() +
    # scale the y axis on a log2 scale
    scale_y_continuous(trans = "log2") +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(x = "D", y = "tau (days)", fill = "fixation rate",
            title = "Optimal Dilution Ratio (constant resource flux)",
            subtitle = "fixation rate (loci per hour per (mutations per generation))") +
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
