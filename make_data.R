source("stochastic.R")

# Total population
total_pop <- function(sol, dt = 1, strains = c("N_S", "N_A", "N_B", "N_AB")) {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep) %>%
    summarise(total = sum(value) * dt) %>%
    pull(total)
}

# Final population
final_pop <- function(sol, strains = c("N_S", "N_A", "N_B", "N_AB")) {
  sol %>%
    filter(time == max(time) & variable %in% strains) %>%
    group_by(rep) %>%
    summarise(final = sum(value)) %>%
    pull(final)
}

# Time to reach target
target_time <- function(sol, target = 1, strains = c("N_S", "N_A", "N_B", "N_AB")) {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep, time) %>%
    summarise(total = sum(value), .groups = "drop") %>%
    group_by(rep) %>%
    summarise(t = time[min(which(diff(sign(total - target)) != 0), Inf)]) %>%
    pull(t)
}

# Whether target was hit
target_hit <- function(sol, target = 1, strains = c("N_S", "N_A", "N_B", "N_AB")) {
  !is.na(target_time(sol, target, strains))
}

log_plot(simulate(first = "A", rep = 10, deterministic = F)[[1]])
# Create a grid of parameter combinations
# summary <- expand.grid(dose_rep = seq(1, 12, 1), kappa = seq(0.5, 3, 0.5))
summary <- expand.grid(cA = seq(0.5, 1, 0.05), m_A = c(2 ^ - seq(1, 6, 1), 0), first = c("A"))
summary$m_A_discrete <- as.factor(round(summary$m_A, 3))
data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- simulate(
        init = c(N_S = 1e9, N_A = 0, N_B = 0, N_AB = 0),
        R0 = 1e9,
        k = 0,
        alpha = 0,
        supply = 1e8,
        mu = 1,
        first = summary$first[i],
        bcidal_A = 2,
        bcidal_B = 2,
        delta = 0,
        time = 120,
        tau = 1e4,
        max_step = 1e-1,
        rep = 1e3,
        HGT = 0,
        dose_rep = 1,
        dose_gap = 8,
        influx = 5 * c(C_A = summary$cA[i], C_B = 1 - summary$cA[i]),
        cycl = TRUE,
        m_A = 1e-9 * summary$m_A[i], m_B = 1e-9 * (1 - summary$m_A[i]),
        d_A = 0.2, d_B = 0.2,
        deterministic = FALSE
    )
    wins <- !target_hit(data[[i]][[1]], target = 1e3, strains = c("N_A", "N_B"))
    summary$wins[i] <- mean(wins)
    ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
    print(i / nrow(summary))
}
log_plot(data[[44]][[1]][data[[1]][[1]]$rep <= 50, ], use = c("N_S", "N_A", "N_B", "N_AB", "R"))
summary[i, ]
# plot with one independent variable and one dependent variable
ggplot(summary[summary$m_A <= 0.5, ], aes(x = cA, y = wins, color = m_A_discrete)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    theme_light() +
    # scale_y_continuous(limits = c(0, 1)) +
    labs(
        x = "proportion of dose that is A",
        y = "P(extinct)",
        # shape = "m_A",
        color = "m_A"
    ) +
    theme(
        plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.position = "bottom"
    )

# plot with two independent variables and one dependent variable
ggplot(summary, aes(x = log2(m_A), y = cA)) +
    geom_tile(aes(fill = wins)) +
    labs(fill = "P(extinct)") +
    theme_minimal() +
    labs(
        x = "Log2 proportion of de novo mutations that are A",
        y = "Proportion of dose that is drug A",
        size = "Probability that MDR emerges"
    ) +
    scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
    theme(
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.spacing.x = unit(1, "cm"),
        strip.text = element_text(size = 25, face = "bold")
    )

summary
sol <- data[[1]][[1]][data[[1]][[1]]$rep == 43, ]
sol2 <- sol[sol$variable =="N_S" & sol$time <= 15,]
plot(sol2$time, sol2$value, type = "l")
single_run(config, 1)[, c("time", "N_S", "N", "A1", "A2")]
sol
# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    summary$wins[i] <- mean(final_pop(data[[i]][[1]], c("N_A", "N_B")) > 1e2)
    rep <- data[[i]][[2]]$rep
    ci <- binom.test(summary$wins[i] * rep, rep, conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}

for (i in seq_len(nrow(summary))) {
    pops <- total_pop(data[[i]][[1]], strains = "N_S")
    summary$pop[i] <- mean(pops)
    ci <- t.test(pops, conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}

for (i in seq_len(nrow(summary))) {
    wins <- target_hit(data[[i]][[1]], target = 1e2, strains = "N_A")
    summary$wins[i] <- mean(wins)
    ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}




# Next I will make a plot similar to that above, but with a colour scale
ggplot(summary, aes(x = D, y = HGT, color = wins)) +
    geom_point(size = 10) +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal() +
    labs(
        title = "Impact of D and HGT on MDR emergence",
        x = "Dilution fraction D",
        y = "Recombination rate",
        color = "Probability that MDR emerges"
    ) +
    theme(
        legend.position = "bottom",
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.size = unit(2, "cm")
    )

# # Save the simulation results to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/31_08_mut_vary_cycl.RData")
# read the data file we just wrote back as data
load("C:/Users/s4528540/Downloads/results/17_07_synergy_reciprocal.RData")