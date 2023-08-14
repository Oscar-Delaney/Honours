source("stochastic.R")

# Total population
total_pop <- function(sol, dt = 1, strains = c("S", "R1", "R2", "R12")) {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep) %>%
    summarise(total = sum(value) * dt) %>%
    pull(total)
}

# Final population
final_pop <- function(sol, strains = c("S", "R1", "R2", "R12")) {
  sol %>%
    filter(time == max(time) & variable %in% strains) %>%
    group_by(rep) %>%
    summarise(final = sum(value)) %>%
    pull(final)
}

# Time to reach target
target_time <- function(sol, target = 1, strains = c("S", "R1", "R2", "R12")) {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep, time) %>%
    summarise(total = sum(value), .groups = "drop") %>%
    group_by(rep) %>%
    summarise(t = time[min(which(diff(sign(total - target)) != 0), Inf)]) %>%
    pull(t)
}

# Whether target was hit
target_hit <- function(sol, target = 1, strains = c("S", "R1", "R2", "R12")) {
  !is.na(target_time(sol, target, strains))
}

# Create a grid of parameter combinations
# summary <- expand.grid(dose_rep = seq(1, 12, 1), kappa = seq(0.5, 3, 0.5))
T <- 1
summary <- expand.grid(cA = seq(0, T, 0.1), m_A = seq(0, 1e-9, 1e-10))
summary$cB <- T - summary$cA
summary$m_A_discrete <- as.factor(summary$m_A)
data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- simulate(
        init = c(S = 1e9, RA = 0, RB = 0, RAB = 0),
        N0 = 5e9,
        k = 1e0,
        alpha = 0,
        supply = 3e8,
        mu = 1,
        bcidal1 = 1,
        bcidal2 = 1,
        time = 40,
        tau = 1e4,
        max_step = 1e-1,
        rep = 1e3,
        HGT = 0,
        dose_rep = 1,
        dose_gap = 1e4,
        influx = 5 * c(A = summary$cA[i], B = summary$cB[i]),
        cycl = FALSE,
        m_A = summary$m_A[i], m_B = 1e-9,
        d1 = 0, d2 = 0,
        deterministic = FALSE
    )
    # log_plot(data[[3]][[1]], use = c("S", "RA", "RB"))
    wins <- target_hit(data[[i]][[1]], target = 1e2, strains = c("RA", "RB"))
    summary$wins[i] <- mean(wins)
    ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
    print(i / nrow(summary))
}
# log_plot(data[[1]][[1]][data[[1]][[1]]$rep <= 10, ])

# plot with one independent variable and one dependent variable
ggplot(summary[summary$m_A == 0e-10,], aes(x = cA / T, y = wins, color = m_A_discrete)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    theme_light() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
        x = "proportion of dose that is A",
        y = "Probability that resistance emerges"
    ) +
    theme(
        plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.position = "bottom"
    )

summary
sol <- data[[1]][[1]][data[[1]][[1]]$rep == 43, ]
sol2 <- sol[sol$variable =="S" & sol$time <= 15,]
plot(sol2$time, sol2$value, type = "l")
single_run(config, 1)[, c("time", "S", "N", "A1", "A2")]
sol
# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    summary$wins[i] <- mean(final_pop(data[[i]][[1]], c("R1", "R2")) > 1e2)
    rep <- data[[i]][[2]]$rep
    ci <- binom.test(summary$wins[i] * rep, rep, conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}

for (i in seq_len(nrow(summary))) {
    pops <- total_pop(data[[i]][[1]], strains = "S")
    summary$pop[i] <- mean(pops)
    ci <- t.test(pops, conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}

for (i in seq_len(nrow(summary))) {
    wins <- target_hit(data[[i]][[1]], target = 1e2, strains = "R1")
    summary$wins[i] <- mean(wins)
    ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}


# plot with two independent variables and one dependent variable
ggplot(summary, aes(x = D, y = HGT, size = wins)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal() +
    labs(
        title = "Impact of D and HGT on MDR emergence",
        x = "Dilution fraction D",
        y = "Recombination rate",
        size = "Probability that MDR emerges"
    ) +
    theme(
        legend.position = "bottom",
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )

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
save(summary, file = "C:/Users/s4528540/Downloads/results/17_07_synergy_2on1.RData")
# read the data file we just wrote back as data
load("C:/Users/s4528540/Downloads/results/17_07_synergy_reciprocal.RData")