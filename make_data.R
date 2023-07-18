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
summary <- expand.grid(bcidal1 = seq(0, T, 0.1), bcidal2 = seq(0, T, 0.1))
summary$bstatic1 <- T - summary$bcidal1
summary$bstatic2 <- T - summary$bcidal2
data <- list()
for (i in seq_len(nrow(summary))) {
    data[[1]] <- simulate(
        init = c(S = 3e7, R1 = 0, R2 = 0, R12 = 0),
        N0 = 1e9,
        k = 1e0,
        alpha = 0,
        supply = 3e8,
        mu = 1,
        bcidal1 = summary$bcidal1[i],
        bcidal2 = summary$bcidal2[i],
        bstatic1 = summary$bstatic1[i],
        bstatic2 = summary$bstatic2[i],
        zeta1 = c(S = 1, R1 = 28, R2 = 1, R12 = 28),
        zeta2 = c(S = 1, R1 = 1, R2 = 28, R12 = 28),
        i12 = 2, i21 = 2,
        delta = 0.1,
        time = 60,
        tau = 1e4,
        max_step = 1e-1,
        rep = 1e3,
        HGT = 0,
        dose_rep = 1,
        dose_gap = 10,
        influx = 3 * c(A1 = 1, A2 = 1),
        cycl = TRUE,
        m1 = 1e-9, m2 = 1e-9,
        d1 = 0.1, d2 = 0.1,
        deterministic = FALSE
    )
    wins <- target_hit(data[[1]][[1]], target = 1e2, strains = c("R1", "R2"))
    summary$wins[i] <- mean(wins)
    ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
    print(i / nrow(summary))
}
# log_plot(data[[1]][[1]][data[[1]][[1]]$rep <= 10, ], use = c("S", "R1", "R2", "N"))

ggplot(summary, aes(x = bcidal1 / T, y = bcidal2 / T)) +
    geom_tile(aes(fill = 1 - wins)) +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(x = "bcidal1", y = "bcidal2", fill = "P(extinct)") +
    theme_minimal() +
    theme(
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm")
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
# plot with one independent variable and one dependent variable
ggplot(summary, aes(x = bcidal1 / T, y = wins)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    theme_light() +
    # scale_x_log10() +
    scale_y_continuous(limits = c(0, 1)) +
    # scale_y_log10(limits = c(1, 1e11)) +
    labs(
        x = "Mode (0 = bacteriostatic, 1 = bactericidal)",
        y = "Probability that resistance emerges in 100 hours"
    ) +
    theme(
        plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25)
    )


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