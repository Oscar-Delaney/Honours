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

sols <- data[[12]][[1]]
log_plot(sols[sols$rep <=10,])

# Create a grid of parameter combinations
summary <- expand.grid(dose_rep = seq(1, 24, 1), kappa = seq(1, 1, 0.1))

data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- simulate(
    init = c(S = 1e9, R1 = 0, R2 = 0, R12 = 0),
    N0 = 1e5,
    k = 1e6,
    supply = 1e8,
    mu = 1.6 * c(1, 0.9, 0.9, 0.8),
    phi1 = 2,
    phi2 = 2,
    time = 60,
    tau = 1e4,
    rep = 1e3,
    HGT = 0,
    dose_rep = summary$dose_rep[i],
    dose_gap = 12 / summary$dose_rep[i],
    influx = 10 * c(A1 = 1, A2 = 1) / summary$dose_rep[i],
    m1 = 1e-9, m2 = 1e-9,
    alpha = 1,
    kappa1 = summary$kappa[i],
    kappa2 = summary$kappa[i],
)
log_plot(data[[i]][[1]], use = c("S", "R1", "R2", "R12", "N"))
print(i / nrow(summary))
}

# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    summary$wins[i] <- mean(final_pop(data[[i]][[1]], c("R1", "R2")) > 1e3)
    rep <- data[[i]][[2]]$rep
    ci <- binom.test(summary$wins[i] * rep, rep, conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}

# plot with one independent variable and one dependent variable
ggplot(summary, aes(x = dose_rep, y = wins)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    theme_light() +
    # scale_x_log10() +
    scale_y_continuous(limits = c(0,1)) +
    labs(
        title = "Resistance Evolution",
        x = "Number of doses per day",
        y = "Probability that MDR emerges in 300 hours"
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

ggplot(summary, aes(x = dose_rep, y = kappa)) +
    geom_tile(aes(fill = wins)) +
    scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
    labs(x = "doses per day", y = "kappa", fill = "proportion",
            title = "Hill shape parameter interacts with dosing frequency",
            subtitle = "proportion of runs where MDR becomes established") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 29, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )


# # Save the simulation results to a file
save(data, file = "C:/Users/s4528540/Downloads/results/data4.RData")
# read the data file we just wrote back as data
load("C:/Users/s4528540/Downloads/results/data4.RData")
