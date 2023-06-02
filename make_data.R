source("stochastic.R")

# A function to interpret a simulation run
metrics <- function(sols, D) {
    # Filter the data frame to keep only the last time step
    final <- sols[sols$time == max(sols$time), ]
    # Calculate the average final double-resistant population size
    R12 <- mean(final$value[final$variable == "R12"])
    # Calculate the average final bacterial population size
    total <- mean(final$value) * 4
    # Calculate the proportion of runs with MDR emergingg=
    R12 <- sols[sols$variable == "R12" & sols$time == max(sols$time), ]$value
    wins <- mean(R12 * D > 1e3)
    return(wins)
}

metrics_plot <- function(summary, metric) {
    title <- ifelse(metric == "R12", "Final R12 Population Size",
        "Final Total Population Size")
    # Plot the summary statistics
    ggplot(summary, aes(x = D, y = HGT, size = log10(get(metric)+1))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal() +
    labs(title = title,
        x = "Dilution fraction D)",
        y = "Recombination rate",
        size = paste("log_10 of", metric)) +
    theme(
        legend.position = "bottom",
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )
}

log_plot(data[[72]][[1]][data[[1]][[1]]$rep <= 10, ])
sol <- data[[3]][[1]]
sol[sol$variable == "N" & sol$time == max(sol$time), ]$value
# Create a grid of parameter combinations
time <- 3e2
# summary <- expand.grid(HGT = HGT_range)
summary <- expand.grid(dose_rep = seq(1, 24, 3), kappa = seq(0.5, 3, 0.3))
# summary$wins <- 0

data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- simulate(
    init = c(S = 1e15, R1 = 0, R2 = 0, R12 = 0),
    N0 = 1e16,
    k = 0,
    alpha = 0,
    mu = 0.15 * c(1, 0.8, 0.8, 0.64),
    phi1 = 0.6,
    phi2 = 0.6,
    time = time,
    rep = 1e2,
    D = 1,
    HGT = 0,
    dose_rep = summary$dose_rep[i],
    dose_gap = 24 / summary$dose_rep[i],
    influx = 15 * c(A1 = 1, A2 = 1) / summary$dose_rep[i],
    m1 = 1.5e-8, m2 = 1.5e-8,
    kappa1 = summary$kappa[i],
    kappa2 = summary$kappa[i],
)
print(i / nrow(summary))
}

# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    sols <- data[[i]][[1]]
    R12 <- sols[sols$variable == "R12" & sols$time == max(sols$time), ]$value
    summary$wins[i] <- mean(R12 * data[[i]][[2]]$D > 1e3)
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
