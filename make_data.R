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

# Create a grid of parameter combinations
time <- 200
freq_range <- seq(1, 30, by = 1)
D_range <- 10 ^ (-0.1 * freq_range)
HGT_range <- 10 ^ seq(-20, -10, by = 0.5)
summary <- expand.grid(HGT = HGT_range)
summary <- expand.grid(dose_rep = 0:20)
summary$wins <- 0

# Run the simulations (case 1)
data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- simulate(
        init = c(S = 1e7, R1=0, R2=0, R12=0),
        N0 = 2e8,
        k = 1e7,
        mu = 1,
        time = 500,
        dose_rep = 10,
        rep = 1000,
        HGT = summary$HGT[i]
        )

}

log_plot(data[[11]][[1]][data[[1]][[1]]$rep >= 90, ], type = "all")

# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    summary$wins[i] <- metrics(data[[i]][[1]], data[[i]][[2]]$D)
}

# Run the simulations (case 2)
data <- list()
for (i in seq_len(nrow(summary))) {
    HGT <- summary$HGT[i]
    data[[i]] <- simulate(
        seed = NULL,
        init = 1e10 * c(S = 1, R1=0, R2=0, R12=0),
        N0 = 1e11,
        k = 1e10,
        mu = 0.6 * c(1, 1, 1, 0.5),
        time = 500,
        rep = 100,
        HGT = HGT
    )
}

# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    sols <- data[[i]][[1]]
    summary$final_pop[i] <- 4 * mean(sols$value[sols$time >= 400])
}

png("HGT_SSWM.png", width = 12, height = 10, units = "in", res = 300)
# plot with one independent variable and one dependent variable
ggplot(summary, aes(x = HGT, y = wins)) +
    geom_point(size = 3) +
    theme_light() +
    scale_x_log10() +
    labs(
        title = "HGT in strong-selection weak-mutation regime",
        x = "recombination rate",
        y = "Probability that MDR emerges in 500 hours"
    ) +
    theme(
        plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25)
    )
dev.off()
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
# rather than size to represent the probability of MDR emergence
# how could I make the color scale in the legend bigger?
# answer: use legend.key.size
# make the graph
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
save(data, file = "results/data4.RData")
# # read the data file we just wrote back as data
# load("results/data.RData")

# Plot the summary statistics
metrics_plot(summary, "R12")
metrics_plot(summary, "total")
log_plot(data[[97]][[1]], type = "all")


# Trying to find regions of paramter space where cycling is best
log_plot(simulate(
      rep = 1, dose_rep = 10, init = c(S = 1e8, R1=0, R2=0, R12=0),
      time = 300, D = 1, HGT = 0, freq = 10, N0 = 1e9, theta = 0,
      d1 = 0, d2 = 0, m1 = 5e-9, m2 = 5e-9, phi1 = 1, phi2 = 1,
      zeta1 = c(S = 1,8,4,8), zeta2 = c(S = 1,4,8,8),
      mu = c(1, 0.9, 0.9, 0.81), alpha = 0, k = 1e7, cycl = F)[[1]],
      type = "all")



# Weak mutation strong selection recombination is less useful
log_plot(simulate(
        init = c(S = 1e7, R1=0, R2=0, R12=0),
        N0 = 2e8,
        k = 1e7,
        mu = 1,
        time = 500,
        dose_rep = 10,
        rep = 4,
        HGT = 1e-10
        )[[1]])
