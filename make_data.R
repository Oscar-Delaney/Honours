source("Stochastic.R")

# A function to interpret a simulation run
metrics <- function(sim_result) {
  # Filter the data frame to keep only the last time step
  final <- sim_result[sim_result$time == max(sim_result$time), ]
  # Calculate the average final double-resistant population size
  R12 <- mean(final$value[final$variable == "R12"])
  # Calculate the average final bacterial population size
  total <- mean(final$value) * 4
  return(list(R12 = R12, total = total))
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
time <- 100
r <- 0.5
bnecks <- seq(3, 10, by = 1)
D_range <- exp(-r * time / bnecks)
HGT_range <- 10 ^ seq(-20, -13, by = 1)
summary <- expand.grid(D = D_range, HGT = HGT_range)
summary$R12 <- 0
summary$total <- 0

# Run the simulations
data <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    HGT <- summary$HGT[i]
    data[[i]] <- simulate(rep = 10, time = time, D = D, HGT = HGT, freq = -log(D)/r)
}

# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    metric_results <- metrics(data[[i]])
    summary$R12[i] <- metric_results$R12
    summary$total[i] <- metric_results$total
}

# Save the simulation results to a file
save(data, file = "results/data.RData")
# Save summary statistics to a CSV file
write.csv(summary, file = "results/summary_stats.csv", row.names = FALSE)

# Plot the summary statistics
metrics_plot(summary, "R12")
metrics_plot(summary, "total")
