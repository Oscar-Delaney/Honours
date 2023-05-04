source("stochastic.R")

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
r <- 0.4
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

# read the file we just wrote back as summary
summary <- read.csv("results/summary_stats.csv")
# read the data file we just wrote back as data
load("results/data.RData")
# Plot the summary statistics
metrics_plot(summary, "R12")
metrics_plot(summary, "total")
log_plot(data[[5]], type = "all")

# Reconstruct Torella et al. 2010
metrics2 <- function(sols) {
    # Find the first time step when total population is below 1
    # noting that sols is in long dataframe format
    t_clear_S <- min(sols$time[sols$variable == "S" & sols$value < 1])
    t_clear_R1 <- min(sols$time[sols$variable == "R1" & sols$value < 1])
    t_clear_R2 <- min(sols$time[sols$variable == "R2" & sols$value < 1])
    t_clear <- max(t_clear_S, t_clear_R1, t_clear_R2)
    # Calculate the double-resistant population size at t_clear
    R12 <- mean(sols$value[sols$time == t_clear & sols$variable == "R12"])
    return(list(R12 = R12, t_clear = t_clear))
}

metrics_plot2 <- function(summary, metric, title, ylab) {
    ggplot(summary, aes(x = theta, y = 1/get(metric))) +
        geom_point(aes(color = factor(N0))) +
        geom_line(aes(color = factor(N0))) +
        theme_minimal() +
        labs(
            title = title,
            x = "Drug interaction parameter (theta)",
            y = ylab
        ) +
        # log-transform y axis
        scale_y_log10() +
        scale_color_discrete(name = "N0") +
        theme(
            legend.position = "bottom",
            plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size = 25),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20)
        )
}


m_rate <- 1e-6
theta = seq(-1, 1, by = 0.5)
N0 = 10 ^ seq(6, 10, by = 1)
summary <- expand.grid(theta = theta, N0 = N0)
summary$R12 <- 0
summary$t_clear <- 0

# Run the simulations
data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- sols <- simulate(
        rep = 1,
        deterministic = TRUE,
        stewardship = "comb",
        time = 50,
        freq = 1000,
        D = 1,
        N0 = summary$N0[i],
        HGT = 0,
        m1 = m_rate,
        m2 = m_rate,
        d1 = 0,
        d2 = 0,
        influx = 7 * c(A1 = 1, A2 = 1),
        init = 3.6e10 * c(S = 1 - 2 * m_rate, R1 = m_rate, R2 = m_rate, R12 = 0),
        psi = 0.94 * c(1, 1, 1, 1),
        zeta1 = c(1, 1e9, 1, 1e9),
        zeta2 = c(1, 1, 1e9, 1e9),
        theta = summary$theta[i] * c(1, 1, 1, 1),
        mu = 0.117 * c(1, 1, 1, 1),
        k = 1e9 * c(1, 1, 1, 1),
    )
}

# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    metric_results <- metrics2(data[[i]])
    summary$R12[i] <- metric_results$R12
    summary$t_clear[i] <- metric_results$t_clear
}

metrics_plot2(summary, "R12", "Public Health Efficacy", "1 / Double-resistant population size")
metrics_plot2(summary, "t_clear", "Clinical Efficacy", "1 / Time to clearance (days)")
# Save the simulation results to a file
save(data, file = "results/data2.RData")
# Save summary statistics to a CSV file
write.csv(summary, file = "results/summary_stats2.csv", row.names = FALSE)

sol = data[[3]] 
sol[sol$variable =="N",]
