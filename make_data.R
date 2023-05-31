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
time <- 5e2
# freq_range <- seq(1, 30, by = 1)
# D_range <- 10 ^ (-0.1 * freq_range)
HGT_range <- 10 ^ seq(-20, -10, by = 5)
summary <- expand.grid(HGT = HGT_range)
# summary <- data.frame(dose_rep = 0:0)
summary$wins <- 0

# Run the simulations (case 1)
data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- 
print(i / nrow(summary))
}
log_plot(sol <- simulate(
    init = c(S = 3e7, R1 = 0, R2 = 0, R12 = 0),
    N0 = 1e12,
    k = 1e10,
    mu = 1.5 * c(1, 0.9, 0.9, 0.81),
    phi1 = 2,
    phi2 = 2,
    time = 60,
    rep = 10,
    D = 1,
    HGT = 0
)[[1]])

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

# add error bars
png("test.png", width = 12, height = 10, units = "in", res = 300)
# plot with one independent variable and one dependent variable
ggplot(summary, aes(x = dose_rep, y = wins)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    theme_light() +
    # scale_x_log10() +
    scale_y_continuous(limits = c(0,1)) +
    labs(
        title = "Cycling in strong-selection weak-mutation regime",
        x = "Doses per cycle (0 = combination therapy)",
        y = "Probability that MDR emerges in 1000 hours"
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
save(data, file = "C:/Users/s4528540/Downloads/results/data4.RData")
# read the data file we just wrote back as data
load("C:/Users/s4528540/Downloads/results/data4.RData")

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
time <- 30
init_S <- 1e6
df <- data.frame(max_step = 10 ^ seq(0, -4, by = -1), run_time = NA, r_error = NA)

for (i in seq_len(nrow(df))) {
    df$run_time[i] <- system.time({
        sol = simulate(
            init = c(S = init_S, R1=0, R2=0, R12=0),
            tau = 1e3,
            mu = 1,
            alpha = 0,
            k = 0,
            time = time,
            dose_rep = 1,
            rep = 1,
            influx = c(A1 = 0, A2 = 0),
            HGT = 0,
            m1 = 0,
            m2 = 0,
            max_step = df$max_step[i],
            deterministic = F
            )[[1]]
    })["elapsed"]
    df$r_error[i] <- 1 - log(sol[sol$time == time & sol$variable == "S",]$value / init_S) / time
    print(i/nrow(df))
}
df
