# Load the matrixStats package
library(matrixStats)

# Create a matrix from the list
value <- c(13.76145, 69.00000, 65.00000, 22.64049)
print(
    matrixStats::rowQuantiles(matrix(value, nrow = 1), probs = c(0.25))[[1]]
)
data_matrix <- matrix(data_list, nrow = 1)

# Calculate the 25th and 75th percentiles using matrixStats::rowQuantiles
percentiles <- matrixStats::rowQuantiles(data_matrix, probs = c(0.25))

# Access the 25th and 75th percentiles
quantile_25 <- percentiles[1]
quantile_75 <- percentiles[2]

# Print the 25th and 75th percentiles
print(quantile_25)
print(quantile_75)


# Unused code
IQR_bounds = list(colQuantiles(as.matrix(value), probs = c(0.25, 0.5, 0.75), na.rm = TRUE))

num_cores <- detectCores()
solutions <- tryCatch(
{
    mclapply(1:rep, function(x) {single_run(config, x)}, mc.cores = num_cores)
},
error = function(e) {
    lapply(1:rep, function(x) {single_run(config, x)})
}
)


  print(length(solution$time))
  duplicated_time_values <- solution$time[duplicated(solution$time)]
  print(duplicated_time_values)


# adaptivetau investigation

transitions <- list(
    s = c(y1 = +1),
    d = c(y2 = +1),
    u = c(y3 = +1)
)

rates <- function(state,parms,t) {
    return(c(
        s = state["y1"],
        d = state["y2"],
        u = 0
        )
    )
}

state <- c(y1 = 1e0, y2 = 1e0, y3=0)
parms <- 0
sol <- as.data.frame(ssa.adaptivetau(state, transitions, rates,
    parms, tf = 10, deterministic = c(F,T,F)))
# sol$y3 <- exp(sol$time)
ggplot(sol, aes(x = time)) +
    geom_step(aes(y = y1, color = "y1")) +
    geom_step(aes(y = y2, color = "y2")) +
    # geom_step(aes(y = y3, color = "y3")) +
    scale_color_manual(name = "Species", values = c("y1" = "blue", "y2" = "red")) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      breaks = 10^seq(0, 20),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_minimal() +
    labs(title = "Adaptive tau method",
        x = "Time",
        y = "Population size") +
    theme(
        legend.position = "bottom",
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )
# plot(sol$time,sol$y1)
sol


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
        # scale_y_log10() +
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
theta = seq(-2, 2, by = 0.1)
N0 = 10 ^ seq(8, 10, by = 1)
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
        tau = 1000,
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

data_to_plot <- summary[summary$N0 == 1e8,]
metrics_plot2(data_to_plot, "R12", "Public Health Efficacy", "1 / Double-resistant population size")
metrics_plot2(data_to_plot, "t_clear", "Clinical Efficacy", "1 / Time to clearance (days)")
# Save the simulation results to a file
save(data, file = "results/data2.RData")
# Save summary statistics to a CSV file
write.csv(summary, file = "results/summary_stats2.csv", row.names = FALSE)

sol = data[[3]] 
sol[sol$variable =="N",]


# Critique of Wahl papers
source("stochastic.R")
# tau = 0.1
log_plot(wahl2(rep = 1, tau = 0.10, D = 0.9)[[1]], type = "all")
log_plot(simulate(tau = 0.1, time = 300, mu = 1, dose_gap = 1e3, D = 0.9)[[1]])
# No resource constraints
wahl1 <- function(rep, tau, D) {
    sols <- simulate(
        seed = NULL,
        rep = rep,
        tau = tau,
        D = D,
        influx = c(A1 = 0, A2 = 0),
        N0 = 1,
        alpha = 0,
        k = 0,
        m1 = 1e-9,
        m2 = 0,
        mu = c(S = 1.025, R1 = 1.025 + s),
        init = c(S = round(3e8 * D), R1 = 0, R2 = 0, R12 = 0)
    )
}

# Constant resource concentration in dilution media
wahl2 <- function(rep, tau, D) {
    sols <- simulate(
        deterministic = TRUE,
        dt = 0.01,
        seed = NULL,
        time = 400,
        rep = rep,
        tau = tau,
        D = D,
        influx = c(A1 = 0, A2 = 0),
        N0 = 3e8,
        alpha = 1,
        k = 3e7,
        m1 = 1e-99,
        m2 = 0,
        mu = c(S = 1.15, R1 = 1.15 + s),
        init = c(S = round(3e8 * D), R1 = 0, R2 = 0, R12 = 0)
    )
}

# Constant total resource supply per time
wahl3 <- function(rep, tau, D) {
    sols <- simulate(
        seed = NULL,
        rep = rep,
        tau = tau,
        D = D,
        influx = c(A1 = 0, A2 = 0),
        N0 = 1e10 * tau / (1 - D),
        alpha = 1,
        k = 1e9,
        m1 = 1e-9,
        m2 = 0,
        mu = c(S = 1.2, R1 = 1.2 + s),
        init = c(S = round(1e10 * D), R1 = 0, R2 = 0, R12 = 0)
    )
}

log_plot(data[[1]][[1]][data[[1]][[1]]$rep >= 90, ], type = "all")


# Wahl 1
summary <- data.frame(tau = seq(0.1, 3, by = 0.1))
summary$D <- exp(-summary$tau)

# Wahl 2
summary <- expand.grid(tau = seq(0.1, 3, 0.2), D = exp(-seq(0.1, 3, 0.2)))

# General
s <- 0.1
summary$wins <- 0
data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- wahl2(rep = 100, tau = summary$tau[i], D = summary$D[i])
}

# calculate summary statistics
for (i in seq_len(nrow(summary))) {
    # count the number of runs on which R1 was at least 1e-5 of S at the end
    end <- data[[i]][[1]][data[[i]][[1]]$time == 100,]
    R1 <- end[end$variable == "R1",]$value
    S <- end[end$variable == "S",]$value
    summary$wins[i] <- mean(ifelse(S==0, R1, R1 / S) > 1e-7)
}


# use ggplot to scatterplot D on the x axis and wins on the y axis
ggplot(summary, aes(x = D, y = wins)) +
    geom_point() +
    theme_light() +
    labs(
        title = "Optimal Dilution Ratio (resource unconstrained)",
        x = "D",
        y = "Probability that a mutant fixes in 100 hours"
    ) +
    theme(
        plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25)
    )

ggplot(summary, aes(x = D, y = tau)) +
    geom_tile(aes(fill = wins)) +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(x = "D", y = "tau", fill = "proportion",
            title = "Optimal Dilution Ration",
            subtitle = "proportion of runs where R1 becomes established") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )
