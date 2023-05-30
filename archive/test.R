# adaptivetau investigation
library(ggplot2)
library(adaptivetau)

transitions <- list(
    s = c(y1 = +1) # stochastic
    # d = c(y2 = +1) # deterministic
)

rates <- function(state, parms, t) {
    return(c(
        s = state["y1"]
        # d = state["y2"]
        )
    )
}

c <- 1e3
state <- c * c(y1 = 1) #, y2 = 1)
sol <- as.data.frame(ssa.adaptivetau(state, transitions, rates,
    c, tf = 10,tl.params = list(epsilon = 1e-3)))
sol$y3 <- c * exp(sol$time)
# sol$y3 <- 1
# for (i in seq_len(nrow(sol))[-1]) {
#     sol$y3[i] <- sol$y3[i-1] * (1 + sol$time[i] - sol$time[i-1])
# }

ggplot(sol, aes(x = time)) +
  geom_step(aes(y = y1, color = "Stochastic")) +
#   geom_step(aes(y = y2, color = "Deterministic")) +
  geom_line(aes(y = y3, color = "exp(t)")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = 10^seq(0, 20),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_minimal() +
  labs(title = "Adaptive tau method investigation",
       x = "Time",
       y = "Population size",
       color = "method") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )


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
save(data, file = "C:/Users/s4528540/Downloads/results/data2.RData")
# Save summary statistics to a CSV file
write.csv(summary, file = "C:/Users/s4528540/Downloads/results/summary_stats2.csv", row.names = FALSE)

sol = data[[3]] 
sol[sol$variable =="N",]

sol = wahl2(rep = 1, tau = 0.10, D = 0.9)[[1]]
sol[sol$time >= 99, ]
