# Replicating Nyhoegen and Uecker 2023
g_max <- 0.088
g_min <- -6.5
death_intrinsic <- 0.01
patient <- T
log_plot(simulate(
        init = c(S = 1e10, R1=0, R2=0, R12=0),
        mu = (g_max - death_intrinsic) * c(1, 0.9, 0.9, 0.81),
        phi1 = g_max - g_min,
        phi2 = g_max - g_min,
        dose_gap = 12,
        dose_rep = 1,
        tau = 1e5,
        influx = 5 * g_max / -g_min * c(A1 = 1, A2 = 1),
        theta = 0,
        zeta1 = c(S = 1, 28, 1, 28),
        zeta2 = c(S = 1, 1, 28, 28),
        cycl = T,
        time = 100,
        d1 = 0.05 * patient,
        d2 = 0.05 * patient,
        keep_old_drugs = patient,
        rep = 10
    )[[1]], type = "all")

# for a given dose_rep, find the influx that causes all strains to go extinct
# within 1000 hours

check <- function(dose_rep, influx) {
    g_max <- 0.088
    g_min <- -6.5
    death_intrinsic <- 0.01
    patient <- FALSE
    sols <- simulate(
        init = c(S = 1e10, R1=0, R2=0, R12=0),
        mu = (g_max - death_intrinsic) * c(1, 0.9, 0.9, 0.81),
        phi1 = g_max - g_min,
        phi2 = g_max - g_min,
        dose_gap = 12,
        dose_rep = dose_rep,
        tau = 1e5,
        influx = influx * g_max / -g_min * c(A1 = 1, A2 = 1),
        time = 500,
        d1 = log(2) / 3.5 * patient,
        d2 = log(2) / 3.5 * patient,
        keep_old_drugs = patient,
        rep = 10
    )[[1]]
    final <- sols[sols$time == max(sols$time), ]
    return(sum(final[final$variable %in% c("R1", "R2", "R12"), "value"]) == 0)
}

binary_search <- function(dose_rep, min_influx, max_influx) {
  while (max_influx - min_influx > 0.01) {  # adjust the precision as needed
    mid_influx <- (min_influx + max_influx) / 2
    if (check(dose_rep, mid_influx)) {
      max_influx <- mid_influx
    } else {
      min_influx <- mid_influx
    }
  }
  return(min_influx)
}
# Create sequences
dose_rep_seq <- seq(1, 30, 1)

# Initialize a result vector
summary <- data.frame(dose_rep = dose_rep_seq)
summary$min_dose <- NA

# Iterate over the possible values
for (i in dose_rep_seq) {
    summary$min_dose[i] <- binary_search(i, 1, 10)
}
summary

ggplot(summary, aes(x = dose_rep / 2, y = min_dose)) +
    geom_point(size = 5) +
    scale_y_log10() +
    theme_minimal() +
    labs(
        title = "Nyhoegen Figure 2a blue replicate",
        x = "Days between drug switches",
        y = "First effective dose"
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

clear_time <- function(sols) {
  sols %>%
    # Filter to only include the populations of interest
    filter(variable %in% c("S", "R1", "R2", "R12")) %>%
    # Group by time and rep, and calculate total population at each time and rep
    group_by(time, rep) %>%
    summarise(total_population = sum(value), .groups = 'drop') %>%
    # Arrange by rep and time
    arrange(rep, time) %>%
    # Group by rep
    group_by(rep) %>%
    # Add a column for the first time when total population is 0 for each rep
    mutate(first_zero_time = ifelse(any(total_population == 0), min(time[total_population == 0]), NA)) %>%
    # For each replicate, get the first row (which will have the first_zero_time)
    summarise(first_zero_time = first(first_zero_time), .groups = 'drop') %>%
    # Get the maximum first_zero_time across all reps
    summarise(max_first_zero_time = max(first_zero_time)) %>%
    # Pull out the maximum first zero time
    pull(max_first_zero_time)
}

# Create a grid of parameter combinations
influx_range <- seq(5, 20, by = 1)
dose_rep_range <- c(1, 2, 4, 8, 16)
kappa_range <- c(1, 2)
patient_range <- c(TRUE, FALSE)
cs_range <- c(1, 12 / 17, 7 / 17)
interaction_range <- c(-1, 0, 1)
summary <- expand.grid(influx = influx_range, cs = cs_range, dose_rep = dose_rep_range)
summary$t_clear <- 0

# Run the simulations
data <- list()
for (i in seq_len(nrow(summary))) {
    print(i)
    data[[i]] <- nyhoegen(
        dose_rep = summary$dose_rep[i],
        influx = summary$influx[i],
        kappa = 1,
        patient = TRUE,
        cs = summary$cs[i],
        interaction = 0,
        rep = 10
    )
}

# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    summary$t_clear[i] <- clear_time(data[[i]][[1]])
}
summary

plot_data <- summary %>%
    filter(dose_rep >= 1)
ggplot(plot_data, aes(x = influx, y = t_clear,
    color = as.factor(round(cs, 2)),
    linetype = as.factor(round(dose_rep, 2)))) +
    geom_line(size = 3) +
    theme_minimal() +
    labs(
        title = "Nyhoegen Figure 4 replication",
        x = "Drug dose in multiples of MIC",
        y = "Time to clearance",
        color = "Collateral sensitivity",
        linetype = "dose_rep"
    ) +
    theme(
        legend.position = "bottom",
        legend.box = "vertical",
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.size = unit(2, "cm")
    )


nyhoegen <- function(rep, dose_rep, influx, kappa, patient, cs, interaction) {
    g_max <- 0.088
    g_min <- -6.5
    death_intrinsic <- 0.01
    sols <- simulate(
        init = c(S = 1e10, R1=0, R2=0, R12=0),
        mu = (g_max - death_intrinsic) * c(1, 0.9, 0.9, 0.81),
        phi1 = g_max - g_min,
        phi2 = g_max - g_min,
        dose_gap = 12,
        dose_rep = dose_rep,
        tau = 1e5,
        kappa1 = kappa,
        kappa2 = kappa,
        zeta1 = c(S = 1, R1 = 28, R2 = 1 * cs, R12 = 28),
        zeta2 = c(S = 1, R1 = 1 * cs, R2 = 28, R12 = 28),
        influx = influx * (g_max / -g_min) ^ (1 / kappa) * c(A1 = 1, A2 = 1),
        time = 500,
        d1 = log(2) / 3.5 * patient,
        d2 = log(2) / 3.5 * patient,
        keep_old_drugs = patient,
        rep = rep,
        theta = interaction,
    )
    return(sols)
}
