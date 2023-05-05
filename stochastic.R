library(adaptivetau)
library(deSolve)
library(tidyverse)
library(ggnewscale)
library(future)
library(future.apply)

# a pharmacodynamic function for singele-antibiotic induced killing of bacteria
hill <- function(A, params) {
  # code to extract the individual parameters
  with(as.list(sapply(c("psi", "phi", "zeta", "kappa"), function(x) {
    x <- unname(params[grep(x, names(params))])
  })), {
    # compute the actual hill function
    return(phi * (A / zeta)^kappa / ((A / zeta)^kappa - (psi - phi) / psi))
  })
}

# interaction between the two drugs
interaction <- function(A1_death, A2_death, phi1, phi2, theta) {
  return(theta * A1_death * A2_death / sqrt(phi1 * phi2))
}

# total death rate from two antibiotics with interactions
death <- function(A1, A2, params) {
  death1 <- hill(A1, params[c("psi", "phi1", "zeta1", "kappa1")])
  death2 <- hill(A2, params[c("psi", "phi2", "zeta2", "kappa2")])
  interaction <- interaction(death1, death2,
    params["phi1"], params["phi2"], params["theta"])
  return(death1 + death2 + interaction)
}

# a growth rate function for nutrient-limited growth
monod <- function(N, params) {
  with(as.list(params), {
    return(mu * N / (N + k))
  })
}

# a function that reduces all populations by a factor of D, in expectation
bottleneck <- function(state, config) {
  with(config, {
    pops <- state[rownames(params)]
    if (deterministic) {
      diluted <- pops * D
    } else {
      diluted <- setNames(rbinom(length(pops), pops, D), names(pops))
    }
    new_state <- c(diluted, state["N"] * D + N0 * (1 - D),
    state[c("A1", "A2")] * D + influx * pattern * (1 - D))
    return(new_state)
  })
}

# a function outputting the transitions that can occur in the model
make_transitions <- function() {
  return(list(
  S_growth = c(S = +1),
  R1_growth = c(R1 = +1),
  R2_growth = c(R2 = +1),
  R12_growth = c(R12 = +1),
  S_death = c(S = -1),
  R1_death = c(R1 = -1),
  R2_death = c(R2 = -1),
  R12_death = c(R12 = -1),
  HGT_MDR_loss = c(S = -1, R1 = +1, R2 = +1, R12 = -1),
  HGT_MDR_gain = c(S = +1, R1 = -1, R2 = -1, R12 = +1),
  N_depletion = c(N = -1),
  A1_depletion = c(A1 = -1),
  A2_depletion = c(A2 = -1)
))
}

# compute the rate at which each transition occurs
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # Calculate replication rates
    replication_rates <- sapply(rownames(params), function(var) {
      get(var) * monod(N, params[var, c("mu", "k")])
    })
    # chance of a replication in row i resulting in a strain j cell
    mutation <- matrix(c(
      S  = c((1 - m1) * (1 - m2), m1 * (1 - m2), (1 - m1) * m2, m1 * m2),
      R1 = c(0, 1 - m2, 0, m2),
      R2 = c(0, 0, 1 - m1, m1),
      R12 = c(0, 0, 0, 1)
    ), nrow = 4, byrow = TRUE)
    # Calculate growth rates including mutations
    growth_rates <- replication_rates %*% mutation
    # Calculate death rates
    death_rates <- sapply(rownames(params), function(var) {
      get(var) * death(A1 = A1, A2 = A2, params[var, ])
    })
    # Calculate other rates
    HGT_MDR_loss <- HGT * R12 * S
    HGT_MDR_gain <- HGT * R1 * R2
    N_depletion <- replication_rates %*% params[, "alpha"]
    A1_depletion <- d1 * A1
    A2_depletion <- d2 * A2
    # Combine all rates and return
    rate_list <- c(growth_rates, death_rates, HGT_MDR_loss, HGT_MDR_gain,
      N_depletion, A1_depletion, A2_depletion)
    return(setNames(rate_list, names(make_transitions())))
  })
}

# a function to convert the transition rates into an ODE function
ode_rates <- function(t, state, config) {
  adaptivetau_rates = rates(state, config, t)
  ode_rates_xyz <- with(as.list(adaptivetau_rates),{
    c(
      S = S_growth - S_death + HGT_MDR_gain - HGT_MDR_loss,
      R1 = R1_growth - R1_death - HGT_MDR_gain + HGT_MDR_loss,
      R2 = R2_growth - R2_death - HGT_MDR_gain + HGT_MDR_loss,
      R12 = R12_growth - R12_death + HGT_MDR_gain - HGT_MDR_gain,
      N = -N_depletion,
      A1 = -A1_depletion,
      A2 = -A2_depletion
    )
  })
  return(list(ode_rates_xyz))
}

# a function to implement one run of the model
single_run <- function(config, x) {
  with(config, {
    # Define the transitions of the model
    transitions <- make_transitions()
    # Initialise the state variables
    state <- c(init, N = N0, influx * pattern)
    t <- 0
    while (t < max(time_grid)) {
      # Run the model between bottlenecks, deterministically or stochastically
      if (deterministic) {
        times <- time_grid[time_grid <= freq]
        new <- ode(state, times, ode_rates, config)
      } else {
        new <- ssa.adaptivetau(
          state, transitions, rates, config, tf = freq,
          deterministic = grep("deplet", names(transitions))
        )
      }
      # Make the time column reflect the overall time accurately
      new[, "time"] <- new[, "time"] + t
      # Update the solution
      solution <- if (t == 0) new else rbind(solution, new[-1, ])
      # Update the time
      t <- t + freq
      # Run the bottleneck and update the state
      if (stewardship == "cycl" && (t / freq) %% dose_rep == 0) {
        config$pattern <- 1 - config$pattern
      }
      state <- bottleneck(new[nrow(new), ], config)
    }
    # Interpolate the solution to the common time grid
    variables <- c("S", "R1", "R2", "R12", "N", "A1", "A2")
    approx_vars <- lapply(variables, function(var) {
      approx(solution[, "time"], solution[, var], xout = time_grid)$y
    })
    solution_interpolated <- data.frame(
      time = time_grid,
      setNames(approx_vars, variables),
      rep = x
    )
    return(solution_interpolated)
  })
}

# a function to simulate the model
simulate <- function(
  rep = 1,
  dose_rep = 1,
  deterministic = FALSE, # should be either TRUE or FALSE
  stewardship = "cycl", #  "cycl" or "comb" or "1_only" or "2_only"
  time = 100, # time to simulate, in hours
  dt = 0.1, # time step, in hours
  freq = 10, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  N0 = 1e15, # initial nutrient concentration
  HGT = 0, # rate of horizontal gene transfer
  m1 = 1e-9, # rate of mutations conferring resistance to drug 1
  m2 = 1e-9, # rate of mutations conferring resistance to drug 2
  d1 = log(2) / 3.5, # rate of drug 1 elimination
  d2 = log(2) / 3.5, # rate of drug 2 elimination
  influx = 7 * c(A1 = 1, A2 = 1), # drug influx concentrations
  # lists of genotype-specific parameters, in the order S, R1, R2, R12
  init = c(S = 1e12, R1 = 0, R2 = 0, R12 = 0), # initial population sizes
  psi = 0.3 * c(1, 1, 1, 1), # growth rate with no drugs
  phi1 = 2 * psi, # maximum reduction in fitness for drug 1
  phi2 = 2 * psi, # maximum reduction in fitness for drug 2
  zeta1 = c(1, 28, 1, 28), # MIC drug 1, units scales so zeta1_S = 1
  zeta2 = c(1, 1, 28, 28), # MIC drug 2, units scales so zeta2_S = 1
  kappa1 = c(1, 1, 1, 1), # Hill function steepness parameter drug 1
  kappa2 = c(1, 1, 1, 1), # Hill function steepness parameter drug 2
  theta = c(0, 0, 0, 0), # drug interaction term
  mu = 0.88 * c(1, 0.9, 0.9, 0.81), # growth rate with infinite resources
  k = 1e14 * c(1, 1, 1, 1), # resource concentration at half-maximal growth
  alpha = c(1, 1, 1, 1) # resources used per unit growth
  ) {
  # Define the parameters of the model
  config <- list(
    D = D,
    N0 = N0,
    init = init,
    HGT = HGT,
    m1 = m1,
    m2 = m2,
    d1 = d1,
    d2 = d2,
    influx = influx,
    stewardship = stewardship,
    deterministic = deterministic,
    dose_rep = dose_rep,
    freq = freq,
    time_grid = seq(0, time, by = dt) # a common time grid for all runs
  )
  if (stewardship == "2_only") {
    config$pattern <- c(0, 1)
  } else if (stewardship == "comb") {
    config$pattern <- c(1, 1)
  } else {
    config$pattern <- c(1, 0)
  }
  config$params <- cbind(psi, phi1, phi2, zeta1, zeta2,
    kappa1, kappa2, theta, mu, k, alpha)
  rownames(config$params) <- c("S", "R1", "R2", "R12")
  # Run the simulation rep number of times, using parallelisation if possible
  plan(multisession) # compatible with both unix and Windows
  solutions <- bind_rows(future_lapply(1:rep, function(x) {
    single_run(config, x)},
    future.seed = TRUE))
  return(pivot_longer(solutions, cols = -c(time, rep), names_to = "variable"))
}

# A function to summarise the output of the simulation
summarise <- function(solutions) {
  summary <- solutions %>%
    group_by(time, variable) %>%
    reframe(
      mean = mean(value),
      sd = sd(value),
      se = sd / sqrt(n()),
      ci_lower = max(0, mean - 1.96 * se),
      ci_upper = mean + 1.96 * se,
      IQR_bounds = list(quantile(value, c(0.25, 0.5, 0.75)))
    ) %>%
    # convert the list of quantiles to individual columns
    mutate(
      IQR_lower = map_dbl(IQR_bounds, 1),
      median = map_dbl(IQR_bounds, 2),
      IQR_upper = map_dbl(IQR_bounds, 3)
    ) %>%
    select(-IQR_bounds) # Remove the IQR_bounds column
  return(summary)
}

log_plot <- function(solutions, type = "mean") {
  # Choose the type of central tendency and range to plot
  if (type == "median") {
    summary <- summarise(solutions)
    summary$central <- summary$median
    summary$lower <- summary$IQR_lower
    summary$upper <- summary$IQR_upper
    summary$rep <- 1
  } else if (type == "mean") {
    summary <- summarise(solutions)
    summary$central <- summary$mean
    summary$lower <- summary$ci_lower
    summary$upper <- summary$ci_upper
    summary$rep <- 1
  } else if (type == "all") {
    summary <- solutions
    summary$central <- summary$value
  } else {
    stop("type must be either 'median' or 'mean'")
  }
  # filtered <- filter(summary, variable %in% c("S", "R1", "R2", "R12"))
  filtered <- summary %>%
    filter(variable %in% c("S", "R1", "R2", "R12")) %>%
    mutate(variable = factor(variable, levels = c("S", "R1", "R2", "R12"))) %>%
    arrange(variable)
  # Initialise the colours
  colors <- c("black", "navy", "#800000", "#008000")
  # Create antibiotic concentrations data frame
  background_df <- solutions %>%
    filter(rep == 1 & variable %in% c("A1", "A2")) %>%
    group_by(variable) %>%
    mutate(value = value / max(value)) %>%
    ungroup() %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(xmin = lag(time), xmax = time) %>%
    filter(!is.na(xmin)) %>%
    select(xmin, xmax, A1, A2)
  peak <- max(solutions$value, na.rm = TRUE)
  # Create the plot
  plot <- ggplot() +
    # Add the gradient backgrounds
    geom_rect(data = background_df,
      aes(xmin = xmin, xmax = xmax, ymin = peak * 10^0.2,
         ymax = peak * 10^0.4, fill = A1), color = NA, alpha = 1) +
    scale_fill_gradient(low = "white", high = colors[2],
      limits = c(0, 1), name = "A1", labels = NULL) +
    new_scale_fill() +
    geom_rect(data = background_df,
      aes(xmin = xmin, xmax = xmax, ymin = peak * 10^0.4,
        ymax = peak * 10^0.6, fill = A2), color = NA, alpha = 1) +
    scale_fill_gradient(low = "white", high = colors[3],
      limits = c(0, 1), name = "A2", labels = NULL) +
    # Add the lines
    new_scale_fill() +
    geom_line(data = filtered, aes(x = time, y = central, color = variable, linetype = factor(rep)),
      linewidth = 1 + 0.5 / max(filtered$rep)) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = rep("solid", max(filtered$rep))) +
    guides(linetype = "none") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      breaks = 10^seq(0, 20),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(
      title = "Bacterial growth over time",
      x = "Time (hours)",
      y = "Population Size",
      color = "Strain",
      fill = "Strain"
    ) +
    theme_light() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
  # Add the confidence intervals
  if (type %in% c("mean", "median")) {
    plot <- plot +
      geom_ribbon(data = filtered, alpha = 0.3,
        aes(x = time, ymin = lower, ymax = upper, fill = variable)) +
      scale_fill_manual(values = colors)
  }
  # Display the plot
  print(plot)
}

system.time(log_plot(simulate(dose_rep = 7, deterministic = F, time = 100, rep = 1), type = "all"))
system.time(simulate(dose_rep = 7, deterministic = T, time = 100, rep = 1000))
sol$value[7001:7007]
# old
# > sol$value[7001:7007]
# [1] 2.008323e+14 7.518872e+13 1.713794e+10 7.624686e+09 1.822959e+12
# [6] 2.317572e-06 8.816358e-01
# new
# > sol$value[7001:7007]
# [1] 2.008323e+14 7.518872e+13 1.713794e+10 7.624686e+09 1.822959e+12
# [6] 2.317572e-06 8.816358e-01