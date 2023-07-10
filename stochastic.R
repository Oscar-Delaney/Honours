library(adaptivetau)
library(deSolve)
library(tidyverse)
library(ggnewscale)
library(future)
library(future.apply)

# floating point compatible version of x%%y = 0
is.multiple <- function(denominator, numerator) {
  ratio <- numerator / denominator
  near_int <- round(ratio)
  return(near(ratio, near_int))
}

# a growth rate function for nutrient-limited growth
monod <- function(N, mu, k) {
  if (N == 0) {return(0)}
  return(mu * ifelse(k == 0, 1, 1 / (1 + k / N)))
}

# implement a population bottleneck and/or drug dosing
event <- function(state, config) {
  with(config, {
    t <- state["time"] # extract the time
    pops <- state[names(init)] # extract just the cell counts
    N <- state["N"] # extract the nutrient concentration
    drugs <- state[c("A1", "A2")] # extract the drug concentrations
    if (is.multiple(tau, t)) {
      if (deterministic) {
        pops <- pops * D
      } else {
        pops <- setNames(rbinom(length(pops), pops, D), names(pops))
        }
      N <- N * D + N0 * (1 - D)
      drugs <- drugs * D
    }
    if (is.multiple(dose_gap, t)) {
      drugs <- drugs * keep_old_drugs + influx * pattern
    }
    return(c(pops, N, drugs))
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
  N_antidepletion = c(N = +1),
  A1_depletion = c(A1 = -1),
  A2_depletion = c(A2 = -1)
))
}

# compute the rate at which each transition occurs
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # Calculate death rates per cell
    A1e <- 1 / (1 + (A1 / zeta1 * 2 ^ (i21 * (1 - 2 ^ -A2))) ^ -kappa1)
    A2e <- 1 / (1 + (A2 / zeta2 * 2 ^ (i12 * (1 - 2 ^ -A1))) ^ -kappa2)
    deaths <- bcidal1 * A1e + bcidal2 * A2e
    statics <- 1 - pmax(1, bstatic1 * A1e + bstatic2 * A2e)
    # Calculate replication rates
    replication_rates <- c(S, R1, R2, R12) * monod(N, mu, k) * statics
    if (near(min(statics), 0)) {print("statics out of bounds")}
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
    death_rates <- c(S, R1, R2, R12) * (deaths + delta)
    # Calculate other rates
    HGT_MDR_loss <- HGT * R12 * S
    HGT_MDR_gain <- HGT * R1 * R2
    N_depletion <- sum(replication_rates * alpha)
    N_antidepletion <- supply
    A1_depletion <- d1 * A1
    A2_depletion <- d2 * A2
    # Combine all rates and return
    rate_list <- c(growth_rates, death_rates, HGT_MDR_loss, HGT_MDR_gain,
      N_depletion, N_antidepletion, A1_depletion, A2_depletion)
    return(setNames(rate_list, names(make_transitions())))
  })
}

# a function to convert the transition rates into an ODE function
ode_rates <- function(t, state, config) {
  with(as.list(rates(state, config, t)), {
    return(list(c(
      S = S_growth - S_death + HGT_MDR_gain - HGT_MDR_loss,
      R1 = R1_growth - R1_death - HGT_MDR_gain + HGT_MDR_loss,
      R2 = R2_growth - R2_death - HGT_MDR_gain + HGT_MDR_loss,
      R12 = R12_growth - R12_death + HGT_MDR_gain - HGT_MDR_gain,
      N = -N_depletion,
      A1 = -A1_depletion,
      A2 = -A2_depletion
    )))
  })
}

# a function to implement one run of the model
single_run <- function(config, x) {
  with(config, {
    # Define the transitions of the model
    transitions <- make_transitions()
    # Initialise the state variables
    state <- c(init, N = N0, influx * pattern)
    time_grid <- seq(0, time, by = dt) # a common time grid for all runs
    event_times <- sort(unique(round(c(time, seq(0, time, tau),
      seq(0, time, dose_gap)), 10)))
    for (t in event_times[-length(event_times)]) {
      # Determine the time until the next bottleneck or dose
      end <- min(event_times[event_times > t] - t)
      # Run the model between events, deterministically or stochastically
      if (deterministic) {
        times <- c(time_grid[time_grid <= end], end) # ensures length(times) > 1
        new <- ode(state, times, ode_rates, config)
      } else {
        # set the seed for reproducibility
        if (is.numeric(seed)) set.seed(round(seed + (x * time + t) / tau))
        new <- ssa.adaptivetau(
          state, transitions, rates, config, tf = end,
          tl.params = list(maxtau = max_step),
          deterministic = grep("depletion", names(transitions))
        )
      }
      # Make the time column reflect the overall time accurately
      new[, "time"] <- new[, "time"] + t
      # Avoid duplicate times by adding a small increment after the bottleneck
      new[1, "time"] <- new[1, "time"] * (1 + 1e-6)
      # Update the solution
      solution <- if (t == 0) new else rbind(solution, new)
      # flip the pattern after the appropriate number of infusions of the drug
      if (cycl && is.multiple(dose_gap * dose_rep, t + end)) {
        config$pattern <- 1 - config$pattern
      }
      # Run the bottleneck and/or dose and update the state
      state <- event(new[nrow(new), ], config)
    }
    # Interpolate the solution to the common time grid
    approx_vars <- lapply(colnames(solution), function(var) {
      approx(solution[, "time"], solution[, var], xout = time_grid)$y
    })
    solution_interpolated <- data.frame(
      setNames(approx_vars, colnames(solution)),
      rep = x
    )
    # return(solution)
    return(solution_interpolated)
  })
}

# a function to simulate the model
simulate <- function(
  rep = 1, # number of runs of the simulation
  seed = NULL, # seed for reproducibility
  deterministic = FALSE, # should be either TRUE or FALSE
  cycl = TRUE, # should be either TRUE or FALSE
  keep_old_drugs = TRUE, # should be either TRUE or FALSE
  time = 100, # time to simulate, in hours
  dt = 0.1, # time step, in hours
  max_step = Inf, # SSA max step parameter, only used if deterministic = FALSE
  tau = 10, # frequency of bottlenecks, in hours
  dose_rep = 1, # number of doses of the same drug before switching
  dose_gap = 10, # time between doses of the same drug, in hours
  D = 0.1, # dilution ratio at bottlenecks
  N0 = 1e15, # initial nutrient concentration
  HGT = 0, # rate of horizontal gene transfer
  m1 = 1e-9, # rate of mutations conferring resistance to drug 1
  m2 = 1e-9, # rate of mutations conferring resistance to drug 2
  d1 = log(2) / 3.5, # rate of drug 1 elimination
  d2 = log(2) / 3.5, # rate of drug 2 elimination
  i12 = 0, # interaction effect of drug 1 on drug 2
  i21 = 0, # interaction effect of drug 2 on drug 1
  bcidal1 = 0.6, # maximum drug 1 death rate
  bcidal2 = 0.6, # maximum drug 2 death rate
  bstatic1 = 0, # maximum reduction in growth rate due to drug 1
  bstatic2 = 0, # maximum reduction in growth rate due to drug 2
  influx = 7 * c(A1 = 1, A2 = 1), # drug influx concentrations, units of zeta_s
  init = c(S = 1e12, R1 = 0, R2 = 0, R12 = 0), # initial population sizes
  zeta1 = c(S = 1, R1 = 28, R2 = 1, R12 = 28), # [drug 1] at half-max effect
  zeta2 = c(S = 1, R1 = 1, R2 = 28, R12 = 28), # [drug 2] at half-max effect
  kappa1 = c(S = 1, R1 = 1, R2 = 1, R12 = 1), # Hill function shape parameter
  kappa2 = c(S = 1, R1 = 1, R2 = 1, R12 = 1), # Hill function shape parameter
  mu = 0.88 * c(S = 1, R1 = 0.9, R2 = 0.9, R12 = 0.81), # maximum growth rate
  delta = 0 * c(S = 1, R1 = 1, R2 = 1, R12 = 1), # intrinsic death rate
  k = 1e14 * c(S = 1, R1 = 1, R2 = 1, R12 = 1), # [N] at half-max growth rate
  alpha = c(S = 1, R1 = 1, R2 = 1, R12 = 1), # nutrients used per replication
  supply = 0 # nutrient supply rate
  ) {
  # Define the parameters of the model
  config <- as.list(environment())
  config$influx <- influx * c(zeta1["S"], zeta2["S"]) # normalise drug units
  config$pattern <- if (cycl) c(1, 0) else c(1, 1) # pattern of drug application
  # return(config)
  # Run the simulation rep number of times, using parallelisation if possible
  plan(multisession) # compatible with both unix and Windows
  solutions <- bind_rows(future_lapply(1:rep, function(x) {
    single_run(config, x)},
    future.seed = TRUE))
  # Convert the solutions to long format
  long <- pivot_longer(solutions, cols = -c(time, rep), names_to = "variable")
  return(list(long, config))
}

log_plot <- function(solutions, type = "all", use = c("S", "R1", "R2", "R12")) {
  # Choose the type of central tendency and range to plot
  if (type == "mean") {
    summary <- solutions %>%
      group_by(time, variable) %>%
      reframe(
        central = mean(value),
        se = sd(value) / sqrt(n()),
        lower = max(0, central - 1.96 * se),
        upper = central + 1.96 * se,
        rep = 1
      )
  } else if (type == "median") {
    summary <- solutions %>%
      group_by(time, variable) %>%
      reframe(
        IQR_bounds = list(quantile(value, c(0.25, 0.5, 0.75))),
        rep = 1
      ) %>%
      # convert the list of quantiles to individual columns
      mutate(
        lower = map_dbl(IQR_bounds, 1),
        central = map_dbl(IQR_bounds, 2),
        upper = map_dbl(IQR_bounds, 3)
      ) %>%
      select(-IQR_bounds)
  } else if (type == "all") {
    summary <- solutions
    colnames(summary)[colnames(summary) == "value"] <- "central"
  } else {
    stop("type must be 'all', 'median' or 'mean'")
  }
  filtered <- summary %>%
    filter(variable %in% use) %>%
    mutate(variable = factor(variable, levels = use)) %>%
    arrange(variable)
  # Initialise the colours
  colors <- c("black", "navy", "#800000", "#008000", "gold")
  # Create antibiotic concentrations data frame
  background_df <- solutions %>%
    filter(rep == min(rep) & variable %in% c("A1", "A2")) %>%
    group_by(variable) %>%
    mutate(value = value / max(value)) %>%
    mutate(across(value, ~replace_na(.x, 0))) %>%
    ungroup() %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(xmin = lag(time), xmax = time) %>%
    filter(!is.na(xmin)) %>%
    select(xmin, xmax, A1, A2)
  peak <- max(solutions[solutions$variable %in% use, "value"], na.rm = TRUE)
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
    geom_line(data = filtered, aes(x = time, y = central, color = variable,
      linetype = factor(rep)), linewidth = 1 + 0.5 / max(filtered$rep)) +
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
  if (type %in% c("mean", "median") && max(solutions$rep) > 1) {
    plot <- plot +
      geom_ribbon(data = filtered, alpha = 0.3,
        aes(x = time, ymin = lower, ymax = upper, fill = variable)) +
      scale_fill_manual(values = colors)
  }
  # Display the plot
  print(plot)
}

system.time(log_plot(simulate()[[1]]))
# simulate(deterministic = TRUE)[[1]]$value[7001:7007]
# [1] 3.051810e+14 2.515378e+11 1.181122e+12 1.717214e+09 6.784406e+12
# [6] 1.333550e-02 9.662624e-01