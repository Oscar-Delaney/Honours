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

# a growth rate function for resource-limited growth
monod <- function(R, mu, k) {
  if (R == 0) {return(0)}
  return(mu * ifelse(k == 0, 1, 1 / (1 + k / R)))
}

# implement a population bottleneck and/or drug dosing
event <- function(state, config) {
  with(config, {
    t <- state["time"] # extract the time
    pops <- state[names(init)] # extract just the cell counts
    R <- state["R"] # extract the resource concentration
    drugs <- state[c("C_A", "C_B")] # extract the drug concentrations
    if (is.multiple(tau, t)) {
      if (deterministic) {
        pops <- pops * D
      } else {
        pops <- setNames(rbinom(length(pops), pops, D), names(pops))
        }
      R <- R * D + R0 * (1 - D) + N_add
      drugs <- drugs * D
    }
    if (is.multiple(dose_gap, t)) {
      drugs <- drugs * keep_old_drugs + influx * pattern
    }
    return(c(pops, R, drugs))
  })
}

# a function outputting the transitions that can occur in the model
make_transitions <- function() {
  return(list(
  S_growth = c(N_S = +1),
  N_A_growth = c(N_A = +1),
  N_B_growth = c(N_B = +1),
  N_AB_growth = c(N_AB = +1),
  S_death = c(N_S = -1),
  N_A_death = c(N_A = -1),
  N_B_death = c(N_B = -1),
  N_AB_death = c(N_AB = -1),
  HGT_MDR_loss = c(N_S = -1, N_A = +1, N_B = +1, N_AB = -1),
  HGT_MDR_gain = c(N_S = +1, N_A = -1, N_B = -1, N_AB = +1),
  N_depletion = c(R = -1),
  N_antidepletion = c(R = +1),
  A_depletion = c(C_A = -1),
  B_depletion = c(C_B = -1)
))
}

# compute the rate at which each transition occurs
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # Calculate death rates per cell
    Ae <- 1 / (1 + (C_A / zeta_A * 2 ^ (i_B_A * (1 - 2 ^ -C_B))) ^ -kappa_A)
    Be <- 1 / (1 + (C_B / zeta_B * 2 ^ (i_A_B * (1 - 2 ^ -C_A))) ^ -kappa_B)
    deaths <- bcidal_A * Ae + bcidal_B * Be
    statics <- 1 - pmin(1, bstatic_A * Ae + bstatic_B * Be)
    statics <- (1 - bstatic_A * Ae) * (1 - bstatic_B * Be)
    # Calculate replication rates
    replication_rates <- c(N_S, N_A, N_B, N_AB) * monod(R, mu, k) * statics
    # if (near(min(statics), 0)) {print("statics out of bounds")}
    # chance of a replication in row i resulting in a strain j cell
    mutation <- matrix(c(
      N_S  = c((1 - m_A) * (1 - m_B), m_A * (1 - m_B), (1 - m_A) * m_B, m_A * m_B),
      N_A = c(0, 1 - m_B, 0, m_B),
      N_B = c(0, 0, 1 - m_A, m_A),
      N_AB = c(0, 0, 0, 1)
    ), nrow = 4, byrow = TRUE)
    # Calculate growth rates including mutations
    growth_rates <- replication_rates %*% mutation
    # Calculate death rates
    death_rates <- c(N_S, N_A, N_B, N_AB) * (deaths + delta)
    # Calculate other rates
    HGT_MDR_loss <- HGT * N_AB * N_S
    HGT_MDR_gain <- HGT * N_A * N_B
    N_depletion <- sum(replication_rates * alpha)
    N_antidepletion <- supply
    A_depletion <- d_A * C_A
    B_depletion <- d_B * C_B
    # Combine all rates and return
    rate_list <- c(growth_rates, death_rates, HGT_MDR_loss, HGT_MDR_gain,
      N_depletion, N_antidepletion, A_depletion, B_depletion)
    return(setNames(rate_list, names(make_transitions())))
  })
}

# a function to convert the transition rates into an ODE function
ode_rates <- function(t, state, config) {
  with(as.list(rates(state, config, t)), {
    return(list(c(
      N_S = S_growth - S_death + HGT_MDR_gain - HGT_MDR_loss,
      N_A = N_A_growth - N_A_death - HGT_MDR_gain + HGT_MDR_loss,
      N_B = N_B_growth - N_B_death - HGT_MDR_gain + HGT_MDR_loss,
      N_AB = N_AB_growth - N_AB_death + HGT_MDR_gain - HGT_MDR_gain,
      R = -N_depletion,
      C_A = -A_depletion,
      C_B = -B_depletion
    )))
  })
}

# a function to implement one run of the model
single_run <- function(config, x) {
  with(config, {
    # Define the transitions of the model
    transitions <- make_transitions()
    # Initialise the state variables
    state <- c(init, R = R0, influx * pattern)
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
  R0 = 1e15, # initial resource concentration
  HGT = 0, # rate of horizontal gene transfer
  m_A = 1e-9, # rate of mutations conferring resistance to drug A
  m_B = 1e-9, # rate of mutations conferring resistance to drug B
  d_A = log(2) / 3.5, # rate of drug A elimination
  d_B = log(2) / 3.5, # rate of drug B elimination
  i_A_B = 0, # interaction effect of drug A on drug B
  i_B_A = 0, # interaction effect of drug B on drug A
  bcidal_A = 0.6, # maximum drug A death rate
  bcidal_B = 0.6, # maximum drug B death rate
  bstatic_A = 0, # maximum reduction in growth rate due to drug A
  bstatic_B = 0, # maximum reduction in growth rate due to drug B
  influx = 7 * c(C_A = 1, C_B = 1), # drug influx concentrations, units of zeta_s
  init = c(N_S = 1e12, N_A = 0, N_B = 0, N_AB = 0), # initial population sizes
  zeta_A = c(N_S = 1, N_A = 28, N_B = 1, N_AB = 28), # [drug A] at half-max effect
  zeta_B = c(N_S = 1, N_A = 1, N_B = 28, N_AB = 28), # [drug B] at half-max effect
  kappa_A = c(N_S = 1, N_A = 1, N_B = 1, N_AB = 1), # Hill function shape parameter
  kappa_B = c(N_S = 1, N_A = 1, N_B = 1, N_AB = 1), # Hill function shape parameter
  mu = 0.88 * c(N_S = 1, N_A = 0.9, N_B = 0.9, N_AB = 0.81), # maximum growth rate
  delta = 0 * c(N_S = 1, N_A = 1, N_B = 1, N_AB = 1), # intrinsic death rate
  k = 1e14 * c(N_S = 1, N_A = 1, N_B = 1, N_AB = 1), # [R] at half-max growth rate
  alpha = c(N_S = 1, N_A = 1, N_B = 1, N_AB = 1), # resources used per replication
  supply = 0, # resource supply rate
  N_add = 0 # additional resource added at each bottleneck, above base media
  ) {
  # Define the parameters of the model
  config <- as.list(environment())
  config$influx <- influx * c(zeta_A["N_S"], zeta_B["N_S"]) # normalise drug units
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

log_plot <- function(solutions, type = "all", use = c("N_S", "N_A", "N_B", "N_AB")) {
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
    filter(rep == min(rep) & variable %in% c("C_A", "C_B")) %>%
    group_by(variable) %>%
    mutate(value = value / max(value)) %>%
    mutate(across(value, ~replace_na(.x, 0))) %>%
    ungroup() %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(xmin = lag(time), xmax = time) %>%
    filter(!is.na(xmin)) %>%
    select(xmin, xmax, C_A, C_B)
  peak <- max(solutions[solutions$variable %in% use, "value"], na.rm = TRUE)
  # Create the plot
  plot <- ggplot() +
    # Add the gradient backgrounds
    geom_rect(data = background_df,
      aes(xmin = xmin, xmax = xmax, ymin = peak * 10^0.2,
         ymax = peak * 10^0.4, fill = C_A), color = NA, alpha = 1) +
    scale_fill_gradient(low = "white", high = colors[2],
      limits = c(0, 1), name = expression(C[A]), labels = NULL) +
    new_scale_fill() +
    geom_rect(data = background_df,
      aes(xmin = xmin, xmax = xmax, ymin = peak * 10^0.4,
        ymax = peak * 10^0.6, fill = C_B), color = NA, alpha = 1) +
    scale_fill_gradient(low = "white", high = colors[3],
      limits = c(0, 1), name = expression(C[B]), labels = NULL) +
    # Add the lines
    new_scale_fill() +
    geom_line(data = filtered, aes(x = time, y = central, color = variable,
      linetype = factor(rep)), linewidth = 1 + 0.5 / max(filtered$rep)) +
    scale_color_manual(
      values = colors, 
      breaks = c("N_S", "N_A", "N_B", "N_AB", "R"),
      labels = expression(N[S], N[A], N[B], N[AB], R)
    ) +
    # scale_color_manual(values = colors) +
    scale_linetype_manual(values = rep("solid", max(filtered$rep))) +
    guides(linetype = "none") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      breaks = 10^seq(0, 20),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(
      # title = "Bacterial growth over time",
      x = "Time (hours)",
      y = "Population size",
      color = "Strain",
      fill = "Strain"
    ) +
    theme_light() +
    theme(
      legend.position = "bottom",
      # plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 35),
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
  return(plot)
}

system.time(log_plot(simulate()[[1]]))
# simulate(deterministic = TRUE)[[1]]$value[7001:7007]
# [1] 3.051810e+14 2.515378e+11 1.181122e+12 1.717214e+09 6.784406e+12
# [6] 1.333550e-02 9.662624e-01