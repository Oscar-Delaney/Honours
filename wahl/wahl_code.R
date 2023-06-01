library(adaptivetau)
library(deSolve)
library(tidyverse)
library(ggnewscale)
library(future)
library(future.apply)

# a growth rate function for nutrient-limited growth
monod <- function(N, mu, k) {
  return(mu * ifelse(k == 0, 1, 1 / (1 + k / N)))
}

# a function that reduces all populations by a factor of D, in expectation
bottleneck <- function(state, config) {
  with(config, {
    pops <- state[names(init)] # extract just the cell counts
    N <- state["N"] # extract the nutrient concentration
    if (deterministic) {
        pops <- pops * D
    } else {
        pops <- setNames(rbinom(length(pops), pops, D), names(pops))
    }
    N <- N * D + N0 * (1 - D)
    return(c(pops, N))
  })
}

# a function outputting the transitions that can occur in the model
make_transitions <- function() {
  return(list(
  W_growth = c(W = +1), # wild type growth
  M_growth = c(M = +1), # mutant growth
  N_depletion = c(N = -1) # nutrient depletion
))
}

# compute the rate at which each transition occurs
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # Calculate replication rates
    replication_rates <- c(W, M) * monod(N, mu, k)
    # chance of a replication in row i resulting in a strain j cell
    mutation <- matrix(c(
      W  = c((1 - m1), m1),
      M = c(0, 1)
      ), nrow = 2, byrow = TRUE)
    # Calculate growth rates including mutations
    growth_rates <- replication_rates %*% mutation
    # Calculate nutrient depletion rate
    N_depletion <- sum(replication_rates * alpha)
    # Combine all rates and return
    rate_list <- c(growth_rates, N_depletion)
    return(setNames(rate_list, names(make_transitions())))
  })
}

# a function to convert the transition rates into an ODE function
ode_rates <- function(t, state, config) {
  with(as.list(rates(state, config, t)), {
    return(list(c(
      W = W_growth,
      M = M_growth,
      N = -N_depletion
      )))
  })
}

# a function to implement one run of the model
single_run <- function(config, x) {
  with(config, {
    # Define the transitions of the model
    transitions <- make_transitions()
    # Initialise the state variables
    state <- c(init, N = N0)
    time_grid <- seq(0, time, by = dt) # a common time grid for all runs
    bottlenecks <- unique(round(c(seq(0, time, tau), time), 10))
    for (t in bottlenecks[-length(bottlenecks)]) {
      # Determine the time until the next bottleneck or dose
      end <- min(bottlenecks[bottlenecks > t] - t)
      # Run the model between bottlenecks, deterministically or stochastically
      if (deterministic) {
        times <- c(time_grid[time_grid <= end], end) # ensures length(times) > 1
        new <- ode(state, times, ode_rates, config)
      } else {
        if (is.numeric(seed)) set.seed(seed) # set the seed for reproducibility
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
      # Run the bottleneck and/or dose and update the state
      state <- bottleneck(new[nrow(new), ], config)
    }
    # Interpolate the solution to the common time grid
    approx_vars <- lapply(colnames(solution), function(var) {
      approx(solution[, "time"], solution[, var], xout = time_grid)$y
    })
    solution_interpolated <- data.frame(
      setNames(approx_vars, colnames(solution)),
      rep = x
    )
    return(solution_interpolated)
  })
}

# a function to simulate the model
simulate <- function(
  rep = 1, # number of runs of the simulation
  seed = NULL, # seed for reproducibility
  deterministic = FALSE, # should be either TRUE or FALSE
  time = 100, # time to simulate, in hours
  dt = 0.1, # time step, in hours
  max_step = Inf, # SSA max step parameter, only used if deterministic = FALSE
  tau = 10, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  N0 = 1e15, # initial nutrient concentration
  m1 = 1e-9, # rate of mutations conferring resistance to drug 1
  init = c(W = 1e12, M = 0), # initial population sizes
  mu = 0.88 * c(W = 1, M = 1.1), # maximum growth rate
  k = 1e14 * c(W = 1, M = 1), # [N] at half-max growth rate
  alpha = c(W = 1, M = 1) # nutrients used per replication
  ) {
  # Define the parameters of the model
  config <- as.list(environment())
  # Run the simulation rep number of times, using parallelisation if possible
  plan(multisession) # compatible with both unix and Windows
  solutions <- bind_rows(future_lapply(1:rep, function(x) {
    single_run(config, x)},
    future.seed = TRUE))
  # Convert the solutions to long format
  long <- pivot_longer(solutions, cols = -c(time, rep), names_to = "variable")
  return(list(long, config))
}

log_plot <- function(solutions, type = "all") {
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
    filter(variable %in% c("W", "M")) %>%
    mutate(variable = factor(variable, levels = c("W", "M"))) %>%
    arrange(variable)
  # Initialise the colours
  colors <- c("black", "navy")
  peak <- max(solutions$value, na.rm = TRUE)
  # Create the plot
  plot <- ggplot() +
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

system.time(log_plot(simulate()[[1]], type = "all"))
# simulate(deterministic = TRUE)[[1]]$value[7001:7007]