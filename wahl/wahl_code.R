library(adaptivetau)
library(deSolve)
library(tidyverse)
library(ggnewscale)
library(future)
library(future.apply)

# a growth rate function for nutrient-limited growth
monod <- function(N, r, k) {
  return(r * ifelse(k == 0, 1, 1 / (1 + k / N)))
}

# a function that reduces all populations by a factor of D, in expectation
bottleneck <- function(state, config) {
  with(config, {
    pops <- state[names] # extract just the cell counts
    N <- state["N"] # extract the nutrient concentration
    if (deterministic) {
        pops <- pops * D
    } else {
        pops <- setNames(rbinom(length(pops), pops, D), names)
    }
    N <- N * D + N0 * (1 - D)
    return(c(pops, N))
  })
}

# a function outputting the transitions that can occur in the model
make_transitions <- function(names) {
    # Create the initial list
    rate_list <- list()
    # Add the growth rates for each mutant
    for (i in names) {
      rate_list[[paste0(i, "_growth")]] <- setNames(c(+1), i)
    }
    # Add the nutrient depletion rate
    rate_list[["N_depletion"]] <- setNames(c(-1), "N")
  return(rate_list)
}

# compute the rate at which each transition occurs
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # Calculate replication rates
    replication_rates <- state[names] * monod(N, r_max, k)
    # Calculate growth rates including mutations
    growth_rates <- t(replication_rates) %*% mutation_matrix
    # Calculate nutrient depletion rate
    N_depletion <- sum(replication_rates * alpha)
    # Combine all rates and return
    rate_list <- c(growth_rates, N_depletion)
    return(setNames(rate_list, c(names, "N")))
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
    transitions <- make_transitions(names)
    # Define the fitness of each genotype
    fitness <- rexp(loci, 1 / s)
    config$r_max <- r * (1 + as.matrix(genotypes) %*% fitness)
    # Initialise the state variables
    state <- init
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
        if (is.numeric(seed)) set.seed(seed + t) # set the seed for reproducibility
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
      # Run the bottleneck and update the state
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
  tau = 3, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  N0 = 1e9, # initial nutrient concentration
  m1 = 1e-9, # rate of mutations conferring resistance to drug 1
  init_W = 1e8, # initial wild type population size
  init_M = 0, # initial mutant population sizes
  loci = 1, # number of loci where beneficial mutations may occur
  r = 1, # wild type growth rate with infinite resources
  s = 0.1, # mean fitness effect size of beneficial mutation
  k = 1e8, # [N] at half-max growth rate
  alpha = 1 # nutrients used per replication
  ) {
  if(loci > 10) {
    stop("Too many loci")
  }
  # Define the starting state
  names <- paste0("G",intToBin(1:2^loci - 1))
  init <- setNames(c(init_W, rep(init_M, 2^loci - 1), N0), c(names, "N"))
  # Generate all possible binary strings of length loci
  genotypes <- expand.grid(replicate(loci, c(0, 1), simplify = FALSE))
  genotypes <- setNames(rev(genotypes), paste0("M", 1:loci))
  # Initialize the mutation matrix
  mutation_matrix <- diag(1 - m1, 2 ^ loci)
  # Loop over all genotypes
  for (i in 1:2^loci) {
    for (j in 1:2^loci) {
      # Compute the Hamming distance between the two genotypes
      hamming_distance <- sum(genotypes[i, ] != genotypes[j, ])
      # hamming_distance <- sum(strsplit(names[i], "")[[1]] != strsplit(names[j], "")[[1]])
      # If the Hamming distance is 1, set the mutation rate
      if (hamming_distance == 1) {
        mutation_matrix[i, j] <- m1 / loci
      }
    }
  }
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
    filter(!(variable %in% c("N"))) %>%
    mutate(variable = factor(variable, levels = unique(variable))) %>%
    arrange(variable)
  # Initialise the colours
  colors <- c("black", hcl.colors(length(levels(filtered$variable)) - 1, "Dark 2"))
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

system.time(log_plot(simulate(loci = 4)[[1]]))
# simulate(deterministic = TRUE)[[1]]$value[3001:3003]