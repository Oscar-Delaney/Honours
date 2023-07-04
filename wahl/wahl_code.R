library(adaptivetau)
library(deSolve)
library(tidyverse)
library(ggnewscale)
library(future)
library(future.apply)
library(R.utils)

# a growth rate function for nutrient-limited growth
monod <- function(R, r, k) {
  return(r * ifelse(k == 0, 1, 1 / (1 + k / R)))
}

# fitness advantages of new mutants
mutant_fitness <- function(num_mutants, r, w, names) {
  return(setNames(r * (1 + c(0, rexp(num_mutants, 1 / w))), names))
}

# a function that reduces all populations by a factor of D, in expectation
bottleneck <- function(state, config) {
  with(config, {
    pops <- setNames(rbinom(length(names), state[names], D), names)
    R <- state["R"] * D + R0 * (1 - D)
    return(c(pops, R))
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
    rate_list[["N_depletion"]] <- setNames(c(-1), "R")
  return(rate_list)
}

# compute the rate at which each transition occurs
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # Calculate replication rates
    replication_rates <- state[names] * monod(R, r_vec, k)
    if (is.numeric(num_mutants)) {
        # find which mutant category we should be filling
        to_mutate <- min(which(state == 0), num_mutants + 2) - 1
        # chance of a replication in row i resulting in a strain j cell
        mutation <- diag(num_mutants + 1)
        if (to_mutate <= num_mutants) {
          mutation[1, to_mutate + 1] <- mu
          mutation[1, 1] <- 1 - mu
        } else {
          print("Out of mutant spots! :(")
        }
    }
    # Calculate growth rates including mutations
    growth_rates <- replication_rates %*% mutation
    # Calculate nutrient depletion rate
    N_depletion <- sum(replication_rates * alpha)
    # Combine all rates and return
    rate_list <- c(growth_rates, N_depletion)
    return(setNames(rate_list, c(names, "R")))
  })
}

# a function to implement one run of the model
single_run <- function(config, x) {
  with(config, {
    # Define the transitions of the model
    transitions <- make_transitions(names)
    # Initialise the state variables
    state <- init
    if (is.numeric(num_mutants)) {
      config$r_vec <- mutant_fitness(num_mutants, r, w, names)
    } else {
      config$r_vec <- r * t(1 + as.matrix(genotypes) %*% rexp(loci, 1 / w))
    }
    time_grid <- seq(0, time, by = dt) # a common time grid for all runs
    bottlenecks <- unique(round(c(seq(0, time, tau), time), 10))
    for (t in bottlenecks[-length(bottlenecks)]) {
      # Determine the time until the next bottleneck or dose
      end <- min(bottlenecks[bottlenecks > t] - t)
      # Create a new vector of growth rates for mutants, if locus-agnostic
      if (is.numeric(num_mutants)) {
        config$r_vec <- ifelse(state[0:num_mutants + 1] == 0,
          mutant_fitness(num_mutants, r, w, names), config$r_vec)
      }
      # set the seed for reproducibility
      if (is.numeric(seed)) set.seed(round(seed + (x * time + t) / tau))
      # Run the model between bottlenecks
      new <- ssa.adaptivetau(
        state, transitions, rates, config, tf = end,
        tl.params = list(maxtau = max_step),
        deterministic = grep("depletion", names(transitions))
      )
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
    # add a new row at the end with the r_vec
    solution_interpolated[nrow(solution_interpolated) + 1, ] <-
      c(pi * 1e6, (config$r_vec / r) - 1, 0, x)
    return(solution_interpolated)
  })
}

# a function to simulate the model
simulate <- function(
  rep = 1, # number of runs of the simulation
  seed = NULL, # seed for reproducibility
  time = 100, # time to simulate, in hours
  dt = 0.1, # time step, in hours
  max_step = Inf, # SSA max step parameter
  tau = 3, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  R0 = 1e9, # initial nutrient concentration
  mu = 1e-9, # rate of mutations conferring resistance to drug 1
  N = 1e9, # (1/D) of the initial wild type population
  num_mutants = NULL, # number of mutants
  loci = NULL, # number of mutable loci in the genome
  r = 1, # wild type growth rate with infinite resources
  w = 0.1, # mean fitness effect size of beneficial mutation
  k = 1e8, # [R] at half-max growth rate
  alpha = 1 # nutrients used per replication
  ) {
  # Define the parameters of the model
  if (is.numeric(num_mutants)) {
    names <- c("W", paste0("M", 1:num_mutants))
  } else {
      if (is.null(loci)) stop("Must specify num_mutants XOR loci")
      names <- paste0("G", intToBin(1:2^loci - 1))
      # Generate all possible binary strings of length loci
      genotypes <- expand.grid(replicate(loci, c(0, 1), simplify = FALSE))
      genotypes <- setNames(rev(genotypes), paste0("M", 1:loci))
      # Initialize the mutation matrix
      mutation <- diag(1 - mu, 2 ^ loci)
      # Loop over all genotypes
      for (i in 1:2^loci) {
        for (j in 1:2^loci) {
          # If the Hamming distance is 1, set the mutation rate to be nonzero
          if (sum(genotypes[i, ] != genotypes[j, ]) == 1) {
            mutation[i, j] <- mu / loci
          }
        }
      }
  }
  init <- setNames(c(round(N * D), rep(0, length(names) - 1), R0), c(names, "R"))
  config <- as.list(environment())
  # Run the simulation rep number of times, using parallelisation if possible
  plan(multisession) # compatible with both unix and Windows
  set.seed(seed) # set the seed for reproducibility
  solutions <- bind_rows(future_lapply(1:rep, function(x) {
    single_run(config, x)},
    future.seed = TRUE))
  # Convert the solutions to long format
  long <- pivot_longer(solutions, cols = -c(time, rep), names_to = "variable")
  config$s_all <- long %>% filter(near(time, pi * 1e6)) %>% select(-time)
  long <- long[!near(long$time, pi * 1e6), ]
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
    filter(!(variable %in% c("R"))) %>%
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
