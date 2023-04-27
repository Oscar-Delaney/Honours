library(adaptivetau)
library(deSolve)
library(ggplot2)
library(ggnewscale)
library(reshape2)
library(dplyr)
library(purrr)
library(matrixStats)
library(profvis)
library(parallel)
library(future)
library(future.apply)
library(promises)
N_adj <- 100 # hacky fix to non-negative N requirement in adaptivatau

# a pharmacodynamic function for singele-antibiotic induced killing of bacteria
hill <- function(A, params) {
  psi <- params[1]
  phi <- params[2]
  zeta <- params[3]
  kappa <- params[4]
  return(phi * (A / zeta)^kappa / ((A / zeta)^kappa - (psi - phi) / psi))
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
  mu <- params[1]
  k <- params[2]
  N <- max(0,N - N_adj)
  return(mu * N / (N + k))
}

# a function specifying the amount of nutrients depleted by the bacteria
deplete <- function(state, alpha, monod_params) {
  return(sum(sapply(seq_along(alpha), function(i) {
    alpha[i] * monod(state["N"], monod_params[i, ]) * state[i]
  })))
}

# update the pattern based on the previous drug administered, under cycling
update_pattern <- function(prev = 1) {
  # if the previous drug administered was 1, then the next drug is 2
  if (prev == 1) {
    return(c(0, 1))
  }
  return(c(1, 0))
}

# a function that reduces all populations by a factor of D, in expectation
bottleneck <- function(state,
  pattern = c(1, 1),
  D = 0.1,
  N0 = 100,
  influx = c(10, 10),
  cycl = FALSE,
  pharmacokinetic = FALSE,
  deterministic = FALSE) {
  if (cycl) {
    pattern <- update_pattern(state["prev"])
  }
  populations <- state[c("S", "R1", "R2", "R12")]
  if (deterministic) {
    diluted <- populations * D
  } else {
    diluted <- rbinom(length(populations), populations, D)
  }
  state <- c(diluted, N = N0,
    pharmacokinetic * state[c("A1", "A2")] + influx * pattern,
    state["prev"] %% 2 + 1) # update the previous drug from 2 to 1 or 1 to 2
  state[state < 0] <- 0 # backup - shouldn't be needed
  # set the names of the state variables
  names(state) <- c("S", "R1", "R2", "R12", "N", "A1", "A2", "prev")
  return(state)
}

# a function outputting the transitions that can occur in the model
make_transitions <- function() {
  return(list(
  c(S = +1), # growth in S
  c(R1 = +1), # growth in R1
  c(R2 = +1), # growth in R2
  c(R12 = +1), # growth in R12
  c(S = -1), # drug-induced death in S
  c(R1 = -1), # drug-induced death in R1
  c(R2 = -1), # drug-induced death in R2
  c(R12 = -1), # drug-induced death in R12
  c(R1 = +1), # mutation to R1
  c(R2 = +1), # mutation to R2
  c(R12 = +1), # mutation to R12
  c(S = -1, R1 = +1, R2 = +1, R12 = -1), # HGT loss of MDR
  c(S = +1, R1 = -1, R2 = -1, R12 = +1), # HGT gain of MDR
  c(N = -1), # nutrient depletion
  c(A1 = -1), # drug 1 depletion
  c(A2 = -1) # drug 2 depletion
))
}

# computing transition rates, given the current state and parameters
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # extract the parameters
    params <- config$params
    # extract the current state
    S <- state["S"]
    R1 <- state["R1"]
    R2 <- state["R2"]
    R12 <- state["R12"]
    N <- state["N"]
    A1 <- state["A1"]
    A2 <- state["A2"]
    # compute the rates
    S_growth <- S * (1 - config$m1) * (1 - config$m2) * monod(N, params["S", c("mu", "k")])
    R1_growth <- R1 * (1 - config$m2) * monod(N, params["R1", c("mu", "k")])
    R2_growth <- R2 * (1 - config$m1) * monod(N, params["R2", c("mu", "k")])
    R12_growth <- R12 * monod(N, params["R12", c("mu", "k")])
    S_death <- S * death(A1 = A1, A2 = A2, params["S", ])
    R1_death <- R1 * death(A1 = A1, A2 = A2, params["R1", ])
    R2_death <- R2 * death(A1 = A1, A2 = A2, params["R2", ])
    R12_death <- R12 * death(A1 = A1, A2 = A2, params["R12", ])
    R1_mutation <- S * config$m1 * (1 - config$m2) * monod(N, params["S", c("mu", "k")])
    R2_mutation <- S * config$m2 * (1 - config$m1) * monod(N, params["S", c("mu", "k")])
    R12_mutation <- S * config$m1 * config$m2 * monod(N, params["S", c("mu", "k")]) +
      R1 * config$m2 * monod(N, params["R1", c("mu", "k")]) +
      R2 * config$m1 * monod(N, params["R2", c("mu", "k")])
    HGT_MDR_loss <- config$HGT * R12 * S
    HGT_MDR_gain <- config$HGT * R1 * R2
    N_depletion <- deplete(state, params[, "alpha"], params[, c("mu", "k")])
    A1_depletion <- config$d1 * A1
    A2_depletion <- config$d2 * A2
    # return the rates
    return(c(
      S_growth = unname(S_growth),
      R1_growth = unname(R1_growth),
      R2_growth = unname(R2_growth),
      R12_growth = unname(R12_growth),
      S_death = unname(S_death),
      R1_death = unname(R1_death),
      R2_death = unname(R2_death),
      R12_death = unname(R12_death),
      R1_mutation = unname(R1_mutation),
      R2_mutation = unname(R2_mutation),
      R12_mutation = unname(R12_mutation),
      HGT_MDR_loss = unname(HGT_MDR_loss),
      HGT_MDR_gain = unname(HGT_MDR_gain),
      N_depletion = unname(N_depletion),
      A1_depletion = unname(A1_depletion),
      A2_depletion = unname(A2_depletion)
    ))
  })

}

# a function to convert the transition rates into an ODE function
ode_rates <- function(t, state, config) {
  adaptivetau_rates = rates(state, config, t)
  ode_rates_xyz <- with(as.list(adaptivetau_rates),{
    c(
      S = S_growth - S_death + HGT_MDR_gain - HGT_MDR_loss,
      R1 = R1_growth - R1_death + R1_mutation - HGT_MDR_gain + HGT_MDR_loss,
      R2 = R2_growth - R2_death + R2_mutation - HGT_MDR_gain + HGT_MDR_loss,
      R12 = R12_growth - R12_death + R12_mutation + HGT_MDR_gain - HGT_MDR_gain,
      N = -N_depletion,
      A1 = -A1_depletion,
      A2 = -A2_depletion,
      prev = 0
    )
  })
  return(list(ode_rates_xyz))
}

# a function to implement one run of the model
single_run <- function(config, x) {
  # Define the transitions of the model
  transitions <- make_transitions()
  # Initialise the state variables
  state <- c(config$init, N = config$N0, config$influx * config$pattern,
      prev = sum(config$pattern * c(1, 2)))
  t <- 0
  while (t < max(config$time_grid)) {
    # Run the model between bottlenecks, deterministically or stochastically
    if (config$deterministic) {
      times <- config$time_grid[config$time_grid <= config$freq]
      new_solution <- ode(state, times, ode_rates, config)
    } else {
      new_solution <- ssa.adaptivetau(state, transitions, rates, config,
        tf = config$freq)
    }
    # Make the time column reflect the overall time accurately
    new_solution[, 1] <- new_solution[, 1] + t
    # Run the bottleneck and update the state
    state <- bottleneck(
      state = new_solution[nrow(new_solution), ],
      pattern = config$pattern,
      D = config$D,
      N0 = config$N0,
      influx = config$influx,
      cycl = config$stewardship == "cycl",
      pharmacokinetic = config$pharmacokinetic,
      deterministic = config$deterministic
    )

    # Update the solution
    if (t == 0) {
      solution <- new_solution
    } else {
      solution <- rbind(solution, new_solution[-1, ])
    }
    # Update the time
    t <- t + config$freq
  }
  # Interpolate the solution to the common time grid
  solution <- as.data.frame(solution)
  solution <- data.frame(
    time = config$time_grid,
    S = approx(solution$time, solution$S, xout = config$time_grid)$y,
    R1 = approx(solution$time, solution$R1, xout = config$time_grid)$y,
    R2 = approx(solution$time, solution$R2, xout = config$time_grid)$y,
    R12 = approx(solution$time, solution$R12, xout = config$time_grid)$y,
    N = approx(solution$time, solution$N, xout = config$time_grid)$y,
    A1 = approx(solution$time, solution$A1, xout = config$time_grid,
      method = "constant", f = 1)$y,
    A2 = approx(solution$time, solution$A2, xout = config$time_grid,
      method = "constant", f = 1)$y,
    rep = x
  )
  return(solution)
}

# a function to simulate the model
simulate <- function(
  rep = 1,
  pharmacokinetic = FALSE, # should be either TRUE or FALSE
  deterministic = FALSE, # should be either TRUE or FALSE
  stewardship = "cycl", #  "cycl" or "comb" or "1_only" or "2_only"
  time = 100, # time to simulate, in hours
  dt = 0.01, # time step, in hours
  freq = 10, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  N0 = 1e15, # initial nutrient concentration
  HGT = 0, # rate of horizontal gene transfer
  m1 = 1e-9, # rate of mutations conferring resistance to drug 1
  m2 = 1e-9, # rate of mutations conferring resistance to drug 2
  d1 = log(2) / 3.5 * pharmacokinetic, # rate of drug 1 elimination
  d2 = log(2) / 3.5 * pharmacokinetic, # rate of drug 2 elimination
  influx = 119 * c(A1 = 1, A2 = 1), # drug influx concentrations
  # lists of genotype-specific parameters, in the order S, R1, R2, R12
  init = c(S = 1e10, R1 = 0, R2 = 0, R12 = 0), # initial population sizes
  psi = 0.3 * c(1, 1, 1, 1), # growth rate with no drugs
  phi1 = 2 * psi, # maximum reduction in fitness for drug 1
  phi2 = 2 * psi, # maximum reduction in fitness for drug 2
  zeta1 = 17 * c(1, 28, 1, 28), # MIC drug 1
  zeta2 = 17 * c(1, 1, 28, 28), # MIC drug 2
  kappa1 = c(1, 1, 1, 1), # Hill function steepness parameter drug 1
  kappa2 = c(1, 1, 1, 1), # Hill function steepness parameter drug 2
  theta = c(0, 0, 0, 0), # drug interaction term
  mu = 0.88 * c(1, 0.9, 0.9, 0.81), # growth rate with infinite resources
  k = c(100, 100, 100, 100), # resource concentration at half-maximal growth
  alpha = c(1, 1, 1, 1) # resources used per unit growth
  ) {
  # Define the parameters of the model
  config <- list(
    D = D,
    N0 = N0 + N_adj, # hacky adjustment to avoid adaptivetau nonnegative issue
    init = init,
    HGT = HGT,
    m1 = m1,
    m2 = m2,
    d1 = d1,
    d2 = d2,
    influx = influx,
    stewardship = stewardship,
    pharmacokinetic = pharmacokinetic,
    deterministic = deterministic,
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
  config$params <- cbind(psi, phi1, phi2, zeta1, zeta2, kappa1, kappa2, theta, mu, k, alpha)
  rownames(config$params) <- c("S", "R1", "R2", "R12")
  # Run the simulation rep number of times, using parallelisation if possible
  plan(multisession) # compatible with both unix and Windows
  solutions <- future_lapply(1:rep, function(x) {
    single_run(config, x)},
    future.seed = TRUE)
  return(melt(do.call(rbind, solutions), id.vars = c("time", "rep")))
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
  filtered <- filter(summary, variable %in% c("S", "R1", "R2", "R12"))
  # Initialise the colours
  colors <- c("black", "navy", "#800000", "#008000")
  # Create antibiotic concentrations data frame
  times <- unique(solutions$time)
  background_df <- data.frame(
    xmin = times[-length(times)],
    xmax = times[-1],
    A1 = solutions[solutions$rep == 1 & solutions$variable == "A1", "value"][-1] /
      max(solutions[solutions$variable == "A1", "value"]),
    A2 = solutions[solutions$rep == 1 & solutions$variable == "A2", "value"][-1] /
      max(solutions[solutions$variable == "A2", "value"])
  )
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

system.time(log_plot(simulate(deterministic = T, stewardship = "cycl",time=100,rep = 10,m1=0,m2=0,N0=1e8, init = c(S=1e6,R1=0,R2=0,R12=0)), type = "all"))
log_plot(simulate(), type = "all")
