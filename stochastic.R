library(adaptivetau)
library(ggplot2)
library(reshape2)

# a pharmacodynamic function for antibioti-induced killing of bacteria
hill <- function(A, params) {
  psi <- params[1]
  phi <- params[2]
  zeta <- params[3]
  kappa <- params[4]
  return(phi * (A / zeta)^kappa / ((A / zeta)^kappa - (psi - phi) / psi))
}

# a growth rate function for nutrient-limited growth
monod <- function(N, params) {
  mu <- params[1]
  k <- params[2]
  return(mu * N / (N + k))
}

# a function specifying the amount of nutrients depleted by the bacteria
deplete <- function(state,alpha,monod_params) {
  return(sum(sapply(seq_along(alpha), function(i) {
    alpha[i] * monod(state["N"],monod_params[i,]) * state[i]
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

# Define a bottleneck function that reduces all state by a fixed fraction
bottleneck <- function(state,
  pattern = c(1,1),
  D = 0.1,
  N0 = 100,
  influx = c(10,10),
  cycl = FALSE,
  pharmacokinetic = FALSE,
  deterministic = FALSE) {
  if (cycl) {
    pattern <- update_pattern(state["prev"])
  }
  # if deterministic then dilute the population by a factor of D
  # else use a binomial distribution to model dilution
  populations <- state[c("S", "R1", "R2", "R12")]
  if (deterministic) {
    diluted <- populations * D
  } else {
    diluted <- rbinom(length(populations), populations, D)
    # set the names of the diluted populations
    names(diluted) <- c("S", "R1", "R2", "R12")
  }
  state <- c(diluted, N = N0,
    pharmacokinetic * state[c("A1", "A2")] + influx * pattern,
    state["prev"] %% 2 + 1) # update the previous drug from 2 to 1 or 1 to 2
  state[state < 0] <- 0 # backup - shouldn't be needed
  return(state)
}

# a function returning the transitions that can occur in the model
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
  c(A1 = -1), # drug 1 depletion
  c(A2 = -1), # drug 2 depletion
  c(N = -1) # nutrient depletion
))
}

# a function returning the rate at which each transition occurs,
# given the current state and parameters
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
    S_death <- S* (hill(A = A1, params["S", c("psi", "phi1", "zeta1", "kappa1")])
      + hill(A = A2, params["S", c("psi", "phi2", "zeta2", "kappa2")]))
    R1_death <- R1 * (hill(A = A1, params["R1", c("psi", "phi1", "zeta1", "kappa1")])
      + hill(A = A2, params["R1", c("psi", "phi2", "zeta2", "kappa2")]))
    R2_death <- R2 * (hill(A = A1, params["R2", c("psi", "phi1", "zeta1", "kappa1")])
      + hill(A = A2, params["R2", c("psi", "phi2", "zeta2", "kappa2")]))
    R12_death <- R12 * (hill(A = A1, params["R12", c("psi", "phi1", "zeta1", "kappa1")])
      + hill(A = A2, params["R12", c("psi", "phi2", "zeta2", "kappa2")]))
    R1_mutation <- S * config$m1 * (1 - config$m2) * monod(N, params["S", c("mu", "k")])
    R2_mutation <- S * config$m2 * (1 - config$m1) * monod(N, params["S", c("mu", "k")])
    R12_mutation <- S * config$m1 * config$m2 * monod(N, params["S", c("mu", "k")]) +
      R1 * config$m2 * monod(N, params["R1", c("mu", "k")]) +
      R2 * config$m1 * monod(N, params["R2", c("mu", "k")])
    HGT_MDR_loss <- config$HGT * R12 * S
    HGT_MDR_gain <- config$HGT * R1 * R2
    A1_depletion <- config$d1 * A1
    A2_depletion <- config$d2 * A2
    N_depletion <- deplete(state, params[, "alpha"],params[, c("mu","k")])
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
      A1_depletion = unname(A1_depletion),
      A2_depletion = unname(A2_depletion),
      N_depletion = unname(N_depletion)
    ))
  })
  
}

# Define a function to plot the resulting solution
bacteria_plot <- function(solution) {
  plot(solution[, 1], solution[, 2], type = "l", col = 1,
    main = "Bacterial growth over time",
    xlab = "Time", ylab = "Population size")
  lines(solution[, 1], solution[, 3], col = 2)
  lines(solution[, 1], solution[, 4], col = 3)
  lines(solution[, 1], solution[, 5], col = 4)
  lines(solution[, 1], solution[, 6], col = 5)
  lines(solution[, 1], solution[, 7], col = 6)
  lines(solution[, 1], solution[, 8], col = 7)
  legend("topright",
    legend = c("S", "R1", "R2", "R12", "Nutrient", "A1", "A2"),
    col = c(1, 2, 3, 4, 5, 6, 7), lty = 1)
}

log_plot <- function(solution) {
  df <- reshape2::melt(
    as.data.frame(solution), id.vars = "time")

  # Filter the df to only include columns S, R1, R2, and R12
  df <- df[df$variable %in% c("S", "R1", "R2", "R12"), ]
  solution_df <- as.data.frame(solution)
  background_df <- data.frame(
    xmin = solution_df$time,
    xmax = c(solution_df$time[-1], solution_df$time[length(solution_df$time)]),
    ymin = 0,
    ymax = max(df$value, na.rm = TRUE),
    A1 = solution_df$A1/max(solution_df$A1),
    A2 = solution_df$A2/max(solution_df$A2)
  )

  # Create the plot
  plot <- ggplot() +
    # Add the gradient background
    geom_rect(
      data = background_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                fill = A1 - A2),
      color = NA
    ) +
    scale_fill_gradient2(low = "#bed5ec", mid = "white", high = "#ecbed5", limits = c(-1, 1),
                         name = NULL, breaks = c(-1,1), labels = c("A1","A2")) +
    # # Add the lines
    geom_line(data = df, aes(x = time, y = value, color = variable), linewidth = 1.5) +
    # scale_y_continuous(trans=scales::pseudo_log_trans(sigma = 10, base = 10)) +
    scale_y_log10() +
    labs(
      title = "Bacterial growth over time",
      x = "Time",
      y = "Population Size",
      color = NULL
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


  # Display the plot
  print(plot)
}

# Define a function to simulate the model
simulate_s <- function(
  pharmacokinetic = FALSE, # should be either TRUE or FALSE
  stewardship = "cycl", #  "cycl" or "comb" or "1_only" or "2_only"
  deterministic = c(
    growth = FALSE,
    death = FALSE,
    mutation = FALSE,
    HGT = FALSE,
    drug_depletion = FALSE,
    nutrient_depletion = FALSE,
    dilution = FALSE
  ), # vector of logicals indicating which transitions should be deterministic
  time = 50, # time to simulate, in hours
  dt = 0.01, # time step, in hours
  freq = 10, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  N0 = 100, # initial nutrient concentration
  HGT = 0.000, # rate of horizontal gene transfer
  m1 = 0.1, # rate of mutations conferring resistance to drug 1
  m2 = 0.1, # rate of mutations conferring resistance to drug 2
  d1 = 1 * pharmacokinetic, # rate of drug 1 elimination
  d2 = 1 * pharmacokinetic, # rate of drug 2 elimination
  influx = c(A1 = 10, A2 = 10), # drug influx concentrations
  # lists of genotype-specific parameters, in the order S, R1, R2, R12
  init = c(S = 100, R1 = 0, R2 = 0, R12 = 0), # initial population sizes
  psi = c(0.1, 0.1, 0.1, 0.1), # growth rate with no drugs
  phi1 = c(0.2, 0.2, 0.2, 0.2), # maximum killing rate for drug 1
  phi2 = c(0.2, 0.2, 0.2, 0.2), # maximum killing rate for drug 2
  zeta1 = c(1, 100, 1, 100), # MIC drug 1
  zeta2 = c(1, 1, 100, 100), # MIC drug 2
  kappa1 = c(1, 1, 1, 1), # Hill function steepness parameter drug 1
  kappa2 = c(1, 1, 1, 1), # Hill function steepness parameter drug 2
  mu = c(1, 1, 1, 1), # growth rate with infinite resources
  k = c(100, 100, 100, 100), # resource concentration at half-maximal growth
  alpha = c(1, 1, 1, 1) # resources used per unit growth
  ) {
  # Define the parameters of the model
  config <- list(
    D = D,
    N0 = N0,
    HGT = HGT,
    m1 = m1,
    m2 = m2,
    d1 = d1,
    d2 = d2,
    influx = influx,
    stewardship = stewardship,
    pharmacokinetic = pharmacokinetic,
    deterministic = c(
      S_growth = deterministic["growth"],
      R1_growth = deterministic["growth"],
      R2_growth = deterministic["growth"],
      R12_growth = deterministic["growth"],
      S_death = deterministic["death"],
      R1_death = deterministic["death"],
      R2_death = deterministic["death"],
      R12_death = deterministic["death"],
      R1_mutation = deterministic["mutation"],
      R2_mutation = deterministic["mutation"],
      R12_mutation = deterministic["mutation"],
      HGT_MDR_loss = deterministic["HGT"],
      HGT_MDR_gain = deterministic["HGT"],
      A1_depletion = deterministic["drug_depletion"],
      A2_depletion = deterministic["drug_depletion"],
      N_depletion = deterministic["nutrient_depletion"]
    )
  )
  if (stewardship == "2_only") {
    config$pattern <- c(0, 1)
  } else if (stewardship == "comb") {
    config$pattern <- c(1, 1)
  } else {
    config$pattern <- c(1, 0)
  }
  config$params <- cbind(psi, phi1, phi2, zeta1, zeta2, kappa1, kappa2, mu, k, alpha)
  rownames(config$params) <- c("S", "R1", "R2", "R12")

  # Define the transitions of the model
  transitions <- make_transitions()

  # Run the simulation
  state <- c(init, N = config$N0, config$influx * config$pattern,
    prev = sum(config$pattern * c(1, 2)))
  t <- 0
  while (t < time) {
    # Run the model between bottlenecks
    new_solution <- ssa.adaptivetau(state, transitions, rates, config, tf = freq,
      deterministic = config$deterministic)
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
      deterministic = deterministic["dilution"]
    )

    # Update the solution
    if (t==0) {solution <- new_solution}
    else {solution <- rbind(solution, new_solution)}
    # Update the time
    t <- t + freq
  }

  return(solution)
}

# bacteria_plot(simulate_s(pharmacokinetic = TRUE))
log_plot(simulate_s(N0=1e5))
