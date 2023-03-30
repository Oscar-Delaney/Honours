library(deSolve)

simulate_d <- function(
  pharmacokinetic = FALSE, # should be either TRUE or FALSE
  stewardship = "cycl", #  "cycl" or "comb" or "1_only" or "2_only"
  time = 50, # time to simulate, in hours
  dt = 0.01, # time step, in hours
  freq = 10, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  N0 = 100, # initial nutrient concentration
  HGT = 0.0000, # rate of horizontal gene transfer
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
    pharmacokinetic = pharmacokinetic
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
    return(sum(sapply(1:length(alpha), function(i) {
      alpha[i] * monod(state["N"],monod_params[i,]) * state[i]
    })))
  }
  # create the right pattern based on the previous drug administered
  update_pattern <- function(prev = 1) {
    # if the previous drug administered was 1, then the next drug is 2
    if (prev == 1) {
      return(c(0, 1))
    }
    return(c(1, 0))
  }

  # Define a bottleneck function that reduces all state by a fixed fraction
  bottleneck <- function(t, state, config = NULL) {
    pattern <- config$pattern
    if (config$stewardship == "cycl") {
      pattern <- update_pattern(state["prev"])
    }
    state <- c(state[c("S", "R1", "R2", "R12")] * config$D, N = config$N0,
      config$pharmacokinetic * state[c("A1", "A2")] + config$influx * pattern,
      state["prev"] %% 2 + 1) # update the previous drug from 2 to 1 or 1 to 2
    state[state < 0] <- 0 # backup - shouldn't be needed
    return(state)
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


  # Define the differential equations for the model
  bacterial_growth <- function(t, state, config = NULL) {
    params <- config$params
    with(as.list(c(state)), {
      S <- state[1]
      R1 <- state[2]
      R2 <- state[3]
      R12 <- state[4]
      N <- state[5]
      A1 <- state[6]
      A2 <- state[7]
      net_recombination <- config$HGT * (R1 * R2 - R12 * S)
      dS <- S * ((1 - config$m1) * (1 - config$m2) *
        monod(N, params["S", c("mu", "k")]) -
        hill(A = A1, params["S", c("psi", "phi1", "zeta1", "kappa1")]) -
        hill(A = A2, params["S", c("psi", "phi2", "zeta2", "kappa2")])) +
        net_recombination
      dR1 <- R1 * ((1 - config$m2) * monod(N, params["R1", c("mu", "k")]) -
        hill(A = A1, params["R1", c("psi", "phi1", "zeta1", "kappa1")]) -
        hill(A = A2, params["R1", c("psi", "phi2", "zeta2", "kappa2")])) +
        S * config$m1 * (1 - config$m2) * monod(N, params["S", c("mu", "k")]) -
        net_recombination
      dR2 <- R2 * ((1 - config$m1) * monod(N, params["R2", c("mu", "k")]) -
        hill(A = A1, params["R2", c("psi", "phi1", "zeta1", "kappa1")]) -
        hill(A = A2, params["R2", c("psi", "phi2", "zeta2", "kappa2")])) +
        S * (1 - config$m1) * config$m2 * monod(N, params["S", c("mu", "k")]) -
        net_recombination
      dR12 <- R12 * (monod(N, params["R12", c("mu", "k")]) -
        hill(A = A1, params["R12", c("psi", "phi1", "zeta1", "kappa1")]) -
        hill(A = A2, params["R12", c("psi", "phi2", "zeta2", "kappa2")])) +
        S * config$m1 * config$m2 * monod(N, params["S", c("mu", "k")]) +
        R1 * config$m2 * monod(N, params["R1", c("mu", "k")]) +
        R2 * config$m1 * monod(N, params["R2", c("mu", "k")]) +
        net_recombination
      dN <- -deplete(state, params[, "alpha"],params[, c("mu","k")])
      dA1 <- -config$d1 * A1
      dA2 <- -config$d2 * A2
      return(list(c(S = dS, R1 = dR1, R2 = dR2, R12 = dR12, N = dN, A1 = dA1, A2 = dA2, prev = 0)))
    })
  }

  # Run the simulation
  times <- seq(0, time, dt)
  state <- c(init, N = N0, influx * config$pattern,
    prev = sum(config$pattern * c(1, 2)))
  events <- list(func = bottleneck, times = seq(freq, time, freq))
  solution <- ode(state, times, bacterial_growth, config, events = events)
  bacteria_plot(solution)
}

simulate_d()
