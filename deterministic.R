library(deSolve)

simulate <- function(pharmacokinetic = FALSE, stewardship = "cycl"){
  # pharmacokinetic should be either TRUE or FALSE
  # stewardship should be "cycl" or "comb" or "1_only" or "2_only"

  # Define the parameters of the model
  config <- list(
    D = 0.1, # bottleneck dilution ratio
    N0 = 100, # initial nutrient concentration
    m1 = 0.1, # rate of mutations conferring resistance to drug 1
    m2 = 0.1, # rate of mutations conferring resistance to drug 2
    d1 = 1 * pharmacokinetic, # rate of drug 1 elimination
    d2 = 1 * pharmacokinetic, # rate of drug 2 elimination
    i1 = 10, # drug 1 influx concentration
    i2 = 10, # drug 2 influx concentration
    stewardship = stewardship, # stewardship strategy
    pharmacokinetic = pharmacokinetic # pharmacokinetic model
  )
  if (stewardship == "2_only") {
    config$pattern <- c(0, 1)
  } else if (stewardship == "comb") {
    config$pattern <- c(1, 1)
  } else {
    config$pattern <- c(1, 0)
  }

  # Initialize 'params' matrix with 4 rows and 10 columns
  params <- matrix(nrow = 4, ncol = 12)

  # Name the rows of the matrix
  rownames(params) <- c("S", "R1", "R2", "R12")

  # Name the columns of the matrix
  colnames(params) <- c("psi", "phi1", "zeta1", "kappa1", "phi2", "zeta2", "kappa2", "mu", "k", "alpha", "m_rate", "r_rate")

  # Assign values to the 'psi' column of the matrix
  params[, "psi"] <- c(0.1, 0.1, 0.1, 0.1)

  # Assign values to the 'phi1' and 'phi2' columns of the matrix
  params[, c("phi1", "phi2")] <- c(0.2, 0.2, 0.2, 0.2)

  # Assign values to the 'zeta1' column of the matrix
  params[, "zeta1"] <- c(1, 100, 1, 100)

  # Assign values to the 'zeta2' column of the matrix
  params[, "zeta2"] <- c(1, 1, 100, 100)

  # Assign values to the 'kappa1' and 'kappa2' columns of the matrix
  params[, c("kappa1", "kappa2")] <- c(1, 1, 1, 1)

  # Assign values to the 'mu' column of the matrix
  params[, "mu"] <- c(1, 1, 1, 1)

  # Assign values to the 'k' column of the matrix
  params[, "k"] <- c(100, 100, 100, 100)

  # Assign values to the 'alpha' column of the matrix
  params[, "alpha"] <- c(1, 1, 1, 1)

  # Assign values to the 'm_rate' column of the matrix
  params[, "m_rate"] <- c(0.1, 0.1, 0.1, 0.1)

  # Assign values to the 'r_rate' column of the matrix
  params[, "r_rate"] <- c(0.1, 0.1, 0.1, 0.1)

  config$params = params

  # a pharmacodynamic function for antibioti-induced killing of bacteria
  hill <- function(A, parms) {
    psi = parms[1] ; phi = parms[2] ; zeta = parms[3]; kappa = parms[4]
    return(phi*(A/zeta)^kappa/((A/zeta)^kappa-(psi-phi)/psi))
  }

  # a growth rate function for nutrient-limited growth
  monod <- function(N,parms){
    mu = parms[1] ; k = parms[2]
    return(mu*N/(N+k))
  }

  # a growth rate function using hill and monod functions
  growth <- function(A1,A2,N,parms) {
    return(monod(N,parms[8:9])-hill(A1,parms[1:4])-hill(A2,parms[c(1,5:7)]))
  }

  # a function specifying the amount of nutrients depleted by the bacteria
  deplete <- function(dS=0,dR1=0,dR2=0,dR12=0,alpha=params[,"alpha"]) {
    return(sum(-alpha*max(0,c(dS,dR1,dR2,dR12))))
  }

  # Define a function to update the pattern for the cyclic treatment strategy
  update_pattern <- function(current_pattern, stewardship) {
    if (stewardship == "cycl") {return(1 - current_pattern)}
    return(current_pattern)
  }

  # Define a bottleneck function that reduces all populations by a fixed fraction
  bottleneck <- function(t, populations, parms = params, current_pattern = pattern) {
    populations <- c(populations[1:4] * config$D, config$N0, config$pharmacokinetic * populations[6:7] + c(config$i1, config$i2) * current_pattern)
    populations[populations < 0] <- 0 # backup - shouldn't be needed
    if (stewardship == "cycl") {updated_pattern <- 1 - current_pattern}
    updated_pattern <- current_pattern
    return(list(populations = populations, pattern = updated_pattern))
  }
  # bottleneck <- function(t, populations, parms = params) {
  #   populations <- c(populations[1:4]*D, N0, pharmacokinetic*populations[6:7]+c(i1,i2)*config$pattern)
  #   # names(populations) <- c("S", "R1", "R2", "R12", "N", "A1", "A2")
  #   populations[populations < 0] <- 0 # backup - shouldn't be needed
  #   if (stewardship=="cycl") {config$pattern <<- 1 - config$pattern}
  #   return(populations)
  #   }

  # Define a function to plot the resulting solution
  bacteria_plot <- function(solution) {
    plot(solution[, 1], solution[,2],type = "l", col = 1,
    main = "Bacterial growth over time",
    xlab = "Time", ylab = "Population size")
    lines(solution[, 1],solution[, 3], col = 2)
    lines(solution[, 1],solution[, 4], col = 3)
    lines(solution[, 1],solution[, 5], col = 4)
    lines(solution[, 1],solution[, 6], col = 5)
    lines(solution[, 1],solution[, 7], col = 6)
    lines(solution[, 1],solution[, 8], col = 7)
    legend("topright", legend = c("S", "R1", "R2", "R12", "Nutrient","A1","A2"),
      col = c(1, 2, 3, 4, 5, 6, 7), lty = 1)
  }


  # Define the differential equations for the model
  bacterial_growth <- function(t, populations, config = NULL) {
    parms <- config$params
    with(as.list(c(populations)), {
      S <- populations[1]
      R1 <- populations[2]
      R2 <- populations[3]
      R12 <- populations[4]
      N <- populations[5]
      A1 <- populations[6]
      A2 <- populations[7]
            
      dS <- S * ((1 - config$m1) * (1 - config$m2) * monod(N, parms["S", c("mu", "k")]) -
        hill(A = A1, parms["S", c("psi", "phi1", "zeta1", "kappa1")]) -
        hill(A = A2, parms["S", c("psi", "phi2", "zeta2", "kappa2")]))
      dR1 <- R1 * ((1 - config$m2) * monod(N, parms["R1", c("mu", "k")]) -
        hill(A = A1, parms["R1", c("psi", "phi1", "zeta1", "kappa1")]) -
        hill(A = A2, parms["R1", c("psi", "phi2", "zeta2", "kappa2")])) +
        S * config$m1 * (1 - config$m2) * monod(N, parms["S", c("mu", "k")])
      dR2 <- R2 * ((1 - config$m1) * monod(N, parms["R2", c("mu", "k")]) -
        hill(A = A1, parms["R2", c("psi", "phi1", "zeta1", "kappa1")]) -
        hill(A = A2, parms["R2", c("psi", "phi2", "zeta2", "kappa2")])) +
        S * (1 - config$m1) * config$m2 * monod(N, parms["S", c("mu", "k")])
      dR12 <- R12 * (monod(N, parms["R12", c("mu", "k")]) -
        hill(A = A1, parms["R12", c("psi", "phi1", "zeta1", "kappa1")]) -
        hill(A = A2, parms["R12", c("psi", "phi2", "zeta2", "kappa2")])) +
        S * config$m1 * config$m2 * monod(N, parms["S", c("mu", "k")]) +
        R1 * config$m2 * monod(N, parms["R1", c("mu", "k")]) +
        R2 * config$m1 * monod(N, parms["R2", c("mu", "k")])
      dN <- deplete(dS,dR1,dR2,dR12,params[,"alpha"])
      dA1 <- -d1
      dA2 <- -d2
      return(list(c(dS, dR1, dR2, dR12, dN, dA1, dA2)))
    })
  }

  # Example usage
  times <- seq(0,50,0.01)
  populations <- c(S = 100, R1 = 10, R2 = 10, R12 = 0, N = 100,
    A1 = i1*config$pattern[1], A2 = i2*config$pattern[2])
  # events <- list(func = bottleneck, times = seq(10,40,10))
  events <- list(func = function(t, populations, parms, current_pattern) {
    res <- bottleneck(t = t, populations = populations, parms = parms, current_pattern = current_pattern)
    return(res$populations)
  }, times = seq(10, 40, 10))

# Run the simulation with the updated events list
solution <- ode(y = populations, times = times, func = bacterial_growth, events = events, parms = list(params, pattern = pattern))

  solution <- ode(y = populations, times = times, func = bacterial_growth, events = events, parms = config)
  bacteria_plot(solution)
}

simulate(stewardship = "2_only")
