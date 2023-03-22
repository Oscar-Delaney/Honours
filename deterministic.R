library(deSolve)

pharmacokinetic <- FALSE # should be either TRUE or FALSE
stewardship <- "cycl" # "cycl" or "comb" or "1_only" or "2_only"

# Define the parameters of the model
D <- 0.1 # bottleneck dilution ratio
N0 <- 100 # initial nutrient concentration
m1 <- 0.1 # rate of mutations conferring resistance to drug 1
m2 <- 0.1 # rate of mutations conferring resistance to drug 2
d1 <- 1*pharmacokinetic # rate of drug 1 elimination
d2 <- 1*pharmacokinetic # rate of drug 2 elimination
i1 <- 10 # drug 1 influx concentration
i2 <- 10 # drug 2 influx concentration
if (stewardship=="2_only") {pattern <- c(0,1)
  } else if (stewardship=="comb") {pattern <- c(1,1)
  } else {pattern <- c(1,0)
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

# Define a bottleneck function that reduces all populations by a fixed fraction
bottleneck <- function(t, populations, parms = params) {
  populations <- c(populations[1:4]*D, N0, pharmacokinetic*populations[6:7]+c(i1,i2)*pattern)
  # names(populations) <- c("S", "R1", "R2", "R12", "N", "A1", "A2")
  populations[populations < 0] <- 0 # backup - shouldn't be needed
  if (stewardship=="cycl") {pattern <<- 1 - pattern}
  return(populations)
  }

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
bacterial_growth <- function(t, populations, parms = params) {
  with(as.list(c(populations)), {
    S <- populations[1]
    R1 <- populations[2]
    R2 <- populations[3]
    R12 <- populations[4]
    N <- populations[5]
    A1 <- populations[6]
    A2 <- populations[7]
    dS <- S*((1-m1)*(1-m2)*monod(N,params["S",8:9]) -
      hill(A=A1,params["S",1:4]) - 
      hill(A=A2,params["S",c(1,5:7)]))
    dR1 <- R1*((1-m2)*monod(N,params["R1",8:9]) -
      hill(A=A1,params["R1",1:4]) - 
      hill(A=A2,params["R1",c(1,5:7)])) +
      S*m1*(1-m2)*monod(N,params["S",8:9])
    dR2 <- R2*((1-m1)*monod(N,params["R2",8:9]) -
      hill(A=A1,params["R2",1:4]) - 
      hill(A=A2,params["R2",c(1,5:7)])) +
      S*(1-m1)*m2*monod(N,params["S",8:9])
    dR12 <- R12*(monod(N,params["R12",8:9]) -
      hill(A=A1,params["R12",1:4]) - 
      hill(A=A2,params["R12",c(1,5:7)])) +
      S*m1*m2*monod(N,params["S",8:9]) +
      R1*m2*monod(N,params["R1",8:9]) +
      R2*m1*monod(N,params["R2",8:9])
    dN <- deplete(dS,dR1,dR2,dR12,params[,"alpha"])
    dA1 <- -d1
    dA2 <- -d2
    return(list(c(dS, dR1, dR2, dR12, dN, dA1, dA2)))
  })
}

# Example usage
times <- seq(0,50,0.01)
populations <- c(S = 100, R1 = 10, R2 = 10, R12 = 0, N = 100,
  A1 = i1*pattern[1], A2 = i2*pattern[2])
events <- list(func = bottleneck, times = seq(10,40,10))
solution <- ode(y = populations, times = times, func = bacterial_growth, events = events)
bacteria_plot(solution)
