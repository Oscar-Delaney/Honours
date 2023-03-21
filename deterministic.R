library(deSolve)

D = 0.1 #bottleneck dilution ratio
N0 = 100 #initial nutrient concentration

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

hill <- function(A, parms) {
  psi = parms[1] ; phi = parms[2] ; zeta = parms[3]; kappa = parms[4]
  return(phi*(A/zeta)^kappa/((A/zeta)^kappa-(psi-phi)/psi))
}

monod <- function(N,parms){
  mu = parms[1] ; k = parms[2]
  return(mu*N/(N+k))
}

# a growth rate function using hill and monod functions
growth <- function(A1,A2,N,parms) {
  return(monod(N,parms[8:9])-hill(A1,parms[1:4])-hill(A2,parms[c(1,5:7)]))
}

# antibiotics <- function(A1,A2,parms) {
#   psi = parms[1]
#   return(psi-hill(A1,c(psi,parms[2:4]))-hill(A2,c(psi,parms[5:7])))
# }



deplete <- function(dS=0,dR1=0,dR2=0,dR12=0,alpha=params[,"alpha"]) {
  return(sum(-alpha*max(0,c(dS,dR1,dR2,dR12))))
}

# Define a bottleneck function that reduces all populations by a fixed fraction
bottleneck <- function(t, populations, parms = params, fraction = D, nutrient = N0) {
  populations <- c(populations[-5]*fraction, nutrient)
  populations[populations < 0] <- 0
  return(populations)
}


# Define the differential equations for the model
bacterial_growth <- function(t, populations, parms = params) {
  with(as.list(c(populations)), {
    S <- populations[1]
    R1 <- populations[2]
    R2 <- populations[3]
    R12 <- populations[4]
    N <- populations[5]
    dS <- S*(growth(A1=1,A2=0,N,params["S",1:9]) - params["S","m_rate"])
    dR1 <- R1*(growth(A1=1,A2=0,N,params["R1",1:9])) + (S-R1)*params["R1","m_rate"]
    dR2 <- R2*(growth(A1=1,A2=0,N,params["R2",1:9])) + (S-R2)*params["R2","m_rate"]
    dR12 <- R12*(growth(A1=1,A2=0,N,params["R12",1:9])) + sum(c(R1,R2)*params[c("R1","R2"),"m_rate"])
    dN <- deplete(dS,dR1,dR2,dR12,params[,"alpha"])
    return(list(c(dS, dR1, dR2, dR12, dN)))
  })
}

# Example usage
times <- seq(0,50,0.01)
populations <- c(S = 100, R1 = 10, R2 = 10, R12 = 0, N = 100)
event_times <- seq(10,40,10)
events <- list(func = bottleneck, times = event_times)
solution <- ode(y = populations, times = times, func = bacterial_growth, events = events)
plot(solution[, 1], solution[,2],type = "l", col = 1, main = "Bacterial growth over time", xlab = "Time", ylab = "Population size")
lines(solution[, 1],solution[, 3], col = 2)
lines(solution[, 1],solution[, 4], col = 3)
lines(solution[, 1],solution[, 5], col = 4)
lines(solution[, 1],solution[, 6], col = 5)
legend("topright", legend = c("Susceptible", "Resistant1", "Resistant2", "Resistant12", "Nutrient"), col = c(1, 2, 3, 4, 5), lty = 1)
