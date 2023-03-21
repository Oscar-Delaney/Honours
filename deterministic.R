library(deSolve)

Psi = 0.1
Phi = 0.2
Zeta = 1
Kappa = 1

hill <- function(psi=Psi,A,phi=Phi,zeta=Zeta,kappa=Kappa) {
  return(phi*(A/zeta)^kappa/((A/zeta)^kappa-(psi-phi)/psi))
}

fitness <- function(psi=Psi,A1,phi1=Phi,zeta1=Zeta,kappa1=Kappa,A2,phi2=Phi,zeta2=Zeta,kappa2=Kappa) {
  return(psi-hill(psi,A1,phi1,zeta1,kappa1)-hill(psi,A2,phi2,zeta2,kappa2))
}

growth <- function(mu,R,k){
  return(mu*R/(R+k))
}

print(hill(0,2,1,1,1))
print(fitness(A1=100,A2=1,zeta2=0.01))


# Function to simulate bacterial growth over time
simulate_growth <- function(initial_populations, params, time) {
  
  # Define the differential equations for the model
  bacterial_growth <- function(t, populations, params) {
    with(list(S = populations[1], R = populations[2], rmax_S = params[1], rmax_R = params[2],
              K = params[3], beta = params[4], delta = params[5]), {
                dS <- rmax_S * S * (1 - ((S + R) / K)) - beta * S * R
                dR <- rmax_R * R * (1 - ((S + R) / K)) + beta * S * R - delta * R
                return(list(c(dS, dR)))
              })
  }
  
  # Set the initial populations and parameters
  S0 <- initial_populations[1]
  R0 <- initial_populations[2]
  rmax_S <- params[1]
  rmax_R <- params[2]
  K <- params[3]
  beta <- params[4]
  delta <- params[5]
  
  # Set the times at which to compute the solutions
  times <- seq(from = 0, to = time, by = 0.01)
  
  # Solve the differential equations using deSolve
  sol <- ode(y = c(S = S0, R = R0), times = times, func = bacterial_growth, parms = params)
  
  # Return the solution
  return(sol)
}

# Example usage
initial_populations <- c(S = 100, R = 10)
params <- c(rmax_S = 0.8, rmax_R = 0.6, K = 1000, beta = 0.01, delta = 0.05)
time <- 50
solution <- simulate_growth(initial_populations, params, time)
plot(solution[,1], solution[,3], type = "l", col = 2, main = "Bacterial growth over time", xlab = "Time", ylab = "Population size")
lines(solution[,1], solution[,2], col = 1)
legend("topright", legend = c("Susceptible", "Resistant"), col = c(1, 2), lty = 1)
