library(deSolve)

q <- 0

# Define the ODE system
ode_fun <- function(t, y, parms) {
  dy_dt <- -y # y' = -y
  return(list(dy_dt))
}

# Define the event function
event_fun <- function(t, y, parms) {
  print(t)
  q <<- q +1
  return(y <- y + 1)
}

# Set up initial conditions and time points
y0 <- 0.5 # initial value of y
times <- seq(0, 5, by = 0.1) # time points to solve ODEs

# Set up the parameters for the ODE system
parms <- NULL

# Solve the ODE system with the event function
ode_sol <- ode(y = y0, times = times, func = ode_fun, parms = parms,
               events = list(func = event_fun, times = c(2)))

# Plot the solution
plot(ode_sol, xlab = "time", ylab = "y")
print(q)
