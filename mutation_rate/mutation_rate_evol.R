library(tidyverse)

# Reproduction rate function
rate <- function(w, s, r, r_max) {
  return(1 / ((((1 + s) ^ -w) / r) + 1 / r_max))
}
# Reproduction function
reproduce <- function(pop, s, r, r_max, det) {
  growth <- rate(pop$w, s, r, r_max)
  if (det) {
    n <- pop$n * growth
    pop$n <- n / sum(n)
  } else {
  pop$n <- rpois(nrow(pop), lambda = pop$n * growth)
  }
  return(pop)
}

do_mutations <- function(offspring, p_vec, jump, det) {
  # Calculate the total number of mutants and non-mutants
  if (det) {
    mutants <- offspring$n * pmin(1, 10^offspring$mu)
    non_mutants <- offspring$n - mutants
    m_counts <- outer(p_vec, mutants, "*")
  } else {
    mutants <- rbinom(nrow(offspring), offspring$n, pmin(1, 10^offspring$mu))
    non_mutants <- offspring$n - mutants
    m_counts <- sapply(mutants, function(x) rmultinom(1, size = x, prob = p_vec))
  }

  # Generate counts for each mutation type, and combine with non-mutants
  class_counts <- rbind(non_mutants, m_counts)

  # Create a data frame for the next generation
  next_gen <- data.frame(
    mu = round(rep(offspring$mu, each = 5) + jump * c(0, 1, -1, 0, 0),2),
    w = rep(offspring$w, each = 5) + c(0, 0, 0, 1, -1),
    n = as.vector(class_counts)
  )

  # Merge classes with same mu and fitness
  next_gen <- next_gen %>%
    group_by(mu, w) %>%
    summarize(n = sum(n), .groups = "drop")
  return(next_gen)
}

# Selection function
selection <- function(next_gen, size, sensitivity, det) {
  if (!det) {
    # Keep pop size constant
    current_pop <- sum(next_gen$n)
    next_gen$n <- rbinom(nrow(next_gen), next_gen$n, min(1, size / current_pop))
  }
  # Remove classes with zero n
  next_gen <- next_gen %>%
    filter(n > sensitivity)
  return(next_gen)
}

# Function to run the simulation
evolve <- function(
  size = 1e2, # Number of individuals in the population
  init_mu = -2, # Since we are storing it in log10 format
  init_w = 0, # Represented by the net number of beneficial mutations
  a = 0.1, # Proportion of mutations that affect mutation rate
  p_mu_up = 0.5, # Proportion of mutations that increase mutation rate
  jump = 1e-1, # Size of mutation rate change, increase or decrease
  p_w_up = 0.1, # Proportion of mutations that increase fitness
  generations = 1e1, # Number of generations to simulate
  s = 1e-2, # Parameter for fitness advantage/disadvantage
  r = 1, # Parameter for baseline reproduction rate
  r_max = 1e5, # Parameter for maximum reproduction rate
  sensitivity = 1e-10, # Parameter for minimum detectable population size
  det = FALSE # Whether the model is deterministic
) {
  # Initialize pop
  pop <- data.frame(
    mu = init_mu,   # log10 of initial mutation rate
    w = init_w,    # initial fitness
    n = size ^ !det    # initial pop size
  )
  # Initialize data frame for statistics
  stats <- data.frame(generation = 0:generations, mu = mean(init_mu), w = init_w)
  # Calculate probabilities of each type of mutation
  p_vec <- c(a * c(p_mu_up, 1 - p_mu_up), (1 - a) * c(p_w_up, 1 - p_w_up))

  # Simulation loop
  for (i in 1:generations) {
    # Step 1: Reproduction
    offspring <- reproduce(pop, s, r, r_max, det)

    # Step 2: Mutation
    next_gen <- do_mutations(offspring, p_vec, jump, det)

    # Step 3: Selection
    pop <- selection(next_gen, size, sensitivity, det)
    if (nrow(pop) == 0) {
      print("Population went extinct!")
      break
    }

    # Save statistics for this generation
    stats[stats$generation == i, "mu"] <- weighted.mean(pop$mu, pop$n)
    stats[stats$generation == i, "w"] <- weighted.mean(pop$w, pop$n)
  }

  return(list(pop = pop, stats = stats))
}

# Function to graph mutation rate or fitness
plot_univariate <- function(stats, var, y_label) {
  ggplot(stats, aes(x = generation)) +
    geom_line(aes(y = get(var)), color = "blue") +
    labs(x = "Generations", y = y_label) +
    theme_minimal() +
  theme(
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
    )
}

# Function to plot heatmap of population
plot_heatmap <- function(pop) {
  ggplot(pop, aes(x = mu, y = w)) +
    geom_tile(aes(fill = n)) +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(x = "Mutation Rate", y = "Fitness",
         title = "Final Population Distribution") +
    theme_minimal() +
  theme(
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
}

# Save the pop and statistics
system.time({results <- evolve(generations = 1e3, init_mu = seq(-2,-1,by=0.1), sensitivity = 1e-25, r_max = Inf, det = T, size = 1e6, r = 1, s = 0.1, jump = 0, p_mu_up = 0.5, p_w_up = 1e-6)})
pop <- results$pop
stats <- results$stats
tail(stats)

# Plot the results
plot_univariate(stats, "mu", "log_10 Average Mutation Rate")
plot_univariate(stats, "w", "Average Fitness Score")
plot_heatmap(pop)
