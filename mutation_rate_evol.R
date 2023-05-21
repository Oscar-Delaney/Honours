library(tidyverse)

# Reproduction function
reproduce <- function(pop, s) {
  pop$n <- rpois(nrow(pop), lambda = pop$n * (1 + s) ^ pop$w)
  return(pop)
}

do_mutations <- function(offspring, p_vec, jump) {
  # Calculate the total number of mutants and non-mutants
  total_mutants <- rbinom(nrow(offspring), offspring$n, 10^offspring$mu)
  non_mutants <- offspring$n - total_mutants

  # Generate counts for each mutation type
  mutation_counts <- cbind(t(non_mutants),matrix(rmultinom(nrow(offspring), size = total_mutants, prob = p_vec), nrow = nrow(offspring)))

  # Create a data frame for the next generation
  next_gen <- data.frame(
    mu = c(rep(offspring$mu, each = 5)) + jump * c(0, 1, -1, 0, 0),
    w = c(rep(offspring$w, each = 5)) + c(0, 0, 0, 1, -1),
    n = as.vector(t(mutation_counts))
  )
  
  # Merge classes with same mu and fitness
  next_gen <- next_gen %>%
    group_by(mu, w) %>%
    summarize(n = sum(n), .groups = "drop")
  return(next_gen)
}

# Selection function
selection <- function(next_gen, size) {
  # Keep pop size constant
  current_pop <- sum(next_gen$n)
  next_gen$n <- rbinom(nrow(next_gen), next_gen$n, min(1, size / current_pop))
  # Remove classes with zero n
  next_gen <- next_gen %>%
    filter(n > 0)
  return(next_gen)
}

evolve <- function(
  size = 1e2, # Number of individuals in the population
  init_mu = -2, # Since we are storing it in log10 format
  init_w = 0, # Represented by the net number of beneficial mutations
  a = 0.1, # Proportion of mutations that affect mutation rate
  p_mu_up = 0.5, # Proportion of mutations that increase mutation rate
  jump = 1e-1, # Size of mutation rate change, increase or decrease
  p_w_up = 0.1, # Proportion of mutations that increase fitness
  generations = 1e1, # Number of generations to simulate
  s = 1e-2 # Parameter for fitness advantage/disadvantage
) {
  # Initialize pop
  pop <- data.frame(
    mu = init_mu,   # log10 of initial mutation rate
    w = init_w,    # initial fitness
    n = size    # initial pop size
  )
  # Initialize data frame for statistics
  stats <- data.frame(generation = 0:generations, mu = init_mu, w = init_w)
  # Calculate probabilities of each type of mutation
  p_vec <- c(a * c(p_mu_up, 1 - p_mu_up), (1 - a) * c(p_w_up, 1 - p_w_up))
  print(p_vec)
  # Simulation loop
  for (i in 1:generations) {
    # Step 1: Reproduction
    offspring <- reproduce(pop, s)

    # Step 2: Mutation
    next_gen <- do_mutations(offspring, p_vec, jump)

    # Step 3: Selection
    pop <- selection(next_gen, size)

    # Save statistics for this generation
    stats[stats$generation == i, "mu"] <- weighted.mean(pop$mu, pop$n)
    stats[stats$generation == i, "w"] <- weighted.mean(pop$w, pop$n)
  }
  return(list(pop = pop, stats = stats))
}

# Save the pop and statistics
system.time({results <- evolve(size = 1e6, generations = 1e2, p_mu_up = 1, p_w_up = 0)})
pop <- results$pop
stats <- results$stats

# Plot average mutation rate over generations
ggplot(stats, aes(x = generation)) +
  geom_line(aes(y = mu), color = "blue") +
  labs(x = "Generation", y = "Average Mutation Rate",
       title = "Average Mutation Rate Over Generations") +
  theme_minimal() +
  theme(
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )

# Plot average fitness score over generations
ggplot(stats, aes(x = generation)) +
  geom_line(aes(y = w), color = "red") +
  labs(x = "Generation", y = "Average Fitness Score",
       title = "Average Fitness Score Over Generations") +
  theme_minimal() +
  theme(
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )

# Plot a heatmap of the pop
ggplot(pop, aes(x = mu, y = w)) +
  geom_tile(aes(fill = n)) +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Mutation Rate", y = "Fitness",
       title = "pop Heatmap") +
  theme_minimal() +
  theme(
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
