library(tidyverse)

# Setting Parameters
population_size <- 1e2 # Number of individuals in the population
init_mutation_rate <- -2 # Since we are storing it in log10 format
init_fitness <- 1 # Represented by the net number of beneficial mutations
alpha <- 0.02 # Proportion of mutations that affect mutation rate
p_increase <- 0.5 # Proportion of mutations that increase mutation rate
mutation_jump_size <- 0.1 # Size of mutation rate change, increase or decrease
p_beneficial <- 0.1 # Proportion of mutations that increase fitness
generations <- 50 # Number of generations to simulate
s <- 0.1 # Parameter for fitness advantage/disadvantage

# Reproduction function
reproduce <- function(population) {
  population$individuals <- rpois(nrow(population), lambda = population$individuals * population$fitness)
  return(population)
}

# Mutation function
mutate <- function(offspring, alpha, p_increase, mutation_jump_size, p_beneficial) {
  # Calculate the number of mutants in each class
  offspring$mutants <- rbinom(nrow(offspring), offspring$individuals, 10^offspring$mutation_rate)
  offspring$non_mutants <- offspring$individuals - offspring$mutants
  # Calculate probabilities of each type of mutation
  p_vec <- c(alpha * c(p_increase, 1 - p_increase), (1 - alpha) * c(p_beneficial, 1 - p_beneficial))
  # Initialise the next generation with the non-mutants
  next_gen <- offspring %>%
    select(-c(individuals, mutants)) %>%
    rename(individuals = non_mutants)
  # Add rows for each mutant
  for (i in seq_len(nrow(offspring))) {
    mutate_into <- table(sample(c("a","b","c","d"), size = offspring$mutants[i], prob = p_vec))
    next_gen <- rbind(next_gen, data.frame(mutation_rate = offspring$mutation_rate[i] + mutation_jump_size * c(1, -1, 0, 0), fitness = offspring$fitness[i] + c(0, 0, 1, -1), individuals = mutate_into))
  }
  # Merge classes with same mutation_rate and fitness
  next_gen <- next_gen %>%
    group_by(mutation_rate, fitness) %>%
    summarize(individuals = sum(individuals), .groups = "drop")
  return(next_gen)
}

# Selection function
selection <- function(next_gen, population_size) {
  # Keep population size constant
  current_pop <- sum(next_gen$individuals)
  next_gen$individuals <- rbinom(length(next_gen), next_gen$individuals, min(1,population_size / current_pop))
  # Remove classes with zero individuals
  next_gen <- next_gen %>%
    filter(individuals > 0)
  return(next_gen)
}

# Initialize population
population <- data.frame(
  mutation_rate = init_mutation_rate,   # log10 of initial mutation rate
  fitness = init_fitness,    # initial fitness
  individuals = population_size    # initial population size
)
# Initialize data frame for statistics
stats <- data.frame(generation = integer(), avg_mutation_rate = numeric(), avg_fitness = numeric())

# Simulation loop
system.time(for (i in 1:generations) {
  # Step 1: Reproduction
  offspring <- reproduce(population)

  # Step 2: Mutation
  next_gen <- mutate(offspring, alpha, p_increase, mutation_jump_size, p_beneficial)

  # Step 3: Selection
  population <- selection(next_gen, population_size)

  # Save statistics for this generation
  stats <- rbind(stats, data.frame(generation = i, avg_mutation_rate = 10^(mean(population$mutation_rate)), avg_fitness_score = (1 + s)^mean(population$fitness_score)))
})

# Plot the results
library(ggplot2)

# Plot average mutation rate over generations
ggplot(stats, aes(x = generation)) +
  geom_line(aes(y = avg_mutation_rate), color = "blue") +
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
  geom_line(aes(y = avg_fitness_score), color = "red") +
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
