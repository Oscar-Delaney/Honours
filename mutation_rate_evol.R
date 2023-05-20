# Setting Parameters
population_size <- 1e3
mutation_rate <- -2 # Since we are storing it in log10 format
fitness_score <- 0
alpha <- 0.02
proportion_increase <- 0.5
mutation_jump_size <- 0.1
proportion_beneficial <- 0.1
generations <- 50
s <- 0.1 # Parameter for fitness advantage/disadvantage

# Initialize population
initialize_population <- function(population_size, mutation_rate, fitness_score) {
  population <- data.frame(
    mutation_rate = rep(mutation_rate, population_size),
    fitness_score = rep(fitness_score, population_size)
  )
  return(population)
}

# Reproduction function
reproduce <- function(population, s) {
  population$offspring_count <- rpois(nrow(population), lambda = (1 + s)^population$fitness_score)
  offspring <- population[rep(seq_len(nrow(population)), population$offspring_count), ]
  return(offspring)
}

# Mutation function
mutate <- function(offspring, mutation_rate, alpha, proportion_increase, mutation_jump_size, proportion_beneficial) {
  # Apply mutations
  for (i in 1:nrow(offspring)) {
    # Determine if mutation affects mutation rate or fitness
    if (runif(1) < alpha) {
      # Mutation affects mutation rate
      if (runif(1) < proportion_increase) {
        # Increase mutation rate
        offspring$mutation_rate[i] <- offspring$mutation_rate[i] + mutation_jump_size
      } else {
        # Decrease mutation rate
        offspring$mutation_rate[i] <- offspring$mutation_rate[i] - mutation_jump_size
      }
    } else {
      # Mutation affects fitness
      if (runif(1) < proportion_beneficial) {
        # Increase fitness
        offspring$fitness_score[i] <- offspring$fitness_score[i] + 1
      } else {
        # Decrease fitness
        offspring$fitness_score[i] <- offspring$fitness_score[i] - 1
      }
    }
  }
  return(offspring)
}

# Selection function
selection <- function(mutated_offspring, population_size) {
  # Keep population size constant
  if (nrow(mutated_offspring) > population_size) {
    selected <- mutated_offspring[sample(1:nrow(mutated_offspring), population_size), ]
  } else {
    selected <- mutated_offspring
  }
  return(selected)
}

# Initialize population
population <- initialize_population(population_size, mutation_rate, fitness_score)

# Initialize data frame for statistics
stats <- data.frame(generation = integer(), avg_mutation_rate = numeric(), avg_fitness_score = numeric())

# Simulation loop
for (i in 1:generations) {
  # Step 1: Reproduction
  offspring <- reproduce(population, s)
  
  # Step 2: Mutation
  mutated_offspring <- mutate(offspring, mutation_rate, alpha, proportion_increase, mutation_jump_size, proportion_beneficial)
  
  # Step 3: Selection
  population <- selection(mutated_offspring, population_size)
  
  # Save statistics for this generation
  stats <- rbind(stats, c(i, 10^(mean(population$mutation_rate)), (1 + s)^mean(population$fitness_score)))
}

# Plot the results
library(ggplot2)

# Plot average mutation rate over generations
ggplot(stats, aes(x = generation)) +
  geom_line(aes(y = avg_mutation_rate), color = "blue") +
  labs(x = "Generation", y = "Average Mutation Rate", 
       title = "Average Mutation Rate Over Generations") +
  theme_minimal()

# Plot average fitness score over generations
ggplot(stats, aes(x = generation)) +
  geom_line(aes(y = avg_fitness_score), color = "red") +
  labs(x = "Generation", y = "Average Fitness Score", 
       title = "Average Fitness Score Over Generations") +
  theme_minimal()
