# Generate and analyse data
run_sims <- function(summary, rep = 1, time = 50, s = 0.1, r = 1, res = TRUE, func = metric) {
    summary$m1 <- summary$D ^ - 0.5 * 1e-9
    data <- list()
    for (i in seq_len(nrow(summary))) {
        D <- summary$D[i]
        m1 <- summary$m1[i]
        tau <- summary$tau[i]
        k_ratio <- res * ifelse("k_ratio" %in% names(summary), summary$k_ratio[i], 1)
        N0 <- 1e9
        if (D >= exp(-r * (1 + k_ratio) * tau) ) {
            data <- simulate(
                seed = i,
                rep = rep,
                time = time,
                dt = 1e-1,
                tau = tau,
                max_step = Inf,
                D = D,
                N0 = N0,
                k = N0 * k_ratio,
                alpha = 1 * res,
                r = 1.023 * r * (1 + k_ratio), # Adaptivetau step size causes observed growth rate to be lower than expected
                s = s,
                m1 = m1,
                init_W = round(N0 * D),
                num_mutants = 1e2,
                # loci = loci,
            )
            summary[i, c("rate", "ci_lower", "ci_upper")] <- func(data)
        }
        print(i / nrow(summary))
    }
    return(summary)
}

# Create your data frame
df <- data.frame(
  id = c(1, 2, 3),
  x = c(1, 2, 3)
)
# Let's say we want to assign a list to the first row of column z
df$z <- vector("list", nrow(df)) # create a list column of the same length as the data frame
df$z[[1]] <- 1:10 # assign a list to the first cell


# Create your list of lists
list_of_lists <- list(seq(1,10), seq(11,20), seq(21,30))

# Assign the list to a new column
df$y <- list_of_lists

# Print the data frame
print(df)

print(df)

s_dist <- function(data) {
    return(list(list("A", "B", "C"), list("A", "B", "C"), 2))
}

# Probability a new mutant at the beginning of a growth phase will go extinct
phi <- function(D, s) {
    return(ifelse(D == 1, 1 / (1 + s), (1 - D) / (D ^ -s - D)))
}

# Hypergeometric function
hyper <- function(s, z) {
    return(Re(hypergeo(1, 1 / (1 + s), 1 + 1 / (1 + s), z)))
}

# theoretical rate function
theory <- function(D, r, s) {
    return(-r / log(D) * (hyper(s, (1 - D) / (D * (D ^ s - 1))) -
        D * hyper(s, (1 - D) * (D ^ s) /  (D ^ s - 1))))
}

# Approximate rate function
approx1_theory <- function(D, r, s) {
    return(r * s * log(D^-1) / (D^-1 -1))
}

# Even more approximate rate function
approx2_theory <- function(D, r, s) {
    return(r * s * sqrt(D))
}

# evaluate a set of simulation results
metric <- function(data) {
    # find the time just after the last bottleneck
    endpoint <- data[[1]] %>%
        filter(variable == "W", rep == 1) %>%
        filter(value - lag(value) < 0) %>%
        tail(1) %>%
        pull(time)
    # find the likelihood of a new mutation at t=0 going extinct
    current_phi <- phi(data[[2]]$D, data[[2]]$s)
    # count the number of each mutant at the endpoint
    final_counts <- data[[1]] %>%
        group_by(rep, variable) %>%
        filter(time == endpoint, !(variable %in% c("W", "N"))) %>%
        summarise(final_value = value, .groups = "keep") %>%
        mutate(p_fix = 1 - current_phi ^ final_value)
    # estimate the number of these mutations that will go on to fix
    fixed <- final_counts %>%
        group_by(rep) %>%
        summarise(n = sum(final_value > 1e1), n_hat = sum(p_fix))
    # estimate the fixation rate and store this
    fixation_rate <- fixed$n_hat / endpoint /
        ((data[[2]]$init_W / data[[2]]$D) * data[[2]]$m1)
    se <- sd(fixation_rate) / sqrt(length(fixation_rate))
    ci <- mean(fixation_rate) + se * qnorm(c(0.5, 0.025, 0.975))
    return(ci)
}

# evaluate a set of simulation results
metric_ci_old <- function(data) {
    # extract some relevant parameters
    m1 <- data[[2]]$m1
    wt_max <- data[[2]]$init_W / data[[2]]$D
    wt <- data[[2]]$names[1]
    loci <- data[[2]]$loci
    current_phi <- phi(data[[2]]$D, data[[2]]$s)
    # find the time just after the last bottleneck
    endpoint <- data[[1]] %>%
        filter(variable == wt, rep == 1) %>%
        filter(value - lag(value) < 0) %>%
        tail(1) %>%
        pull(time)
    # find the number of each genotype at the endpoint
    final <- data[[1]] %>%
        filter(time == endpoint, variable != "N")
    # Note which mutations each genotype has
    for (i in 1:loci) {
        final[[paste0("M", i)]] <- substring(final$variable, i+1, i+1)=="1"
    }
    # Calculate the abundance of each mutation
    counts <- final %>%
        pivot_longer(
            cols = starts_with("M"),
            names_to = "Mutant",
            values_to = "Mutant_value"
        ) %>%
        group_by(rep, Mutant) %>%
        summarise(final_value = sum(value * Mutant_value), .groups = "keep")
    # estimate the fixation rate and store this
    mean <- mean(current_phi ^ counts$final_value)
    ci <- binom.test(x = round(nrow(counts) * c(mean, 1 - mean)))$conf.int
    vec <- c(mean, rev(ci))
    # Note the inverse cdf of the exponential distribution is -log(1-F(x))/lambda
    fixation_rate <- - log(vec) / (endpoint * wt_max * m1 / loci)
    return(fixation_rate)
}

metric_ci <- function(data) {
    # extract some relevant parameters
    m1 <- data[[2]]$m1
    wt_max <- data[[2]]$init_W / data[[2]]$D
    loci <- data[[2]]$loci
    time <- data[[2]]$time
    # find the number of each genotype at the endpoint
    final <- data[[1]] %>%
        filter(time == time, variable != "N")
    # Note which mutatiosn each genotype has
    for (i in 1:loci) {
        final[[paste0("M", i)]] <- substring(final$variable, i+1, i+1)=="1"
    }
    # Calculate the abundance and probability  for each mutation
    counts <- final %>%
        pivot_longer(
            cols = starts_with("M"),
            names_to = "Mutant",
            values_to = "Mutant_value"
        ) %>%
        group_by(rep, Mutant) %>%
        summarise(p = sum(value * Mutant_value) / sum(value), .groups = "keep")
    # estimate the fixation rate and store this
    se <- sd(counts$p) / sqrt(length(counts$p))
    mean <- mean(counts$p)
    vec <- mean + se * qnorm(c(0.5, 0.025, 0.975))
    return(vec)
}

mutation_time <- function(data) {
    # Assuming final_counts and mutation_times have been computed as before...
    joint_data <- final_counts %>%
        inner_join(mutation_times, by = c("rep", "variable"))

    # Extract your data and corresponding probabilities
    X <- joint_data$last_zero
    P <- joint_data$p_fix

    # Normalize the weights
    weights <- P / sum(P)

    # Load the necessary library
    library(KernSmooth)

    # Calculate bandwidth using Wand's rule
    bw <- dpik(X, weights=weights)

    # Compute the density estimate
    kde <- bkde(X, weight=weights, bandwidth=bw)

    # Plot the KDE
    p <- plot(kde, main='Kernel Density Estimate of Mutation Times')
    return(p)
}

s_dist <- function(data) {
    
}