library(adaptivetau)
library(deSolve)
library(tidyverse)
library(ggnewscale)
library(future)
library(future.apply)
library(R.utils)

# a growth rate function for nutrient-limited growth
monod <- function(R, r, k) {
  return(r * ifelse(k == 0, 1, 1 / (1 + k / R)))
}

# fitness advantages of new mutants
mutant_fitness <- function(num_mutants, r, w, names) {
  return(setNames(r * (1 + c(0, rexp(num_mutants, 1 / w))), names))
}

# a function that reduces all populations by a factor of D, in expectation
bottleneck <- function(state, config) {
  with(config, {
    pops <- setNames(rbinom(length(names), state[names], D), names)
    R <- state["R"] * D + R0 * (1 - D)
    return(c(pops, R))
  })
}

# a function outputting the transitions that can occur in the model
make_transitions <- function(names) {
    # Create the initial list
    rate_list <- list()
    # Add the growth rates for each mutant
    for (i in names) {
      rate_list[[paste0(i, "_growth")]] <- setNames(c(+1), i)
    }
    # Add the nutrient depletion rate
    rate_list[["N_depletion"]] <- setNames(c(-1), "R")
  return(rate_list)
}

# compute the rate at which each transition occurs
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # Calculate replication rates
    replication_rates <- state[names] * monod(R, r_vec, k)
    if (is.numeric(num_mutants)) {
        # find which mutant category we should be filling
        to_mutate <- min(which(state == 0), num_mutants + 2) - 1
        # chance of a replication in row i resulting in a strain j cell
        mutation <- diag(num_mutants + 1)
        if (to_mutate <= num_mutants) {
          mutation[1, to_mutate + 1] <- mu
          mutation[1, 1] <- 1 - mu
        } else {
          print("Out of mutant spots! :(")
        }
    }
    # Calculate growth rates including mutations
    growth_rates <- replication_rates %*% mutation
    # Calculate nutrient depletion rate
    N_depletion <- sum(replication_rates * alpha)
    # Combine all rates and return
    rate_list <- c(growth_rates, N_depletion)
    return(setNames(rate_list, c(names, "R")))
  })
}

# a function to implement one run of the model
single_run <- function(config, x) {
  with(config, {
    # Define the transitions of the model
    transitions <- make_transitions(names)
    # Initialise the state variables
    state <- init
    if (is.numeric(num_mutants)) {
      config$r_vec <- mutant_fitness(num_mutants, r, w, names)
    } else {
      config$r_vec <- r * t(1 + as.matrix(genotypes) %*% rexp(loci, 1 / w))
    }
    time_grid <- seq(0, time, by = dt) # a common time grid for all runs
    bottlenecks <- unique(round(c(seq(0, time, tau), time), 10))
    for (t in bottlenecks[-length(bottlenecks)]) {
      # Determine the time until the next bottleneck or dose
      end <- min(bottlenecks[bottlenecks > t] - t)
      # Create a new vector of growth rates for mutants, if locus-agnostic
      if (is.numeric(num_mutants)) {
        config$r_vec <- ifelse(state[0:num_mutants + 1] == 0,
          mutant_fitness(num_mutants, r, w, names), config$r_vec)
      }
      # set the seed for reproducibility
      if (is.numeric(seed)) set.seed(round(seed + (x * time + t) / tau))
      # Run the model between bottlenecks
      new <- ssa.adaptivetau(
        state, transitions, rates, config, tf = end,
        tl.params = list(maxtau = max_step),
        deterministic = grep("depletion", names(transitions))
      )
      # Make the time column reflect the overall time accurately
      new[, "time"] <- new[, "time"] + t
      # Avoid duplicate times by adding a small increment after the bottleneck
      new[1, "time"] <- new[1, "time"] * (1 + 1e-6)
      # Update the solution
      solution <- if (t == 0) new else rbind(solution, new)
      # Run the bottleneck and update the state
      state <- bottleneck(new[nrow(new), ], config)
    }
    # Interpolate the solution to the common time grid
    approx_vars <- lapply(colnames(solution), function(var) {
      approx(solution[, "time"], solution[, var], xout = time_grid)$y
    })
    solution_interpolated <- data.frame(
      setNames(approx_vars, colnames(solution)),
      rep = x
    )
    # add a new row at the end with the r_vec
    solution_interpolated[nrow(solution_interpolated) + 1, ] <-
      c(pi * 1e6, (config$r_vec / r) - 1, 0, x)
    return(solution_interpolated)
  })
}

# a function to simulate the model
simulate <- function(
  rep = 1, # number of runs of the simulation
  seed = NULL, # seed for reproducibility
  time = 100, # time to simulate, in hours
  dt = 0.1, # time step, in hours
  max_step = Inf, # SSA max step parameter
  tau = 3, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  R0 = 1e9, # initial nutrient concentration
  mu = 1e-9, # rate of mutations conferring resistance to drug 1
  N = 1e9, # (1/D) of the initial wild type population
  num_mutants = NULL, # number of mutants
  loci = NULL, # number of mutable loci in the genome
  r = 1, # wild type growth rate with infinite resources
  w = 0.1, # mean fitness effect size of beneficial mutation
  k = 1e8, # [R] at half-max growth rate
  alpha = 1 # nutrients used per replication
  ) {
  # Define the parameters of the model
  if (is.numeric(num_mutants)) {
    names <- c("W", paste0("M", 1:num_mutants))
  } else {
      if (is.null(loci)) stop("Must specify num_mutants XOR loci")
      names <- paste0("G", intToBin(1:2^loci - 1))
      # Generate all possible binary strings of length loci
      genotypes <- expand.grid(replicate(loci, c(0, 1), simplify = FALSE))
      genotypes <- setNames(rev(genotypes), paste0("M", 1:loci))
      # Initialize the mutation matrix
      mutation <- diag(1 - mu, 2 ^ loci)
      # Loop over all genotypes
      for (i in 1:2^loci) {
        for (j in 1:2^loci) {
          # If the Hamming distance is 1, set the mutation rate to be nonzero
          if (sum(genotypes[i, ] != genotypes[j, ]) == 1) {
            mutation[i, j] <- mu / loci
          }
        }
      }
  }
  init <- setNames(c(round(N * D), rep(0, length(names) - 1), R0), c(names, "R"))
  config <- as.list(environment())
  # Run the simulation rep number of times, using parallelisation if possible
  plan(multisession) # compatible with both unix and Windows
  set.seed(seed) # set the seed for reproducibility
  solutions <- bind_rows(future_lapply(1:rep, function(x) {
    single_run(config, x)},
    future.seed = TRUE))
  # Convert the solutions to long format
  long <- pivot_longer(solutions, cols = -c(time, rep), names_to = "variable")
  config$s_all <- long %>% filter(near(time, pi * 1e6)) %>% select(-time)
  long <- long[!near(long$time, pi * 1e6), ]
  return(list(long, config))
}

log_plot <- function(solutions, type = "all") {
  # Choose the type of central tendency and range to plot
  if (type == "mean") {
    summary <- solutions %>%
      group_by(time, variable) %>%
      reframe(
        central = mean(value),
        se = sd(value) / sqrt(n()),
        lower = max(0, central - 1.96 * se),
        upper = central + 1.96 * se,
        rep = 1
      )
  } else if (type == "median") {
    summary <- solutions %>%
      group_by(time, variable) %>%
      reframe(
        IQR_bounds = list(quantile(value, c(0.25, 0.5, 0.75))),
        rep = 1
      ) %>%
      # convert the list of quantiles to individual columns
      mutate(
        lower = map_dbl(IQR_bounds, 1),
        central = map_dbl(IQR_bounds, 2),
        upper = map_dbl(IQR_bounds, 3)
      ) %>%
      select(-IQR_bounds)
  } else if (type == "all") {
    summary <- solutions
    colnames(summary)[colnames(summary) == "value"] <- "central"
  } else {
    stop("type must be 'all', 'median' or 'mean'")
  }
  filtered <- summary %>%
    filter(!(variable %in% c("R"))) %>%
    mutate(variable = factor(variable, levels = unique(variable))) %>%
    arrange(variable)
  # Initialise the colours
  colors <- c("black", hcl.colors(length(levels(filtered$variable)) - 1, "Dark 2"))
  peak <- max(solutions$value, na.rm = TRUE)
  # Create the plot
  plot <- ggplot() +
    geom_line(data = filtered, aes(x = time, y = central, color = variable,
      linetype = factor(rep)), linewidth = 1 + 0.5 / max(filtered$rep)) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = rep("solid", max(filtered$rep))) +
    guides(linetype = "none") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      breaks = 10^seq(0, 20),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(
      title = "Bacterial growth over time",
      x = "Time (hours)",
      y = "Population Size",
      color = "Strain",
      fill = "Strain"
    ) +
    theme_light() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
  # Add the confidence intervals
  if (type %in% c("mean", "median") && max(solutions$rep) > 1) {
    plot <- plot +
      geom_ribbon(data = filtered, alpha = 0.3,
        aes(x = time, ymin = lower, ymax = upper, fill = variable)) +
      scale_fill_manual(values = colors)
  }
  # Display the plot
  print(plot)
}

# Generate and analyse data
run_sims <- function(summary, rep = 1, time = 50, w = 0.1, r = 1, mu = 1e-9,
    res = TRUE, num_mutants = 1e2, loci = NULL) {
    func <- ifelse(is.null(loci), metric, metric_ci)
    data <- list()
    for (i in seq_len(nrow(summary))) {
        D <- summary$D[i]
        tau <- summary$tau[i]
        k_ratio <- res *
            ifelse("k_ratio" %in% names(summary), summary$k_ratio[i], 1)
        data <- simulate(
            seed = i,
            rep = rep,
            time = time,
            dt = 1e-1,
            tau = tau,
            max_step = Inf,
            D = D,
            R0 = 1e9,
            k = 1e9 * k_ratio,
            alpha = 1 * res,
            r = r * (1 + k_ratio), # Adaptivetau step size
            # causes observed growth rate to be lower than expected
            w = w,
            mu = mu,
            N = 1e9,
            num_mutants = num_mutants,
            loci = loci
        )
        summary[i, c("rate", "ci_lower", "ci_upper")] <- func(data)
        print(i / nrow(summary))
    }
    return(summary)
}


# Probability a new mutant at the beginning of a growth phase will go extinct
phi <- function(D, s) {
    return(ifelse(D == 1, 1 / (1 + s), (1 - D) / (D ^ -s - D)))
}

# Hypergeometric function
hyper <- function(w, z) {
    return(Re(hypergeo(1, 1 / (1 + w), 1 + 1 / (1 + w), z)))
}

# theoretical rate function
theory <- function(D, r, w) {
    return(-r / log(D) * (hyper(w, (1 - D) / (D * (D ^ w - 1))) -
        D * hyper(w, (1 - D) * (D ^ w) /  (D ^ w - 1))))
}

# Approximate rate function
approx1_theory <- function(D, r, w) {
    return(r * w * log(D^-1) / (D^-1 -1))
}

# Even more approximate rate function
approx2_theory <- function(D, r, w) {
    return(r * w * sqrt(D))
}

# rate at a given time within [0, tau)
rate_at_t <- function(D, r, w, t) {
  (D * (-1 + D^w) * exp(r * t) * r) /
    (-1 + D^(1 + w) * exp(r * (1 + w) * t) - D^w * (-1 + exp(r * (1 + w) * t)))
}

# Count the number of mutants likely en route to fixation
fixed <- function(data) {
    # find the time just after the last bottleneck
    endpoint <- data[[1]] %>%
        filter(variable == "W", rep == 1) %>%
        filter(value - lag(value) < 0) %>%
        tail(1) %>%
        pull(time)
    # find the likelihood of a new mutation at t=0 going extinct
    current_phi <- phi(data[[2]]$D, data[[2]]$w)
    # count the number of each mutant at the endpoint
    final_counts <- data[[1]] %>%
        group_by(rep, variable) %>%
        filter(time == endpoint, !(variable %in% c("W", "R"))) %>%
        summarise(final_value = value, .groups = "keep") %>%
        mutate(p_fix = 1 - current_phi ^ final_value)
    return(list(final_counts, endpoint))
}

# evaluate a set of simulation results
metric <- function(data) {
    input <- fixed(data)
    final_counts <- input[[1]]
    endpoint <- input[[2]]
    # estimate the number of these mutations that will go on to fix
    fixed <- final_counts %>%
        group_by(rep) %>%
        summarise(n = sum(final_value > 1e1), n_hat = sum(p_fix))
    # estimate the fixation rate and store this
    fixation_rate <- fixed$n_hat / endpoint
    se <- sd(fixation_rate) / sqrt(length(fixation_rate))
    ci <- mean(fixation_rate) + se * qnorm(c(0.5, 0.025, 0.975))
    return(ci)
}

metric_ci <- function(data) {
    # find the number of each genotype at the endpoint
    final <- data[[1]] %>%
        filter(time == max(time), variable != "R")
    # Note which mutations each genotype has
    for (i in 1:data[[2]]$loci) {
        final[[paste0("M", i)]] <- substring(final$variable, i + 1, i + 1) == "1"
    }
    # Calculate the abundance and probability for each mutation
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

r_res <- 1.5
summary <- expand.grid(D = 10 ^ - 0.1, k_ratio = 10 ^ -2)
summary$tau <- - log(summary$D)
summary$k <- as.factor(summary$k_ratio * 1e9)
summary <- run_sims(summary, rep = 1e3, r = 1.023, res = TRUE)
summary


data <- simulate(equilibrate = TRUE, num_mutants = 10, rep = 1, time = 50, tau = 1e4, D = 1, k = 1e9, r = 3, flow = 1)[[1]]
find_W(r = 3, D = 1, c = 1e9, tau = 1e4, k = 1e9, flow = 1)
data <- simulate(equilibrate = TRUE, num_mutants = 10, rep = 1, time = 50, tau = 0.1, D = exp(-0.1), k = 1e9, r = 3 * 1.023, flow = 0)[[1]]
log10(total_pop(data, strains = "W"))
log_plot(data)
find_W(r = 1.5 * (1 + 10), D = exp(-0.1), c = 1e9, tau = 0.1, k = 1e10)

summary <- expand.grid(D = 10 ^ - seq(0, 0.1, by = 0.02), k_ratio = 10 ^ -1)
summary$tau <- - log(summary$D)
summary$k <- as.factor(summary$k_ratio * 1e9)
summary <- run_sims(summary, rep = 3e2, r = r_res, res = TRUE, mu = 0)
summary$theory <- approx1_theory(summary$D, r = 1, w = 0.1) * 1e-9 *
    find_W(r = r_res / r_adj * (1 + summary$k_ratio), D = summary$D, c = 1e9, tau = summary$tau, k = 1e9 * summary$k_ratio)

for(i in seq_len(nrow(summary))) {
    summary$theory[i] <- approx1_theory(summary$D[i], r = - log(summary$D[i]) / summary$tau[i], w = 0.1) * 1e-9 *
        find_W(r = r_res, D = summary$D[i], c = 1e9, tau = summary$tau[i], k = 1e8)
    print(i / nrow(summary))
}

# plot(summary$D, summary$total_pop, log = "xy")

summary$W <- find_W(r = r_res, D = summary$D, c = 1e9, tau = summary$tau, k = 1e7)
data <- simulate(num_mutants = 1, mu = 0, rep = 1, time = 50, tau = 0.3 * log(10), D = 10 ^ - 0.3, k = 1e7, r = r_res, flow = 0)[[1]]
# log10(total_pop(data, strains = "W"))
log_plot(data)
summary

hist(log10(summary$compare))
summary$k <- as.factor(round(summary$tau,2))
base_plot(summary[summary$tau == 24 * 2 ^ -0, ]) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper, color = k),
        linewidth = 0.8) +
    geom_point(aes(x = D, y = theory, color = k), size = 3, shape = 15) +
    geom_point(aes(x = D, y = rate, color = k), size = 3) +
    scale_color_scico_d(palette = "roma")

qw <- summary[summary$D > 0.9, ]
qwe <- qw
qw$old_pop <- find_W(r = r_res / r_adj * (1 + qw$k_ratio), D = qw$D, c = 1e9, tau = qw$tau, k = 1e9 * qw$k_ratio)
qw$new_pop <- find_W(r = r_res * (1 + qw$k_ratio), D = qw$D, c = 1e9, tau = qw$tau, k = 1e9 * qw$k_ratio)
qw$rate2 <- qw$rate * qw$new_pop / qw$old_pop
qwe$rate <- qw$rate2
summary <- rbind(summary, qwe)

find_W(1.53, D, c, k, tau)
RW_ode(0.99 * 1e9, r, D, c, k, tau)
find_W(r, D, c, k, tau)
R0=1e9 * (1 - D)
r=1.5 * (1 + k / c)
D=10 ^ -0.001
c=1e9
k=1e7
tau= - log(D)
print(c(r,D,c,k,tau))

# plot a function comparing R0 and func(R0) on the interval of R0 in (0,c)
df <- data.frame(R0 = 10 ^ seq(1, log10(c - 1), length.out = 1e2))
for (i in seq_len(nrow(df))) {
  df$diff[i] <- func(df$R0[i])
}
plot(df, log = "x", cex.axis = 2, cex.lab = 2)
tail(df)

### testing chemostat model!

summary <- expand.grid(D = 10 ^ - seq(0, 0.1, by = 0.01),
    k_ratio = 10 ^ seq(-2, 1, by = 1))
summary$tau <- - log(summary$D)
summary$k <- as.factor(summary$k_ratio * 1e9)
data <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    tau <- summary$tau[i]
    k_ratio <- res *
        ifelse("k_ratio" %in% names(summary), summary$k_ratio[i], 1)
    data[[i]] <- simulate(
        # seed = i,
        rep = 1e0,
        time = 10,
        dt = 1e-1,
        tau = ifelse(tau == 0, 1e4, tau),
        max_step = Inf,
        D = D,
        flow = ifelse(D == 1, flow, 0),
        R0 = 1e9,
        k = 1e9 * k_ratio,
        alpha = 1,
        r = r * (1 + k_ratio), # * ifelse(D == 1, 1, 1.023),
        # Adaptivetau causes observed growth rate to be lower than expected
        w = 0.1,
        mu = 0,
        N = 1e9,
        num_mutants = 1e0,
        loci = NULL
    )
    summary[i, c("rate", "ci_lower", "ci_upper")] <- func(data[[i]])
    pop <- total_pop(data[[i]][[1]], strains = "W")
    summary[i, "pop"] <- mean(pop)
    # se <- sd(pop) / sqrt(length(pop))
    # summary[i, c("pop", "pop_lower", "pop_upper")] <- mean(pop) + se * qnorm(c(0.5, 0.025, 0.975))
    print(i / nrow(summary))
}
log_plot(data[[12]][[1]][data[[1]][[1]]$rep <= 10, ])
summary$theory <- approx1_theory(summary$D, r = 1, w = 0.1) * 1e-9 *
    find_W(r = r_res / r_adj * (1 + summary$k_ratio), D = summary$D, c = 1e9, tau = summary$tau, k = 1e9 * summary$k_ratio)

base_plot(summary) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper, color = k),
        linewidth = 0.8) +
    # geom_line(aes(x = D, y = theory, color = k), size = 1) +
    geom_point(aes(x = D, y = rate, color = k), size = 3) +
    scale_color_scico_d(palette = "roma")

base_plot(summary) +
    # geom_errorbar(aes(x = D, ymin = pop_lower, ymax = pop_upper, color = k),
    #     linewidth = 0.8) +
    geom_point(aes(x = D, y = pop, color = k), size = 3) +
    scale_color_scico_d(palette = "roma")

kr = 1e0
find_W(r_res * (1 + kr), 10 ^ -0.001, 1e9, 1e9 * kr, 1, 1) -
find_W(r_res / r_adj * (1 + kr), 10 ^ -0.001, 1e9, 1e9 * kr, 1, 1)

D <- 10 ^ - seq(0, 4, 0.1)
W <- find_W(r_res * 2 * (1 + kr), D, 1e9, 1e9 * kr, - log(D), 1)
plot(D, W, log = "xy", cex.axis = 2, cex.lab = 2)

