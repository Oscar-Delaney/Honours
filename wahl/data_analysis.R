# Generate and analyse data
run_sims <- function(summary, rep = 1, time = 50, w = 0.1, r = 1, mu = 1e-9,
    res = TRUE, flow = 1, num_mutants = 1e2, loci = NULL) {
    func <- ifelse(is.null(loci), metric, metric_ci)
    data <- list()
    for (i in seq_len(nrow(summary))) {
        D <- summary$D[i]
        tau <- summary$tau[i]
        k_ratio <- res *
            ifelse("k_ratio" %in% names(summary), summary$k_ratio[i], 1)
        data <- simulate(
            # seed = i,
            rep = rep,
            time = time,
            dt = 1e-2,
            tau = ifelse(tau == 0, 1e4, tau),
            max_step = Inf,
            D = D,
            flow = ifelse(D == 1, flow, 0),
            media = 1e9,
            k = 1e9 * k_ratio,
            alpha = 1 * res,
            r = r * (1 + k_ratio), # * ifelse(D == 1, 1, 1.023),
            # Adaptivetau causes observed growth rate to be lower than expected
            w = w,
            mu = mu,
            N = 1e9,
            num_mutants = num_mutants,
            loci = loci
        )
        summary[i, c("rate", "ci_lower", "ci_upper")] <- func(data)
        # summary[i, "total_pop"] <- log10(mean(total_pop(data[[1]], strains = "W")))
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
    return(r * w / (1 + w) * ifelse(D == 1, 1, log(D^-1) / (D^-1 -1)))
}

# Even more approximate rate function
approx2_theory <- function(D, r, w) {
    return(r * w / (1 + w) * sqrt(D))
}

# theoretical resource constrained adaptation rate
theory_res <- function(D, tau, r, w = 0.1, media = 1e9, k = 1e9, mu = 1e-9, flow = 1) {
    k <- rep(k, length(D))[1:length(D)]
    N <- ifelse(k == 0, media, find_W(r, D, media, k, tau, flow))
    x <- ifelse(D == 1 | tau == 0, flow, log(D) ^ 2 / (1 / D - 1) / tau)
    return(w / (1 + w) * mu * N * x) #  / (1 + w)
}

# rate at a given time within [0, tau)
rate_at_t <- function(D, r, w, t) {
  (D * (-1 + D^w) * exp(r * t) * r) /
    (-1 + D^(1 + w) * exp(r * (1 + w) * t) - D^w * (-1 + exp(r * (1 + w) * t)))
}

# Count the number of mutants likely en route to fixation
fixed <- function(data) {
    # find the time just after the last bottleneck
    endpoint <- min(data[[2]]$time, data[[1]] %>%
        filter(variable == "W", rep == 1) %>%
        filter(value - lag(value) < 0) %>%
        tail(1) %>%
        pull(time))
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
