# Generate and analyse data
run_sims <- function(summary, rep = 1, time = 50, w = 0.1, r = 1, mu = 1e-9, dt = 1e-2,
    res = TRUE, flow = 1, num_mutants = 1e2, loci = NULL, summarize = TRUE) {
    for (i in seq_len(nrow(summary))) {
        D <- summary$D[i]
        tau <- summary$tau[i]
        k_ratio <- res *
            ifelse("k_ratio" %in% names(summary), summary$k_ratio[i], 1)
        data <- simulate(
            seed = i,
            rep = rep,
            time = time,
            dt = dt,
            tau = ifelse(tau == 0, 1e4, tau),
            max_step = Inf,
            D = D,
            flow = ifelse(D == 1, flow, 0),
            media = 1e9,
            k = 1e9 * k_ratio,
            alpha = 1 * res,
            r = r * (1 + k_ratio),
            w = w,
            mu = mu,
            N = 1e9,
            num_mutants = num_mutants,
            loci = loci,
            summarize = summarize
        )
        if (summarize) {
            fixed <- data[[1]] %>%
                group_by(rep) %>%
                summarise(n = sum(final_value > 1e1), n_hat = sum(p_fix))
            # estimate the fixation rate and store this
            fixation_rate <- fixed$n_hat / data[[2]]$endpoint
            se <- sd(fixation_rate) / sqrt(length(fixation_rate))
            ci <- mean(fixation_rate) + se * qnorm(c(0.5, 0.025, 0.975))
        } else {
            func <- ifelse(is.null(loci), metric, metric_ci)
            ci <- func(data)
        }
        summary[i, c("rate", "ci_lower", "ci_upper")] <- ci
        print(i / nrow(summary))
    }
    return(summary)
}


# Probability a new mutant at the beginning of a growth phase will go extinct
phi <- function(D, s) {
    return(ifelse(D == 1, 1 / (1 + s), (1 - D) / (D ^ -s - D)))
}

# theoretical resource constrained adaptation rate
theory <- function(D, tau, r, w = 0.1, media = 1e9, k = 1e9, mu = 1e-9, flow = 1, alpha = 1) {
    k <- rep(k, length(D))[seq_along(length(D))]
    N <- ifelse(k == 0 | alpha == 0, media, find_W(r, D, media, k, tau, flow, alpha))
    x <- ifelse(D == 1 | tau == 0, flow, log(D) ^ 2 / (1 / D - 1) / tau)
    return(w / (1 + w) * mu * N * x)
}

# rate at a given time within [0, tau)
rate_at_t <- function(D, r, w, t) {
    u <- D * exp(r * t)
    theta <- (1 / D - 1) / (1 - D ^ w)
    return(r^2 / -log(D) * (u / (1 + theta * u ^ (1 + w))))
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
