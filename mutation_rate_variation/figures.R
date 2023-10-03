source("stochastic.R")

# Time to reach target
target_time <- function(sol, target = 1, strains = c("N_S", "N_A", "N_B", "N_AB")) {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep, time) %>%
    summarise(total = sum(value), .groups = "drop") %>%
    group_by(rep) %>%
    summarise(t = time[min(which(diff(sign(total - target)) != 0), Inf)]) %>%
    pull(t)
}

# Whether target was hit
target_hit <- function(sol, target = 1, strains = c("N_S", "N_A", "N_B", "N_AB")) {
  !is.na(target_time(sol, target, strains))
}

theory_d <- function(C, zeta = 1) {
    1 / (1 + zeta / C)
}

theory_N <- function(delta, d_A, d_B, mu, init_S, bcidal, bstatic) {
    mu <- mu * (1 - bstatic * d_A) * (1 - bstatic * d_B)
    N <- init_S * mu / (delta + bcidal * (d_A + d_B) - mu)
    ifelse(N < 0, Inf, N)
}

theory_phi <- function(m_A, c_A, c_B, bcidal, bstatic, delta, zeta, mu) {
    d_A <- theory_d(c_A, 1)
    d_A_R <- theory_d(c_A, zeta)
    d_B <- theory_d(c_B, 1)
    d_B_R <- theory_d(c_B, zeta)
    pmin(1, m_A * (delta + bcidal * (d_B + d_A_R)) / (mu * (1 - bstatic * d_A_R) * (1 - bstatic * d_B)) +
    (1 - m_A) * (delta + bcidal * (d_A + d_B_R)) / (mu * (1 - bstatic * d_A) * (1 - bstatic * d_B_R)))
}

theory_extinct <- function(phi, N, m) {
    (1 - m * (1 - phi)) ^ N
}

run_sims <- function(summary, config_only = FALSE, rep = 1e2, zeta = 1e9,
c = 5, kappa = 1) {
    for (i in seq_len(nrow(summary))) {
        m_A <- 1 / (1 + 1 / summary$m_ratio[i])
        data <- simulate(
            init = c(N_S = init_S, N_A = 0, N_B = 0, N_AB = 0),
            R0 = 1e10,
            k = 0,
            alpha = 0,
            supply = 1e8,
            mu = 1,
            bcidal_A = summary$cidal_A[i],
            bcidal_B = summary$cidal_B[i],
            bstatic_A = 1 - summary$cidal_A[i],
            bstatic_B = 1 - summary$cidal_B[i],
            zeta_A = c(N_S = 1, N_A = zeta, N_B = 1, N_AB = zeta),
            zeta_B = c(N_S = 1, N_A = 1, N_B = zeta, N_AB = zeta),
            delta = 1 - 1 / (1 + c ^ -1),
            time = 50,
            tau = 1e4,
            max_step = 1e-1,
            kappa_A = kappa, kappa_B = kappa,
            rep = rep,
            HGT = 0,
            dose_rep = 1,
            dose_gap = 1e4,
            influx = c * c(C_A = summary$cA[i], C_B = 1 - summary$cA[i]),
            cycl = FALSE,
            m_A = m * m_A, m_B = m * (1 - m_A),
            d_A = 0, d_B = 0,
            deterministic = FALSE,
            config_only = config_only
        )
        if (config_only) {
            v <- with(data, as.list(rates(c(N_S = 1, N_A = 1, N_B = 1, N_AB = 1, R = R0, influx * pattern), data, 0)))
            summary$phi[i] <- with(v, m_A * pmin(1, N_A_death / N_A_growth) +
                (1 - m_A) * pmin(1, N_B_death / N_B_growth))
            summary$N[i] <- with(v, init_S / (pmax(1e-30, S_death / S_growth - 1)))
            # summary$extinct[i] <- (1 - m * (1 - summary$phi[i])) ^ summary$N[i]
            summary$extinct[i] <- exp(-m * (1 - summary$phi[i]) * summary$N[i])
        } else {
        wins <- !target_hit(data[[1]], target = 1e3, strains = c("N_A", "N_B"))
        summary$extinct[i] <- mean(wins)
        ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
        summary$ymin[i] <- ci[1]
        summary$ymax[i] <- ci[2]
        print(i / nrow(summary))
        }
    }
    return(summary)
}
# plot with two independent variables and one dependent variable
summary_plot <- function(summary, var) {
    # Find the y-values that maximize 'var' for each x-value
    summary$approx <- 1 / (1 + (summary$m_ratio) ^ 0.5)
    # Rename the factor levels with the desired prefixes
    summary$cidal_A <- ifelse(summary$cidal_A == 0, "Drug A: Bacteriostatic", "Drug A: Bactericidal")
    summary$cidal_B <- ifelse(summary$cidal_B == 0, "Drug B: Bacteriostatic", "Drug B: Bactericidal")
    labels <- expand.grid(m_ratio = 500, cA = 0.95, cidal_A = unique(summary$cidal_A), cidal_B = unique(summary$cidal_B))
    labels$label <- LETTERS[nrow(labels):1]

    max_df <- summary %>%
        group_by(cidal_A, cidal_B, m_ratio) %>%
        arrange(desc(get(var))) %>%
        slice_head(n = 1) %>%
        ungroup()

    ggplot(summary, aes(x = log2(m_ratio), y = cA)) +
        geom_tile(aes(fill = (get(var)))) +
        geom_line(aes(x = log2(m_ratio), y = approx), color = "green", linewidth = 3) +
        geom_point(data = max_df, aes(y = cA, x = log2(m_ratio)), color = "yellow", size = 3) +
        facet_grid(cidal_B ~ cidal_A) +
        geom_text(data = labels, aes(label = label), size = 15, fontface = "bold") +
        labs(fill = "P(extinct)") +
        theme_minimal() +
        labs(
            x = "Log2 odds of A-resistance:B-resistance",
            y = "Proportion of dose that is drug A"
        ) +
        scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
        theme(
            axis.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size = 25),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.position = "bottom",
            legend.key.width = unit(2, "cm"),
            legend.spacing.x = unit(1, "cm"),
            strip.text = element_text(size = 25, face = "bold")
        )
}

# Create a grid of parameter combinations
init_S <- 1e9
m <- 1e-9

### Basic
summary <- expand.grid(cA = seq(0, 1, 0.03), m_ratio = 2 ^ seq(-10, 10, 0.5), cidal_A = c(0, 1), cidal_B = c(0, 1))
basic <- run_sims(summary, config_only = TRUE)

pdf("mutation_rate_variation/basic.pdf", width = 20, height = 20)
summary_plot(basic, "extinct")
dev.off()

### Incomplete resistance
in_res <- run_sims(summary, zeta = 25, c = 5, config_only = TRUE)

pdf("mutation_rate_variation/in_res.pdf", width = 20, height = 20)
summary_plot(in_res, "extinct")
dev.off()

### Kappa variation
kappa3 <- run_sims(summary, kappa = 3, c = 2, config_only = TRUE)

# pdf("mutation_rate_variation/kappa3.pdf", width = 20, height = 20)
summary_plot(kappa3, "extinct")
# dev.off()

kappa0.5 <- run_sims(summary, kappa = 0.1, c = 0.5, config_only = TRUE)
summary_plot(kappa0.5, "extinct")

### Simulations
summary <- expand.grid(cA = seq(0.01, 0.99, length.out = 10), m_ratio = 2 ^ seq(-10, 10, length.out = 10), cidal_A = c(1), cidal_B = c(1))
simulations <- run_sims(summary, config_only = FALSE, rep = 1e2)
summary_plot(simulations, "extinct")
