source("stochastic.R")

# Total population
total_pop <- function(sol, dt = 1, strains = c("N_S", "N_A", "N_B", "N_AB")) {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep) %>%
    summarise(total = sum(value) * dt) %>%
    pull(total)
}

# Final population
final_pop <- function(sol, strains = c("N_S", "N_A", "N_B", "N_AB")) {
  sol %>%
    filter(time == max(time) & variable %in% strains) %>%
    group_by(rep) %>%
    summarise(final = sum(value)) %>%
    pull(final)
}

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

summary$phi <- theory_phi(summary$m_A, C_tot * summary$cA, C_tot * (1 - summary$cA), bcidal, bstatic, delta, zeta, mu)
summary$N <- theory_N(delta, theory_d(C_tot * summary$cA), theory_d(C_tot * (1 - summary$cA)), mu, init_S, bcidal, bstatic)
summary$extinct <- theory_extinct(summary$phi, summary$N, m)


# plot with two independent variables and one dependent variable
summary_plot <- function(summary, var) {
  # Find the y-values that maximize 'var' for each x-value
  summary$approx <- 1 / (1 + (1 / summary$m_A - 1) ^ -0.5)
  max_df <- summary %>%
    group_by(log2(m_A)) %>%
    arrange(desc(get(var))) %>%
    slice_head(n = 1) %>%
    ungroup()

    ggplot(summary, aes(x = log2(m_A), y = cA)) +
    geom_tile(aes(fill = (get(var)))) +
    geom_line(aes(x = log2(m_A), y = approx), color = "green", linewidth = 2) +
    geom_point(data = max_df, aes(y = cA, x = log2(m_A)), color = "yellow", size = 5) +
    labs(fill = "P(extinct)") +
    theme_minimal() +
    labs(
        x = "Log2 proportion of de novo mutations that are A",
        y = "Proportion of dose that is drug A",
        size = "Probability that MDR emerges"
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
summary_plot(summary, "N")

# Create a grid of parameter combinations
# summary <- expand.grid(dose_rep = seq(1, 12, 1), kappa = seq(0.5, 3, 0.5))
init_S <- 1e9
m <- 1e-9
zeta <- 1e9


summary <- expand.grid(cA = seq(0.5, 1, 0.01), m_A = c(2 ^ -seq(1, 20, 0.2)))

c <- 3
th <- 1.001 + c^-5
# data <- list()
for (i in seq_len(nrow(summary))) {
    data[[i]] <- simulate(
        init = c(N_S = init_S, N_A = 0, N_B = 0, N_AB = 0),
        R0 = 1e10,
        k = 0,
        alpha = 0,
        supply = 1e8,
        mu = 1,
        bcidal_A = th,
        bcidal_B = th,
        bstatic_A = 0,
        bstatic_B = 0,
        zeta_A = c(N_S = 1, N_A = zeta, N_B = 1, N_AB = 1e9),
        zeta_B = c(N_S = 1, N_A = 1, N_B = zeta, N_AB = 1e9),
        delta = 0, # 1 - th / (1 + 1 / c),
        time = 50,
        tau = 1e4,
        max_step = 1e-1,
        kappa_A = 5, kappa_B = 5,
        rep = 1e1,
        HGT = 0,
        dose_rep = 1,
        dose_gap = 1e4,
        influx = c * c(C_A = summary$cA[i], C_B = 1 - summary$cA[i]),
        cycl = FALSE,
        m_A = m * summary$m_A[i], m_B = m * (1 - summary$m_A[i]),
        d_A = 0, d_B = 0,
        deterministic = FALSE,
        config_only = TRUE
    )
    v <- with(data[[i]], as.list(rates(c(N_S = 1, N_A = 1, N_B = 1, N_AB = 1, R = R0, influx * pattern), data[[i]], 0)))
    summary$phi[i] <- with(v, summary$m_A[i] * pmin(1, N_A_death / N_A_growth) +
        (1 - summary$m_A[i]) * pmin(1, N_B_death / N_B_growth))
    summary$N[i] <- with(v, init_S / (pmax(1e-30, S_death / S_growth - 1)))
    summary$extinct[i] <- (1 - m * (1 - summary$phi[i])) ^ summary$N[i]
    summary$extinct2[i] <- exp(-m * (1 - summary$phi[i]) * summary$N[i])
    # wins <- !target_hit(data[[i]][[1]], target = 1e3, strains = c("N_A", "N_B"))
    # summary$wins[i] <- mean(wins)
    # ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
    # summary$ymin[i] <- ci[1]
    # summary$ymax[i] <- ci[2]
    # print(i / nrow(summary))
}
summary_plot(summary, "extinct2")
summary[summary$m_A == 2 ^ - 1.01,]
hist(log10(summary$extinct2 / summary$extinct))
tail(summary)
summary_plot(summary, "wins")
for (i in seq_len(nrow(summary))) {
    wins <- !target_hit(data[[i]][[1]], target = 1e3, strains = c("N_A", "N_B"))
    summary$wins[i] <- mean(wins)
    ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
    # print(mean(!target_hit(data[[i]][[1]], target = 1e1, strains = c("N_A", "N_B"))))
}
final <- final_pop(data[[1]][[1]], c("N_B"))
sort(final[final != 0])
which(final > 0 & final < 1e2)
summary$total_pop / summary$N2
log_plot(data[[1]][[1]][data[[1]][[1]]$rep <= 10, ])
# plot with one independent variable and one dependent variable
summary$m_A_discrete <- as.factor(round(summary$m_A, 3))
ggplot(summary, aes(x = cA, y = wins, color = m_A_discrete)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    theme_light() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
        x = "proportion of dose that is A",
        y = "P(extinct))",
        color = "m_A"
    ) +
    theme(
        plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.position = "bottom"
    )


summary
sol <- data[[1]][[1]][data[[1]][[1]]$rep == 43, ]
sol2 <- sol[sol$variable =="N_S" & sol$time <= 15,]
plot(sol2$time, sol2$value, type = "l")
single_run(config, 1)[, c("time", "N_S", "N", "A1", "A2")]
sol
# Calculate the summary statistics
for (i in seq_len(nrow(summary))) {
    summary$wins[i] <- mean(final_pop(data[[i]][[1]], c("N_A", "N_B")) > 1e2)
    rep <- data[[i]][[2]]$rep
    ci <- binom.test(summary$wins[i] * rep, rep, conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}

for (i in seq_len(nrow(summary))) {
    pops <- total_pop(data[[i]][[1]], strains = "N_S")
    summary$pop[i] <- mean(pops)
    ci <- t.test(pops, conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}

for (i in seq_len(nrow(summary))) {
    wins <- target_hit(data[[i]][[1]], target = 1e2, strains = "N_A")
    summary$wins[i] <- mean(wins)
    ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}




# Next I will make a plot similar to that above, but with a colour scale
ggplot(summary, aes(x = D, y = HGT, color = wins)) +
    geom_point(size = 10) +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal() +
    labs(
        title = "Impact of D and HGT on MDR emergence",
        x = "Dilution fraction D",
        y = "Recombination rate",
        color = "Probability that MDR emerges"
    ) +
    theme(
        legend.position = "bottom",
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.size = unit(2, "cm")
    )

# # Save the simulation results to a file
save(summary, file = "C:/Users/s4528540/Downloads/results/17_07_synergy_2on1.RData")
# read the data file we just wrote back as data
load("C:/Users/s4528540/Downloads/results/17_07_synergy_reciprocal.RData")