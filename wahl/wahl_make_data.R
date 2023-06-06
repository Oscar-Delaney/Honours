# Critique of Wahl papers
source("wahl/wahl_code.R")
library(survival)
library(hypergeo)

# calculate the theoretical solution
solve_phi <- Vectorize(
    function(D, s) {
        # Define the function whose root we want to find
        f <- function(phi) {
            return(phi - (D * phi + (1 - D)) ^ (D ^ -(1 + s)))
        }
        # Use uniroot to find the root
        phi <- uniroot(f, lower = 0, upper = 1-1e-10)
        return(phi$root)
    }
)

# theoretical rate at which beneficial mutations fix
fix_rate <- function(D, r, s, N, mu) {
    return(ifelse(D == 1, s / (s + 1), 
    (-N * mu * r * (D ^ -s - 1) ^ 2) /
    (s * log(D) * (D ^ -(s + 1) - 1))))
}

# Hypergeometric function
hyper <- function(s, z) {
    return(Re(hypergeo(1, 1 / (1 + s), 1 + 1 / (1 + s), z)))
}

# theoretical rate function
theory <- function(D, r, s, N, mu) {
    return(-r * N * mu / log(D) * (hyper(s,(1 - D) / (D * (D ^ s - 1))) -
        D * hyper(s, (1 - D) * (D ^ s) /  (D ^ s - 1))))
}


phi <- function(D, s) {
    return(ifelse(D == 1, 1 / (1 + s), (1 - D) / (D ^ -s - D)))
}

theta <- function(D, s) {
    return(D * (D ^ -s - 1) / s)
}

# evaluate a set of simulation results
metric <- function(summary, data, threshold = 1e2) {
    for (i in seq_len(nrow(summary))) {
        # define the minimum population size for a mutant to be ~fixed
        c <- threshold / data[[i]][[2]]$D
        dt <- data[[i]][[2]]$dt # time step
        max_pop <- data[[i]][[2]]$init["W"] / data[[i]][[2]]$D
        mut_rate <- data[[i]][[2]]$m1
        mut_per_gen <- 1 # mut_rate * max_pop
        # find the time at which the mutant first appears en route to fixation
        survival <- data[[i]][[1]] %>%
            group_by(rep) %>%
            summarise(
                fix = any(value[variable == "M"] > c),
                t = max(time[variable == "M" & near(value, 0)], 0) + dt
            )

        # fit the model
        exp_fit <- survreg(Surv(t, fix) ~ 1, data = survival, dist = "exp")

        # extract the log-scale estimate and standard error
        log_estimate <- exp_fit$coefficient[1]
        log_se <- sqrt(exp_fit$var)[1]

        # compute a 95% confidence interval on the log scale
        log_ci <- log_estimate + log_se * qnorm(c(0.025, 0.975))

        # exponentiate to get a confidence interval for the mean on the original scale
        ci <- exp(log_ci)

        # store in summary
        summary$rate[i] <- 1 / exp(log_estimate) / mut_per_gen
        summary$ci_lower[i] <- 1 / ci[2] / mut_per_gen
        summary$ci_upper[i] <- 1 / ci[1] / mut_per_gen
        summary$naive_mean[i] <- 1 / mean(survival$t) / mut_per_gen
        summary$wins[i] <- mean(survival$fix)

        print(i / nrow(summary))
    }
    return(summary)
}

summary <- metric(summary, data)
summary$theory <- theory(summary$D, r, s, summary$N0, m1)
summary
log_plot(data[[1]][[1]][data[[1]][[1]]$rep == 6 & data[[1]][[1]]$time < 100, ])
sols <- data[[4]][[1]]
sols2 <- sols[sols$variable == "W" & sols$time <= 0.24 & sols$rep == 1,]
print(sols2, n=50)
s <- 0.1
N0 <- 1e8
time <- 100
m1 <- 1e-9

# Wahl 1 No resource constraints
r <- 1.024
summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 4.1), N0 = 10 ^ seq(9, 9, 1))
data <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    N0 <- summary$N0[i]
    # m1 <- summary$m1[i]
    data[[i]] <- simulate(
        seed = NULL,
        rep = 1e3,
        time = time,
        dt = 1e-2,
        tau = - log(D),
        D = D,
        N0 = N0,
        k = 0,
        alpha = 0,
        mu = c(W = r, M = r * (1 + s)),
        init_W = round(N0 * D),
        m1 = 1e-9,
        num_mutants = 10
    )
    print(i / nrow(summary))
}


# Wahl 2 Constant resource concentration in dilution media
r <- 1.1
# summary <- expand.grid(tau = seq(0.1, 3, 0.7), D = seq(0.1, 3, 0.7))
summary <- data.frame(D = 10 ^ - seq(0.05, 4, by = 0.2))
data2 <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    data2[[i]] <- simulate(
        seed = NULL,
        rep = 1e3,
        time = time,
        tau = - log(D),
        D = D,
        N0 = N0,
        k = 1e-2 * N0,
        alpha = 1,
        mu = c(W = r, M = r * (1 + s)),
        init = c(W = round(N0 * D), M = 0)
    )
    print(i/nrow(summary))
}

# Wahl 3 Constant total resource supply per time
r <- 1.1
# summary <- expand.grid(tau = seq(0.1, 2, 0.1), D = seq(0.1, 0.95, 0.05))
summary <- data.frame(D = 10 ^ - seq(0.05, 4, by = 0.2))
# summary$D <- exp(-summary$tau)
data3 <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    data3[[i]] <- simulate(
        seed = NULL,
        rep = 1e3,
        time = time,
        tau = - log(D),
        D = D,
        N0 = N0 * log(D) / (D - 1),
        k = 1e-2 * N0,
        alpha = 1,
        mu = c(W = r, M = r * (1 + s)),
        init = c(W = round(N0 * D), M = 0)
    )
    print(i/nrow(summary))
}

# calculate summary statistics
summary$wins <- 0
for (i in seq_len(nrow(summary))) {
    # count the number of runs on which M was at least 1e-5 of W at the end
    end <- data[[i]][[1]][data[[i]][[1]]$time == max(data[[i]][[1]]$time), ]
    M <- end[end$variable == "M",]$value
    summary$theta[i] <- mean(M) * data[[i]][[2]]$D
    # W <- end[end$variable == "W",]$value
    # summary$wins[i] <- mean(ifelse(W==0, M, M / W) > 1e2 / N0)
    # summary$unsure[i] <- 1 - mean(M == 0) - summary$wins[i]
    # rep <- data[[i]][[2]]$rep
    # ci <- binom.test(summary$wins[i] * rep, rep, conf.level = 0.95)$conf.int
    # summary$ci_lower[i] <- ci[1]
    # summary$ci_upper[i] <- ci[2]
}

summary$theory <- N0 * m1 * r * summary$D * (1 - phi(summary$D, s))
summary$theory_old <- fix_rate(summary$D, r, s, N0, m1)
summary$theory_full <- theory(summary$D, r, s, N0, m1)
summary[c("D", "theory", "theory2")]
png("images/test.png", width = 12, height = 10, units = "in", res = 300)
ggplot(summary, aes(x = D, y = rate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)) +
    geom_line(aes(y = theory, color = "new")) +
    geom_line(aes(y = theory_old, color = "old")) +
    geom_line(aes(y = theory_full, color = "full")) +
    scale_color_manual(values = c("new" = "red", "old" = "blue", "full" = "green")) +
    theme_light() +
    scale_x_log10() +
    scale_y_log10() +
    # scale_y_continuous(limits = c(0,1)) +
    labs(
        # title = "Optimal Dilution Ratio (resource unconstrained)",
        x = "mutation rate",
        y = "fixation rate",
        color = "theoretical solution"
    ) +
    theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )
dev.off()

ggplot(summary, aes(x = D, y = wins)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)) +
    geom_line(aes(y = 1 - phi, color = "theory")) +
    theme_light() +
    scale_x_log10() +
    # scale_y_log10() +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values = c("theory" = "red")) +
    labs(
        title = "Optimal Dilution Ratio (resource unconstrained)",
        x = "D",
        y = "1 - phi: probability of a single mutant at t = 0 eventually fixing",
        color = NULL
    ) +
    theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )
# dev.off()

data5 <- list()
r <- 1.2
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    data5[[i]] <- simulate(
        seed = NULL,
        rep = 1e0,
        dt = 0.01,
        time = time,
        tau = - log(D),
        D = D,
        influx = c(A1 = 0, A2 = 0),
        N0 = N0,
        k = 1e-1 * N0,
        alpha = 1,
        m1 = 1e-9,
        m2 = 0,
        mu = c(W = r, M = r * (1 + s)),
        init = c(W = round(N0 * D), M = 0, R2 = 0, R12 = 0)
    )
    print(i/nrow(summary))
}

summary$N_plus <- 0
summary$N_minus <- 0
for (i in seq_len(nrow(summary))) {
    # Find the time at which the final bottleneck occurred
    sols <- data5[[i]][[1]]
    sols <- sols[sols$variable == "N",]
    last_bneck <- sols %>%
        filter(rep == 1) %>%
        filter(value > lag(value)) %>%
        slice_tail(n = 1) %>%
        pull(time)
    dt <- 0.1
    summary$N_plus[i] <- mean(sols[near(sols$time, last_bneck), ]$value)
    summary$N_minus[i] <- mean(sols[near(sols$time, last_bneck - dt), ]$value)
}

png("wahl_nutrients2_v2.png", width = 12, height = 10, units = "in", res = 300)
ggplot(summary) +
    geom_point(aes(x = D, y = N_minus, color = "Before"), size = 3) +
    geom_point(aes(x = D, y = N_plus, color = "After"), size = 3) +
    theme_light() +
    scale_x_log10() +
    # scale_y_log10() +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      breaks = 10^seq(0, 20, by = 1),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = c("Before" = "blue", "After" = "red"), 
                       name = "Nutrient pulse: ") +
    # scale_y_continuous(limits = c(0,1)) +
    labs(
        title = "Nutrient use (constant media resource conc)",
        # title = "Nutrient use (constant resource flux)",
        x = "D",
        y = "Nutrient concentration"
    ) +
    theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )
dev.off()

filtered <- data5[[4]][[1]] %>%
    filter(variable %in% c("N","W")) %>%
    filter(time > 0)

ggplot() +
    geom_point(data = filtered, aes(x = time, y = value, color = variable)) +
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

ggplot(summary, aes(x = D, y = tau)) +
    geom_tile(aes(fill = wins)) +
    scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
    labs(x = "D", y = "tau", fill = "proportion",
            title = "Optimal Dilution Ratio (constant resource flux)",
            subtitle = "proportion of runs where M becomes established") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )

# unused
survival$rexp <- sort(rexp(100, 1 / summary$time[1]))
plot(survival$rexp, sort(survival$t), log = "xy", xlim = c(1e-2, 5e1), ylim = c(1e-2, 5e1))
abline(0, 1)

# # Save the simulation results to a file
save(data, file = "C:/Users/s4528540/Downloads/results/data_05_06.RData")
# read the data file we just wrote back as data
load("C:/Users/s4528540/Downloads/results/data4.RData")
