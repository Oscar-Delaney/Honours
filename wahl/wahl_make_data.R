# Critique of Wahl papers
source("wahl/wahl_code.R")
library(hypergeo)
library(scales)

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
approx_theory <- function(D, r, s) {
    return(r * s * log(D^-1) / (D^-1 -1))
}

# Probability a new mutant at the beginning of a growth phase will go extinct
phi <- function(D, s) {
    return(ifelse(D == 1, 1 / (1 + s), (1 - D) / (D ^ -s - D)))
}

# Define the differential equation
dw_dt <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dW <- mu * W / (1 + k / (influx - W))
    return(list(dW))
  })
}

# Define the function
w_at_tau <- Vectorize(function(w0, influx, tau, D, mu, k) {
  parameters <- c(influx = influx, mu = mu, k = k)
  out <- ode(y = c(W = w0), times = c(0, tau), func = dw_dt, parms = parameters)
  return(out[2, "W"])  # Returns W at tau
})

delta_w <- function(w0, influx, tau, D, mu, k) {
  return(w_at_tau(w0, influx, tau, D, mu, k) * D - w0)
}

# Find root
division_rate <- Vectorize(function(influx, tau, D, mu, k) {
  result <- tryCatch({
    root <- uniroot(delta_w, interval = c(1, influx), influx = influx, tau = tau, D = D, mu = mu, k = k)
    return(root$root * (1 / D - 1) / tau)
  }, 
  error = function(e) {
    return(NA)
  })
  return(result)
})

mu <- 1
influx <- 1e9
k <- 1e9
plot(summary$D, summary$rate, log = "x", type = "l", xlab = "D", ylab = "rate")



# Create a new variable for D^2
summary$D_sq <- summary$D^2

# Fit the quadratic model
fit <- lm(rate ~ D + D_sq, data = summary)

# Print the summary of the fit
summary(fit)


# evaluate a set of simulation results
metric <- function(data) {
    # extract the mutation rate and population size
    m1 <- data[[2]]$m1
    N0 <- data[[2]]$N0
    # find the time just after the last bottleneck
    endpoint <- data[[1]] %>%
        filter(variable == "W", rep == 1) %>%
        filter(value - lag(value) < 0) %>%
        tail(1) %>%
        pull(time)
    # find the likelihood of a new mutation at t=0 going extinct 
    current_phi <- phi(data[[2]]$D, s)
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
    fixation_rate <- fixed$n_hat / endpoint / (N0 * m1)
    se <- sd(fixation_rate) / sqrt(length(fixation_rate))
    ci <- mean(fixation_rate) + se * qnorm(c(0.5, 0.025, 0.975))
    return(ci)
}

# find the time at which the wild-type population is less than half its initial value
metric2 <- function(data) {
    init_W <- data[[2]]$init_W
    m1 <- data[[2]]$m1
    N0 <- data[[2]]$N0
    endpoint <- data[[1]] %>%
        group_by(rep) %>%
        filter(variable == "W") %>%
        filter(value < init_W / 2) %>%
        summarize(time = first(time), .groups = 'drop') %>%
        pull(time)
    fixation_rate <- 1 / endpoint / (N0 * m1)
    se <- sd(fixation_rate) / sqrt(length(fixation_rate))
    ci <- mean(fixation_rate) + se * qnorm(c(0.5, 0.025, 0.975))
    return(ci)
}

metric2(data[[1]])
summary <- metric(summary, data2)
summary$theory <- theory(summary$D, r, s)
summary$approx <- approx_theory(summary$D, r, s)

log_plot(data[[1]][[1]][data[[1]][[1]]$rep == 3 & data[[1]][[1]]$time > 0, ])

data[[1]][[1]] %>%
    filter(rep == 1, variable == "W") %>%
    summarise(max(value)) %>%
    pull() / data[[1]][[2]]$N0


# Wahl 1 No resource constraints
s <- 0.1
time <- 50
r <- 1.023
N0 <- 1e9
summary <- expand.grid(D = 10 ^ - seq(0.1, 2, by = 0.1))
summary$m1 <- summary$D ^ - 0.5 * 1e-9
data <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    m1 <- summary$m1[i]
    data[[i]] <- simulate(
        seed = NULL,
        rep = 1e2,
        time = time,
        dt = 1e-1,
        tau = - log(D),
        D = D,
        N0 = N0,
        k = 0,
        alpha = 0,
        r = r,
        s = s,
        init_W = round(N0 * D),
        num_mutants = 1e2
    )
    summary[i, c("rate", "ci_lower", "ci_upper")] <- metric(data[[i]])
    print(i / nrow(summary))
}

# Wahl 2 Constant resource concentration in dilution media
s <- 0.1
time <- 200
N0 <- 1e9
k_ratio <- 1e-2
r <- 1.1 * (k_ratio + 1)
summary <- expand.grid(D = 10 ^ - seq(0.1, 0.2, by = 0.1))
summary$m1 <- summary$D ^ - 0.5 * 1e-9
data <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    tau <- - log(D)
    m1 <- 1e-9
    k <- k_ratio * N0
    if (D >= exp(-r * tau)) {
        data[[1]] <- simulate(
            seed = NULL,
            rep = 1e1,
            time = time,
            dt = 1e-1,
            tau = tau,
            D = D,
            N0 = N0,
            m1 = m1,
            k = k,
            alpha = 1,
            r = r,
            s = s,
            init_W = round(N0 * D),
            num_mutants = 1e2
        )
        summary[i, c("rate", "ci_lower", "ci_upper")] <- metric(data[[1]])
    }
    print(i / nrow(summary))
}

# Wahl 3 Constant total resource supply per time
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

summary <- metric(summary, data2)
summary <- data.frame(D = 10 ^ - seq(0.01, 5, by = 0.01), s = 0.01, r = 0.1)
summary$rate <- theory(summary$D, summary$r, summary$s)
summary$approx <- summary$r * summary$s * log(1 / summary$D) / (1 / summary$D - 1)


# png("images/test.png", width = 12, height = 10, units = "in", res = 300)
ggplot(summary, aes(x = D, y = rate)) +
    geom_point(size = 3) +
    # geom_line() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)) +
    geom_line(aes(y = theory, color = "full")) +
    geom_line(aes(y = approx, color = "approx")) +
    scale_color_manual(values = c("full" = "red", "approx" = "blue")) +
    theme_light() +
    scale_x_log10() +
    scale_y_log10() +
    # scale_y_continuous(limits = c(0,1)) +
    labs(
        # title = "Optimal Dilution Ratio (resource unconstrained)",
        x = "fitness benefit of mutation",
        y = "fixation rate (loci per hour per (mutations per generation))",
        color = "theory"
    ) +
    theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )


summary <- expand.grid(D = 10 ^ - seq(0.01, 0.2, by = 0.01), tau = 2 ^ - seq(0, 3, by = 0.1))
summary$rate <- division_rate(influx, summary$tau, summary$D, mu, k)

# png("images/test.png", width = 12, height = 10, units = "in", res = 300)
ggplot(summary, aes(x = D, y = log2(tau))) +
    geom_tile(aes(fill = log10(rate))) +
    scale_x_log10() +
    # scale the y axis on a log2 scale
    # scale_y_continuous(trans = "log2") +
    scale_fill_gradient(low = "white", high = "blue",
        breaks = pretty_breaks(n=2)) +
    labs(x = "D", y = "log2 tau (hours)", fill = "log10 fixation rate",
            title = "Optimal Dilution Ratio (constant resource flux)") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
    )
# dev.off()

# # Save the simulation results to a file
save(data, file = "C:/Users/s4528540/Downloads/results/data_05_06.RData")
# read the data file we just wrote back as data
load("C:/Users/s4528540/Downloads/results/data4.RData")
