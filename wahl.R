# Critique of Wahl papers
source("stochastic.R")

log_plot(data3[[11]][[1]][data1[[1]][[1]]$rep >= 990, ], type = "all")

s <- 0.1
N0 <- 3e8
time <- 100

# Wahl 1 No resource constraints
r <- 1.03
summary <- data.frame(D = 10 ^ - seq(0.05, 4, by = 0.2))
data1 <- list()
for (i in seq_len(nrow(summary))) {
    D <- summary$D[i]
    data1[[i]] <- simulate(
        seed = NULL,
        rep = 1e3,
        time = time,
        tau = - log(D),
        D = D,
        influx = c(A1 = 0, A2 = 0),
        N0 = N0,
        k = 1e-2 * N0,
        alpha = 0,
        m1 = 1e-9,
        m2 = 0,
        mu = c(S = r, R1 = r * (1 + s)),
        init = c(S = round(N0 * D), R1 = 0, R2 = 0, R12 = 0)
    )
    print(i/nrow(summary))
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
        influx = c(A1 = 0, A2 = 0),
        N0 = N0,
        k = 1e-2 * N0,
        alpha = 1,
        m1 = 1e-9,
        m2 = 0,
        mu = c(S = r, R1 = r * (1 + s)),
        init = c(S = round(N0 * D), R1 = 0, R2 = 0, R12 = 0)
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
        influx = c(A1 = 0, A2 = 0),
        N0 = N0 * log(D) / (D - 1),
        k = 1e-2 * N0,
        alpha = 1,
        m1 = 1e-9,
        m2 = 0,
        mu = c(S = r, R1 = r * (1 + s)),
        init = c(S = round(N0 * D), R1 = 0, R2 = 0, R12 = 0)
    )
    print(i/nrow(summary))
}

# calculate summary statistics
summary$wins <- 0
for (i in seq_len(nrow(summary))) {
    # count the number of runs on which R1 was at least 1e-5 of S at the end
    end <- data3[[i]][[1]][data3[[i]][[1]]$time == time, ]
    R1 <- end[end$variable == "R1",]$value
    S <- end[end$variable == "S",]$value
    summary$wins[i] <- mean(ifelse(S==0, R1, R1 / S) > 1e3 / N0)
    rep <- data3[[i]][[2]]$rep
    ci <- binom.test(summary$wins[i] * rep, rep, conf.level = 0.95)$conf.int
    summary$ymin[i] <- ci[1]
    summary$ymax[i] <- ci[2]
}

png("wahl3.png", width = 12, height = 10, units = "in", res = 300)
ggplot(summary, aes(x = D, y = wins)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    theme_light() +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,1)) +
    labs(
        title = "Optimal Dilution Ratio (constant resource flux)",
        x = "D",
        y = "Probability that a mutant fixes in 100 hours"
    ) +
    theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25)
    )
dev.off()

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
        mu = c(S = r, R1 = r * (1 + s)),
        init = c(S = round(N0 * D), R1 = 0, R2 = 0, R12 = 0)
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
    filter(variable %in% c("N","S")) %>%
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
            subtitle = "proportion of runs where R1 becomes established") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )
