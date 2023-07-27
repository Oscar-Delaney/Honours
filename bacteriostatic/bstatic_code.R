source("stochastic.R")
library(gridExtra)


# Total population
total_pop <- function(sol, dt = 0.1, strains = "S") {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep) %>%
    summarise(total = sum(value) * dt) %>%
    pull(total)
}

# Final population
final_pop <- function(sol, strains = c("S", "R1", "R2", "R12")) {
  sol %>%
    filter(time == max(time) & variable %in% strains) %>%
    group_by(rep) %>%
    summarise(final = sum(value)) %>%
    pull(final)
}

# Time to reach target
target_time <- function(sol, target = 1, strains = c("S", "R1", "R2", "R12")) {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep, time) %>%
    summarise(total = sum(value), .groups = "drop") %>%
    group_by(rep) %>%
    summarise(t = time[min(which(diff(sign(total - target)) != 0), Inf)]) %>%
    pull(t)
}

# Whether target was hit
target_hit <- function(sol, target = 1, strains = c("S", "R1", "R2", "R12")) {
  !is.na(target_time(sol, target, strains))
}

run_sims <- function(summary, zeta1 = c(S = 1, R1 = 28, R2 = 1, R12 = 28),
zeta2 = c(S = 1, R1 = 1, R2 = 28, R12 = 28), delta = 0.3, rep = 1, dose_gap = 10,
influx = 3 * c(A1 = 1, A2 = 1), m1 = 1e-9, m2 = 1e-9) {
    summary$bstatic1 <- 1 - summary$bcidal1
    summary$bstatic2 <- 1 - summary$bcidal2
    for (i in seq_len(nrow(summary))) {
        cycl <- ifelse(summary$therapy[i] == "Cycling", TRUE, FALSE)
        d <- ifelse(cycl, 0.1, 0.3)
        data <- simulate(
            init = c(S = ifelse(cycl, 3e8, 1e9), R1 = 0, R2 = 0, R12 = 0),
            N0 = ifelse(summary$resources[i] == "Abundant", 1e11, 1e8),
            k = 1e8,
            alpha = 1,
            supply = 2e8,
            mu = 1,
            bcidal1 = summary$bcidal1[i],
            bcidal2 = summary$bcidal2[i],
            bstatic1 = summary$bstatic1[i],
            bstatic2 = summary$bstatic2[i],
            zeta1 = zeta1,
            zeta2 = zeta2,
            delta = delta,
            time = 60,
            tau = 1e4,
            rep = rep,
            dose_gap = dose_gap,
            influx = influx,
            cycl = cycl,
            m1 = m1, m2 = m2,
            d1 = d, d2 = d
        )[[1]]
        wins <- target_hit(data, target = 1e2, strains = c("R1", "R2"))
        summary[i, c("wins", "ymin", "ymax")] <- c(mean(wins),
            binom.test(sum(wins), length(wins))$conf.int)
        total <- total_pop(data)
        summary[i, c("pop", "pop_wmin", "pop_ymax")] <- c(mean(total),
            t.test(total)$conf.int)
        print(i / nrow(summary))
    }
    return(summary)
}

main_plot <- function(summary) {
    labels <- expand.grid(bcidal1 = 0, bcidal2 = 1, therapy = c("Cycling", "Combination"), resources = c("Abundant", "Limiting"))
    labels$label <- c("A", "B", "C", "D")
    p <- ggplot(summary, aes(x = bcidal1, y = bcidal2)) +
        geom_tile(aes(fill = 1 - wins)) +
        scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
        labs(x = expression(theta["BC, A"]), y = expression(theta["BC, B"]), fill = "P(extinct)") +
        facet_grid(rows = vars(resources), cols = vars(therapy)) +
        theme_minimal() +
        theme(
            axis.title = element_text(size = 35, face = "bold"),
            axis.text = element_text(size = 25),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 20),
            legend.position = "bottom",
            legend.key.size = unit(2, "cm"),
            strip.text = element_text(size = 25, face = "bold")  # Increase facet label sizes
        )
    if(nrow(unique(summary[, c("therapy", "resources")])) > 1){
        p <- p + geom_text(data = labels, aes(label = label), vjust = 0.7, hjust = 0.7, size = 20, fontface = "bold")
    }
    return(p)
}

mono_plot <- function(summary, series, lower, upper, ylab, text){
    p <- ggplot(summary, aes(x = bcidal1, y = get(series))) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = get(lower), ymax = get(upper))) +
        theme_light() +
        labs(
            x = expression(theta["BC, A"]),
            y = ylab
        ) +
        theme(
            plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size = 25)
        ) + 
        annotate("text", x = 0, y = Inf, label = text, hjust = 2, vjust = 1,
            size = 20, fontface = "bold")
    if(series == "wins"){
        p <- p + scale_y_continuous(limits = c(0, 1))
    } else {
        p <- p + scale_y_log10(limits = c(1, 1e11))
    }
    return(p)
}

### Figure 1
mono_high_res <- expand.grid(bcidal1 = seq(0, 1, 0.3), bcidal2 = 0,
    therapy = "Combination", resources = "Abundant")
mono_high_res <- run_sims(mono_high_res, delta = 0.6, rep = 1e1,
influx = c(A1 = 10, A2 = 0), m2 = 0)
mono1 <- mono_plot(mono_high_res, "wins", "ymin", "ymax", "Probability that resistance emerges in 100 hours", "A")
mono2 <- mono_plot(mono_high_res, "pop", "pop_wmin", "pop_ymax", "Total population after 100 hours", "B")

# print as a pdf
pdf("fig1.pdf", width = 20, height = 10)
grid.arrange(mono1, mono2, ncol = 2)
dev.off()


### Figure 2
summary <- expand.grid(bcidal1 = seq(0, 1, 1), bcidal2 = seq(0, 1, 1),
    therapy = c("Cycling", "Combination"), resources = c("Abundant", "Limiting"))
summary <- run_sims(summary, rep = 1e2)

pdf("fig2.pdf", width = 20, height = 20)
main_plot(summary)
dev.off()

### Figure S1
summary <- expand.grid(bcidal1 = seq(0, 1, 0.5), bcidal2 = seq(0, 1, 0.5),
    therapy = c("Cycling"), resources = c("Abundant")) 
cs <- run_sims(summary, rep = 1e1, zeta1 = c(S = 1, R1 = 28, R2 = 0.5, R12 = 28),
    zeta2 = c(S = 1, R1 = 0.5, R2 = 28, R12 = 28), delta = 0.25)
main_plot(cs)

### Figure S2
quick_degrade <- run_sims(summary, rep = 1e1, influx = 10 * c(A1 = 1, A2 = 1), d = 0.5, delta = 0.6)
main_plot(quick_degrade)
