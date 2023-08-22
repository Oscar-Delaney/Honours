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
final_pop <- function(sol, strains = c("S", "RA", "RB", "RAB")) {
  sol %>%
    filter(time == max(time) & variable %in% strains) %>%
    group_by(rep) %>%
    summarise(final = sum(value)) %>%
    pull(final)
}

# Time to reach target
target_time <- function(sol, target = 1, strains = c("S", "RA", "RB", "RAB")) {
  sol %>%
    filter(variable %in% strains) %>%
    group_by(rep, time) %>%
    summarise(total = sum(value), .groups = "drop") %>%
    group_by(rep) %>%
    summarise(t = time[min(which(diff(sign(total - target)) != 0), Inf)]) %>%
    pull(t)
}

# Whether target was hit
target_hit <- function(sol, target = 1, strains = c("S", "RA", "RB", "RAB")) {
  !is.na(target_time(sol, target, strains))
}

run_sims <- function(summary, zeta1 = c(S = 1, RA = 28, RB = 1, RAB = 28),
zeta2 = c(S = 1, RA = 1, RB = 28, RAB = 28), delta = 0.3, rep = 1, dose_gap = 10,
influx = 3 * c(A = 1, B = 1), m_A = 1e-9, m_B = 1e-9, d_ = 0, data = FALSE) {
    summary$bstatic1 <- 1 - summary$bcidal1
    summary$bstatic2 <- 1 - summary$bcidal2
    for (i in seq_len(nrow(summary))) {
        cycl <- ifelse(summary$therapy[i] == "Cycling", TRUE, FALSE)
        d <- d_ + ifelse(cycl, 0.1, 0.35)
        sol <- simulate(
            init = c(S = ifelse(cycl, 5e8, 1e10), RA = 0, RB = 0, RAB = 0),
            N0 = ifelse(summary$resources[i] == "Abundant", 1e11, 1e8),
            k = 1e8,
            alpha = 1,
            supply = 1e8,
            mu = 1,
            bcidal1 = summary$bcidal1[i],
            bcidal2 = summary$bcidal2[i],
            bstatic1 = summary$bstatic1[i],
            bstatic2 = summary$bstatic2[i],
            zeta1 = zeta1,
            zeta2 = zeta2,
            delta = ifelse(summary$resources[i] == "Abundant", 0.4, 0.3),
            time = 60,
            tau = 1e4,
            rep = rep,
            dose_gap = dose_gap,
            influx =  influx * (1 + !cycl),
            cycl = cycl,
            m_A = m_A, m_B = m_B,
            d1 = d, d2 = d
        )[[1]]
        wins <- 1 - target_hit(sol, target = 1e2, strains = c("RA", "RB"))
        summary[i, c("wins", "ymin", "ymax")] <- c(mean(wins),
            binom.test(sum(wins), length(wins))$conf.int)
        # total <- total_pop(sol)
        # summary[i, c("pop", "pop_wmin", "pop_ymax")] <- c(mean(total),
        #     t.test(total)$conf.int)
        print(i / nrow(summary))
    }
    if(data){
        return(sol)
    } else {
        return(summary)
    }
}

main_plot <- function(summary) {
    labels <- expand.grid(bcidal1 = 0, bcidal2 = 1, therapy = c("Cycling", "Combination"), resources = c("Abundant", "Limiting"))
    labels$label <- c("A", "B", "C", "D")
    p <- ggplot(summary, aes(x = bcidal1, y = bcidal2)) +
        geom_tile(aes(fill = wins)) +
        scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
        labs(x = expression(theta["A"]), y = expression(theta["B"]), fill = "P(extinct)") +
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
        p <- p + geom_text(data = labels, aes(label = label), vjust = 0.7, hjust = 0.7, size = 15, fontface = "bold")
    }
    return(p)
}

mono_plot <- function(summary, series, lower, upper, ylab, text){
    p <- ggplot(summary, aes(x = bcidal1, y = get(series))) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = get(lower), ymax = get(upper))) +
        theme_light() +
        labs(
            x = expression(theta["A"]),
            y = ylab
        ) +
        theme(
            axis.title = element_text(size = 35),
            axis.text = element_text(size = 25),
            plot.margin = unit(c(0, 0, 0, 2), "cm")
        ) + 
        annotate("text", x = 0, y = Inf, label = text, hjust = 1, vjust = 1,
            size = 15, fontface = "bold")
    if(series == "wins"){
        p <- p + scale_y_continuous(limits = c(0, 1))
    } else {
        p <- p + scale_y_log10(limits = c(1, 1e11))
    }
    return(p)
}

data_dir <- "C:/Users/s4528540/Downloads/results/static"

### Figure 1
summary <- expand.grid(bcidal1 = seq(0, 1, 0.05), bcidal2 = 0,
    therapy = "Cycling", resources = "Abundant")
mono_high_res <- run_sims(summary, delta = 0.4, rep = 1e3,
    influx = c(A = 6, B = 0), dose_gap = 5, m_B = 0)
sol <- run_sims(summary[nrow(summary), ], delta = 0.4, rep = 1e1,
    influx = c(A = 6, B = 0), dose_gap = 5, m_B = 0, data = TRUE)
dynamics <- log_plot(sol, use = c("S", "RA", "RB", "RAB", "N")) +
    annotate("text", x = 0, y = Inf, label = "A", hjust = 0.5, vjust = 1,
        size = 15, fontface = "bold")
# mono1 <- mono_plot(mono_high_res, "pop", "pop_wmin", "pop_ymax", "Total S population", "B")
mono2 <- mono_plot(mono_high_res, "wins", "ymin", "ymax", "P(extinct)", "B")

# print as a pdf
pdf("bacteriostatic/fig1.pdf", width = 20, height = 10)
grid.arrange(dynamics, mono2, ncol = 2)
dev.off()

save(mono_high_res, file = paste0(data_dir, "/fig1.rdata"))

### Figure 2
summary <- expand.grid(bcidal1 = seq(0, 1, 0.05), bcidal2 = seq(0, 1, 0.05),
    therapy = c("Cycling", "Combination"), resources = c("Abundant", "Limiting"))
summary <- run_sims(summary, rep = 1e3)

pdf("bacteriostatic/fig2.pdf", width = 20, height = 20)
main_plot(summary)
dev.off()

save(summary, file = paste0(data_dir, "/fig2.rdata"))

### Figure S1
summary <- expand.grid(bcidal1 = seq(0, 1, 0.05), bcidal2 = seq(0, 1, 0.05),
    therapy = c("Cycling"), resources = c("Abundant")) 
cs <- run_sims(summary, rep = 1e3, zeta1 = c(S = 1, RA = 28, RB = 0.5, RAB = 28),
    zeta2 = c(S = 1, RA = 0.5, RB = 28, RAB = 28), delta = 0.25)

pdf("bacteriostatic/figS1.pdf", width = 10, height = 10)
main_plot(cs)
dev.off()

save(cs, file = paste0(data_dir, "/figS1.rdata"))

### Figure S2
quick_degrade <- run_sims(summary, rep = 1e3, influx = 10 * c(A = 1, B = 1), d_ = 0.4, delta = 0.6)

pdf("bacteriostatic/figS2.pdf", width = 10, height = 10)
main_plot(quick_degrade)
dev.off()

save(quick_degrade, file = paste0(data_dir, "/figS2.rdata"))