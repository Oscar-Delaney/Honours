source("stochastic.R")
library(gridExtra)
library(patchwork)
library(cowplot)

# Whether target was hit
target_hit <- function(sol, target = 1e2, strains = c("N_A", "N_B")) {
    target_times <- sol %>%
      filter(variable %in% strains) %>%
      group_by(rep, time) %>%
      summarise(total = sum(value), .groups = "drop") %>%
      group_by(rep) %>%
      summarise(t = time[min(which(diff(sign(total - target)) != 0), Inf)]) %>%
      pull(t)
    !is.na(target_times)
}

run_sims <- function(summary, zeta_A = c(N_S = 1, N_A = 28, N_B = 1, N_AB = 28),
zeta_B = c(N_S = 1, N_A = 1, N_B = 28, N_AB = 28), delta = 0.25, rep = 1, dose_gap = 10,
influx = 3 * c(C_A = 1, C_B = 1), m_A = 1e-9, m_B = 1e-9, d_ = 0, init_A = 0,
init_B = 0, R0 = 1e8, data = FALSE) {
    summary$bstatic_A <- 1 - summary$bcidal_A
    summary$bstatic_B <- 1 - summary$bcidal_B
    for (i in seq_len(nrow(summary))) {
        cycl <- ifelse(summary$therapy[i] == "Cycling", TRUE, FALSE)
        res <- switch(as.character(summary$resources[i]),
            "Abundant" = 3, "Intermediate" = 1.5, "Limiting" = 0)
        d <- d_ + ifelse(cycl, 0.1, 0.35)
        sol <- simulate(
            seed = i * 1e6,
            init = c(N_S = ifelse(cycl, 5e8, 1e10),
              N_A = init_A, N_B = init_B, N_AB = 0),
            R0 = R0 * 10 ^ res,
            k = 1e8,
            alpha = 1,
            supply = 1e8,
            mu = 1,
            bcidal_A = summary$bcidal_A[i],
            bcidal_B = summary$bcidal_B[i],
            bstatic_A = summary$bstatic_A[i],
            bstatic_B = summary$bstatic_B[i],
            zeta_A = zeta_A,
            zeta_B = zeta_B,
            delta = delta + res * 0.05,
            time = 60,
            tau = 1e4,
            rep = rep,
            dose_gap = dose_gap,
            influx =  influx * (1 + !cycl),
            cycl = cycl,
            m_A = m_A, m_B = m_B,
            d_A = d, d_B = d
        )[[1]]
        wins <- 1 - target_hit(sol)
        summary[i, c("wins", "ymin", "ymax")] <- c(mean(wins),
            binom.test(sum(wins), length(wins))$conf.int)
        print(i / nrow(summary))
    }
    if (data) {
        return(sol)
    } else {
        return(summary)
    }
}

bottom_plot <- ggplot() +
  annotate("segment", x = 0, xend = 0.25, y = 0, yend = 0, arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches")), size = 0.5) +
  annotate("segment", x = 1, xend = 0.75, y = 0, yend = 0, arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches")), size = 0.5) +
  annotate("text", x = 0.15, y = 0, label = "bacteriostatic", vjust = 2, size = 6) +
  annotate("text", x = 0.85, y = 0, label = "bactericidal", vjust = 2, size = 6) +
  annotate("text", x = 0.5, y = 0, label = expression(theta["A"]), size = 10) +
  theme_void()

side_plot <- ggplot() +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 0.25, arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches")), size = 0.5) +
  annotate("segment", x = 0, xend = 0, y = 1, yend = 0.75, arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches")), size = 0.5) +
  annotate("text", x = 0, y = 0.15, label = "bacteriostatic", angle = 90, vjust = -1, size = 6) +
  annotate("text", x = 0, y = 0.85, label = "bactericidal", angle = 90, vjust = -1, size = 6) +
  annotate("text", x = 0, y = 0.5, label = expression(theta["B"]), angle = 90, size = 10) +
  theme_void()

blank_plot <- ggplot() + theme_void()

# Create a 2D tile plot of the results
main_plot <- function(summary) {
    labels <- expand.grid(bcidal_A = 0, bcidal_B = 1, therapy = unique(summary$therapy), resources = unique(summary$resources))
    labels$label <- LETTERS[1:nrow(labels)]

    p <- ggplot(summary, aes(x = bcidal_A, y = bcidal_B)) +
        geom_tile(aes(fill = wins)) +
        labs(fill = "P(extinct)") +
        theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.text = element_text(size = 25),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 20),
            legend.position = "bottom",
            legend.key.size = unit(2, "cm"),
            strip.text = element_text(size = 25, face = "bold")
        )

    if (nrow(unique(summary[, c("therapy", "resources")])) > 1) {
        all_side_plots <- side_plot / side_plot / side_plot / blank_plot +
            plot_layout(heights = c(1, 1, 1, 0.22))

        all_bottom_plots <- bottom_plot + bottom_plot

        p <- p +
            facet_grid(rows = vars(resources), cols = vars(therapy)) +
            geom_text(data = labels, aes(label = label), vjust = 1, hjust = 0, size = 15, fontface = "bold") +
            scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1))

    } else {
        all_side_plots <- side_plot / blank_plot +
            plot_layout(heights = c(1, 0.1))

        all_bottom_plots <- bottom_plot

        p <- p + scale_fill_gradient(low = "white", high = "blue")
    }
    
    data_plot <- p + theme(legend.position = "none")

    final_plot <- ((all_side_plots | ((data_plot / all_bottom_plots) +
    plot_layout(heights = c(20, 1)))) + plot_layout(widths = c(1, 20))) /
              get_legend(p) + plot_layout(heights = c(15, 1))

    return(final_plot)
}

### Figure 1
summary <- expand.grid(bcidal_A = seq(0, 1, 0.05), bcidal_B = 0,
    therapy = "Cycling", resources = "Abundant")
sol <- run_sims(summary[nrow(summary), ], rep = 1e1,
    influx = c(C_A = 6, C_B = 0), dose_gap = 5, m_B = 0, data = TRUE)
dynamics <- log_plot(sol, use = c("N_S", "N_A", "R")) +
    annotate("text", x = 0, y = Inf, label = "A", hjust = 0.5, vjust = 1.5,
        size = 15, fontface = "bold")
mono_high_res <- run_sims(summary, rep = 1e3,
    influx = c(C_A = 6, C_B = 0), dose_gap = 5, m_B = 0)
mono <- ggplot(mono_high_res, aes(x = bcidal_A, y = wins)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_light() +
    labs(
        x = NULL,
        y = "P(extinct)"
    ) +
    theme(
        axis.title = element_text(size = 35),
        axis.text = element_text(size = 25),
        plot.margin = unit(c(0, 0, 0, 2), "cm")
    ) +
    annotate("text", x = 0, y = Inf, label = "B", hjust = 1, vjust = 1.5,
        size = 15, fontface = "bold")
right <- mono / bottom_plot + plot_layout(heights = c(10, 1))
left <- dynamics / blank_plot + plot_layout(heights = c(1, 0))

# print as a pdf
pdf("bacteriostatic/fig1.pdf", width = 20, height = 10)
left | right
dev.off()

save(mono_high_res, file = "bacteriostatic/fig1.rdata")

### Figure 2
summary <- expand.grid(bcidal_A = seq(0, 1, 0.05), bcidal_B = seq(0, 1, 0.05),
    therapy = c("Combination", "Cycling"), resources = c("Abundant", "Intermediate", "Limiting"))
multi <- run_sims(summary, rep = 1e3)

pdf("bacteriostatic/fig2.pdf", width = 20, height = 25)
main_plot(multi)
dev.off()

save(multi, file = "bacteriostatic/fig2.rdata")

### Figure S1
summary <- expand.grid(bcidal_A = seq(0, 1, 0.05), bcidal_B = seq(0, 1, 0.05),
    therapy = c("Cycling"), resources = c("Abundant"))
cs <- run_sims(summary, rep = 1e3, zeta_A = c(N_S = 1, N_A = 28, N_B = 0.5, N_AB = 28),
    zeta_B = c(N_S = 1, N_A = 0.5, N_B = 28, N_AB = 28), delta = -0.05, influx = 7 * c(C_A = 1, C_B = 1))
pdf("bacteriostatic/figS1.pdf", width = 10, height = 10)
main_plot(cs)
dev.off()

save(cs, file = "bacteriostatic/figS1.rdata")

### Figure S2
quick_degrade <- run_sims(summary, rep = 1e3, influx = 30 * c(C_A = 1, C_B = 1), d_ = 0.4, delta = 0.3)

pdf("bacteriostatic/figS2.pdf", width = 10, height = 10)
main_plot(quick_degrade)
dev.off()

save(quick_degrade, file = "bacteriostatic/figS2.rdata")

### Figure S3
pre_existing <- run_sims(summary, rep = 1e3, m_A = 0, m_B = 0, init_B = 5, delta = 0.1)

pdf("bacteriostatic/figS3.pdf", width = 10, height = 10)
main_plot(pre_existing)
dev.off()

save(pre_existing, file = "bacteriostatic/figS3.rdata")