# setwd("/home/s4528540/Honours")
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

# Run many simulations in sequence with set parameters
run_sims <- function(summary, config_only = FALSE, rep = 1e3, zeta = 1e9,
c = 2, kappa = 1, cost = 0, net = 0, d = 0, gap = 1e4, zeta_rand = FALSE) {
    pb <- txtProgressBar(min = 0, max = nrow(summary), style = 3)
    for (i in seq_len(nrow(summary))) {
        m_A <- 1 / (1 + 1 / summary$m_ratio[i])
        c_A <- 1 / (1 + 1 / summary$c_ratio[i])
        data <- simulate(
            init = c(N_S = init_S, N_A = 0, N_B = 0, N_AB = 0),
            R0 = 1e10,
            k = 0,
            alpha = 0,
            supply = 1e8,
            mu = c(1, 1 - cost, 1 - cost, 1 - 2 * cost),
            bcidal_A = summary$cidal_A[i],
            bcidal_B = summary$cidal_B[i],
            bstatic_A = 1 - summary$cidal_A[i],
            bstatic_B = 1 - summary$cidal_B[i],
            zeta_rand = zeta_rand,
            zeta_A = c(N_S = 1, N_A = zeta, N_B = 1, N_AB = zeta),
            zeta_B = c(N_S = 1, N_A = 1, N_B = zeta, N_AB = zeta),
            delta = 1 - 1 / (1 + c ^ -kappa) - net,
            time = 100 * 5 ^ ((1 - summary$cidal_A[i]) * (1 - summary$cidal_B[i])),
            tau = 1e4,
            max_step = 1e-1,
            kappa_A = kappa, kappa_B = kappa,
            rep = rep,
            HGT = 0,
            dose_rep = 1,
            dose_gap = gap,
            influx = c * c(C_A = c_A, C_B = 1 - c_A),
            cycl = FALSE,
            m_A = m * m_A, m_B = m * (1 - m_A),
            d_A = d, d_B = d,
            deterministic = FALSE,
            config_only = config_only
        )
        if (config_only) {
            n <- rep ^ zeta_rand
            phi_values <- numeric(n)
            N_values <- numeric(n)
            extinct_values <- numeric(n)
            zeta_A <- data$zeta_A
            zeta_B <- data$zeta_B
            for (j in 1:n) {
                if (zeta_rand == TRUE) {
                    data$zeta_A["N_A"] <- 1 + qexp(p = j / (n + 1), rate = 1 / (zeta_A["N_A"] - 1))
                    data$zeta_B["N_B"] <- 1 + qexp(p = j / (n + 1), rate = 1 / (zeta_B["N_B"] - 1))
                    # zeta_randomised <- randomise_zeta(zeta_A, zeta_B)
                    # data$zeta_A <- zeta_randomised[[1]]
                    # data$zeta_B <- zeta_randomised[[2]]
                }
                v <- with(data, as.list(rates(c(N_S = 1, N_A = 1, N_B = 1,
                    N_AB = 1, R = R0, influx * pattern), data, 0)))
                phi_values[j] <- with(v, m_A * pmin(1, N_A_death / N_A_growth) +
                    (1 - m_A) * pmin(1, N_B_death / N_B_growth))
                N_values[j] <- with(v, init_S / (pmax(1e-30, S_death / S_growth - 1)))
                extinct_values[j] <- exp(-m * (1 - phi_values[j]) * N_values[j])
            }
            summary$phi[i] <- mean(phi_values)
            summary$N[i] <- mean(N_values)
            summary$extinct[i] <- mean(extinct_values)
        } else {
        wins <- !target_hit(data[[1]], target = 1e2, strains = c("N_A", "N_B"))
        summary$extinct[i] <- mean(wins)
        ci <- binom.test(sum(wins), length(wins), conf.level = 0.95)$conf.int
        summary$ymin[i] <- ci[1]
        summary$ymax[i] <- ci[2]
        }
        if (i %% 50 == 0) {setTxtProgressBar(pb, i)}
    }
    return(summary)
}

# Create a grid of parameter values to test
create_grid <- function(length = 10) {
  expand.grid(c_ratio = 2 ^ seq(-5, 5, length.out = length),
              m_ratio = 2 ^ seq(-10, 10, length.out = length),
              cidal_A = c(0, 1),
              cidal_B = c(0, 1))
}

# plot with two independent variables and one dependent variable
summary_plot <- function(fine, coarse, var = "extinct", theory = TRUE) {
    # Find the y-values that maximize 'var' for each x-value
    max_df <- fine
    if (theory) {
        max_df <- fine %>%
            group_by(cidal_A, cidal_B, m_ratio) %>%
            arrange(desc(get(var))) %>%
            slice_head(n = 1) %>%
            ungroup()
    }
    max_df$theory <- max_df$m_ratio ^ -0.5

    # Rename the factor levels with the desired prefixes
    max_df$cidal_A <- ifelse(max_df$cidal_A == 0, "Drug A: Bacteriostatic", "Drug A: Bactericidal")
    max_df$cidal_B <- ifelse(max_df$cidal_B == 0, "Drug B: Bacteriostatic", "Drug B: Bactericidal")
    coarse$cidal_A <- ifelse(coarse$cidal_A == 0, "Drug A: Bacteriostatic", "Drug A: Bactericidal")
    coarse$cidal_B <- ifelse(coarse$cidal_B == 0, "Drug B: Bacteriostatic", "Drug B: Bactericidal")
    labels <- expand.grid(m_ratio = 500, c_ratio = 19, cidal_A = unique(coarse$cidal_A), cidal_B = unique(coarse$cidal_B))
    labels$label <- LETTERS[nrow(labels):1]

    p <- ggplot(coarse, aes(x = log2(m_ratio), y = log2(c_ratio))) +
        geom_tile(aes(fill = (get(var)))) +
        geom_line(data = max_df, aes(x = log2(m_ratio), y = log2(theory)), color = "green", linewidth = 3) +
        facet_grid(cidal_B ~ cidal_A) +
        geom_text(data = labels, aes(label = label), size = 15, fontface = "bold") +
        labs(fill = "P(extinct)") +
        theme_minimal() +
        labs(
            x = "Log2 odds of A-resistance:B-resistance",
            y = "Log2 odds of A-molecule:B-molecule"
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

    if (theory) {
        p <- p + geom_line(data = max_df, aes(x = log2(m_ratio), y = log2(c_ratio)), color = "yellow", linewidth = 3)
    }

    return(p)
}

# Automate the data creation and visualisation workflow
run_and_save <- function(name, args = list(), theory = TRUE, sims = TRUE) {
  fine_sims <- fine
  # Run simulations
  if (theory) {
    fine_sims <- do.call(run_sims, c(list(fine, config_only = TRUE), args))
  }
  coarse_sims <- fine_sims
  if (sims) {
    coarse_sims <- do.call(run_sims, c(list(coarse, config_only = FALSE), args))
  }
  # Save the plot and data
  p <- summary_plot(fine_sims, coarse_sims, theory = theory)
  ggsave(p, filename = paste0(dir, name, ".pdf"), width = 20, height = 20)
  save(fine_sims, coarse_sims, file = paste0(dir, name, ".rdata"))
}

# defines some static objects
dir <- "mutation_rate_variation/figs/"
init_S <- 1e9
m <- 1e-9
fine <- create_grid(200)
coarse <- create_grid(20)

# Run the simulations and create graphs for each scenario
run_and_save("basic")
run_and_save("in_res", args = list(zeta = 5, zeta_rand = TRUE))
# run_and_save("kappa_high", args = list(c = 5, kappa = 3))
run_and_save("kappa_low", args = list(kappa = 0.2))
run_and_save("costs", args = list(cost = 0.1))
# run_and_save("net", args = list(net = -0.1))
run_and_save("net_kappa", args = list(net = -0.1, kappa = 3))
run_and_save("pk", args = list(gap = 12, d = 0.15, net = -0.2), theory = FALSE)