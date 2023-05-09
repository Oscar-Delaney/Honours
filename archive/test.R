# Load the matrixStats package
library(matrixStats)

# Create a matrix from the list
value <- c(13.76145, 69.00000, 65.00000, 22.64049)
print(
    matrixStats::rowQuantiles(matrix(value, nrow = 1), probs = c(0.25))[[1]]
)
data_matrix <- matrix(data_list, nrow = 1)

# Calculate the 25th and 75th percentiles using matrixStats::rowQuantiles
percentiles <- matrixStats::rowQuantiles(data_matrix, probs = c(0.25))

# Access the 25th and 75th percentiles
quantile_25 <- percentiles[1]
quantile_75 <- percentiles[2]

# Print the 25th and 75th percentiles
print(quantile_25)
print(quantile_75)


# Unused code
IQR_bounds = list(colQuantiles(as.matrix(value), probs = c(0.25, 0.5, 0.75), na.rm = TRUE))

num_cores <- detectCores()
solutions <- tryCatch(
{
    mclapply(1:rep, function(x) {single_run(config, x)}, mc.cores = num_cores)
},
error = function(e) {
    lapply(1:rep, function(x) {single_run(config, x)})
}
)


  print(length(solution$time))
  duplicated_time_values <- solution$time[duplicated(solution$time)]
  print(duplicated_time_values)


# adaptivetau investigation

transitions <- list(
    s = c(y1 = +1),
    d = c(y2 = +1),
    u = c(y3 = +1)
)

rates <- function(state,parms,t) {
    return(c(
        s = state["y1"],
        d = state["y2"],
        u = 0
        )
    )
}

state <- c(y1 = 1e0, y2 = 1e0, y3=0)
parms <- 0
sol <- as.data.frame(ssa.adaptivetau(state, transitions, rates,
    parms, tf = 10, deterministic = c(F,T,F)))
# sol$y3 <- exp(sol$time)
ggplot(sol, aes(x = time)) +
    geom_step(aes(y = y1, color = "y1")) +
    geom_step(aes(y = y2, color = "y2")) +
    # geom_step(aes(y = y3, color = "y3")) +
    scale_color_manual(name = "Species", values = c("y1" = "blue", "y2" = "red")) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      breaks = 10^seq(0, 20),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_minimal() +
    labs(title = "Adaptive tau method",
        x = "Time",
        y = "Population size") +
    theme(
        legend.position = "bottom",
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
    )
# plot(sol$time,sol$y1)
sol
