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

log_plot <- function(solutions) {
  long_data <- melt(solutions, id.vars = c("time","rep"), 
                  measure.vars = c("S", "R1", "R2", "R12"),
                  variable.name = "variable", value.name = "value")
  # Initialise the colours 
  colors <- c("black", "navy", "#800000", "#008000")
  # Create antibiotic concentrations data frame
  times <- unique(solutions$time)
  background_df <- data.frame(
    xmin = times[-length(times)],
    xmax = times[-1],
    A1 = solutions[solutions$rep == 1, ]$A1[-1] /
      max(solutions[solutions$rep == 1, ]$A1),
    A2 = solutions[solutions$rep == 1, ]$A2[-1] /
      max(solutions[solutions$rep == 1, ]$A2)
  )
  peak <- max(long_data$value, na.rm = TRUE)
  # Create the plot
  plot <- ggplot() +
    # Add the gradient backgrounds
    geom_rect(data = background_df,
      aes(xmin = xmin, xmax = xmax, ymin = peak * 10^0.2,
         ymax = peak * 10^0.4, fill = A1), color = NA, alpha = 1) +
    scale_fill_gradient(low = "white", high = colors[2],
      limits = c(0, 1), name = "A1", labels = NULL) +
    new_scale_fill() +
    geom_rect(data = background_df,
      aes(xmin = xmin, xmax = xmax, ymin = peak * 10^0.4,
        ymax = peak * 10^0.6, fill = A2), color = NA, alpha = 1) +
    scale_fill_gradient(low = "white", high = colors[3],
      limits = c(0, 1), name = "A2", labels = NULL) +
    # Add the lines
    new_scale_fill() +
    # long_data, aes(x = time, y = value, color = variable, linetype = factor(rep))
    geom_line(data = long_data, aes(x = time, y = value, color = variable, linetype = factor(rep)),
      linewidth = 1.5) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = rep("solid",max(solutions$rep))) +
    if (TRUE) {
    guides(linetype = "none") + }