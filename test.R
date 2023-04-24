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
