# Step 1: Load necessary libraries
library(tidyverse)
library(car)
library(lmtest)

# Step 2: Load the dataset
data <- read.csv("OD_data_wrangled.csv")
data <- data[!(data$Treatment == "No drug" & data$Replicate %in% c(7,8)),]
data$transfer <- as.integer(sub("^T", "", data$transfer))
data$Strain <- as.factor(data$Strain)
data$Treatment <- as.factor(data$Treatment)
data$Environment <- as.factor(data$Environment)
data$growth <- ifelse(data$OD >= 0.3, 1, 0)
# Step 3: Preliminary data exploration
head(data)
summary(data)

# Check for any missing values
if(any(is.na(data))){
  print("Data contains missing values")
} else {
  print("No missing values found")
}

# Step 4: Build a linear regression model
# Note: This is a basic linear regression model. If you need more complex models 
# like mixed effects, interactions, or transformations, those will need to be added.
model <- lm(growth ~ Strain + transfer + Treatment + Replicate + Environment, data = data)
model <- lm(growth ~ Strain * Environment + Treatment + Replicate + transfer, data = data)
model <- glm(growth ~ Treatment, data = data[data$Environment == "B", ], family = binomial)
# plot(data$Replicate, data$OD, col = data$Strain)
# Step 5: Check assumptions of the regression
# Residuals vs Fitted values
plot(model, which = 1)

# Normal Q-Q plot
plot(model, which = 2)

# Scale-Location plot
plot(model, which = 3)

# Residuals vs Leverage
plot(model, which = 5)

# Check for heteroscedasticity
bptest(model)

# Check for multicollinearity
vif(model) # Variance Inflation Factor

# Step 6: Summarize results
summary(model)

data[data$Treatment == "Combination" & data$Environment == "B", "growth"]
