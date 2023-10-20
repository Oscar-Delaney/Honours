library(tidyverse)

# Load and clean the dataset
data <- read.csv("OD_data_wrangled.csv")
# there is actually no data, all just 0s, so remove those rows
data <- data[!(data$Treatment == "No drug" & data$Replicate %in% 5:8),]
# also remove the one aberrent replicate that had unexpected growth
data <- data[!(data$Treatment == "Combination" & data$Strain =="AB13" & data$Replicate == 6),]

# ensure correct data types
data$transfer <- as.integer(sub("^T", "", data$transfer))
data$Strain <- as.factor(data$Strain)
data$Treatment <- as.factor(data$Treatment)
data$Environment <- as.factor(data$Environment)
data$shift <- as.factor(substr(data$plateID, nchar(data$plateID) - 6, nchar(data$plateID)))
data$rep_factor <- as.factor(data$Replicate)

# make growth a binary variable
hist(data$OD[data$Environment != "!B" & data$Treatment != "!Combination"], breaks = 50)
# based on this histogram, choose an appropriate cutoff for growth/no growth
data$growth <- ifelse(data$OD >= 0.3, 1, 0)
# check the mean of the two groups
data %>% group_by(growth) %>% summarise(mean_growth = mean(OD))

# Check whether shift is a significant predictor of growth
t.test(growth ~ shift, data = data)

# Check whether replicate is a significant predictor of growth
model_rep <- glm(growth ~ rep_factor, data = data, family = binomial)
summary(model_rep)

# Check growth in each env
data %>% group_by(growth, Environment) %>% summarise(n = n())

# simple model for data of interest
model_str <- glm(growth ~ Strain + transfer, data = sub_data, family = binomial)
summary(model_str)