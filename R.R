library(ggplot2)
library(reshape2)

# Set parameters for the simulation
true_theta <- 0.2
sample_sizes <- c(50, 100, 200)
num_simulations <- 1000
num_bootstraps <- 100

rmodel <- function(n, theta) {
  floor(log(runif(n)) / log(1 - theta))
}

# Method of Moments estimator
mom_estimator <- function(data) {
  return(1 / mean(data))
}

# Bias-Corrected Maximum Likelihood Estimator
bcmle_estimator <- function(data) {
  n <- length(data)
  mle <- 1 / mean(data)
  bcmle <- n * mle / (n + 1)
  return(bcmle)
}

# Method of Moments estimator with bootstrapping
mom_bootstrap_estimator <- function(data, num_bootstraps) {
  n <- length(data)
  
  # Function to calculate the Method of Moments estimator
  mom_estimator <- function(sample) {
    return(1 / mean(sample))
  }
  
  # Bootstrap procedure
  bootstrap_estimates <- numeric(num_bootstraps)
  for (b in 1:num_bootstraps) {
    # Generate bootstrap sample
    bootstrap_sample <- sample(data, replace = TRUE)
    
    # Calculate Method of Moments estimator for the bootstrap sample
    bootstrap_estimates[b] <- mom_estimator(bootstrap_sample)
  }
  
  # Calculate the mean of bootstrap estimates as the final estimator
  final_estimator <- mean(bootstrap_estimates)
  
  return(final_estimator)
}

# Function to conduct a simulation study
monte_carlo_simulation <- function(estimator_function, true_theta, sample_sizes, num_simulations) {
  results <- matrix(NA, nrow = length(sample_sizes), ncol = num_simulations)
  
  for (i in 1:length(sample_sizes)) {
    n <- sample_sizes[i]
    for (j in 1:num_simulations) {
      data <- rmodel(n, true_theta)
      theta_hat <- estimator_function(data)
      results[i, j] <- theta_hat
    }
  }
  
  return(results)
}

# Function to conduct a simulation study with bootstrapping
monte_carlo_bootstrap_simulation <- function(estimator_function, true_theta, sample_sizes, num_simulations, num_bootstraps) {
  results <- matrix(NA, nrow = length(sample_sizes), ncol = num_simulations)
  
  for (i in 1:length(sample_sizes)) {
    n <- sample_sizes[i]
    for (j in 1:num_simulations) {
      data <- rmodel(n, true_theta)
      theta_hat <- estimator_function(data, num_bootstraps)
      results[i, j] <- theta_hat
    }
  }
  
  return(results)
}

# Run simulation for each estimator
mom_results <- monte_carlo_simulation(mom_estimator, true_theta, sample_sizes, num_simulations)
bcmle_results <- monte_carlo_simulation(bcmle_estimator, true_theta, sample_sizes, num_simulations)
mom_bootstrap_results <- monte_carlo_bootstrap_simulation(mom_bootstrap_estimator, true_theta, sample_sizes, num_simulations, num_bootstraps)

# Function to calculate bias
calculate_bias <- function(estimates, true_value) {
  return(mean(estimates) - true_value)
}

# Function to calculate mean squared error (MSE)
calculate_mse <- function(estimates, true_value) {
  return(mean((estimates - true_value)^2))
}

# Calculate bias and MSE for each estimator
bias_mom <- apply(mom_results, 1, calculate_bias, true_value = true_theta)
bias_bcmle <- apply(bcmle_results, 1, calculate_bias, true_value = true_theta)
bias_momboot <- apply(mom_bootstrap_results, 1, calculate_bias, true_value = true_theta)

mse_mom <- apply(mom_results, 1, calculate_mse, true_value = true_theta)
mse_bcmle <- apply(bcmle_results, 1, calculate_mse, true_value = true_theta)
mse_momboot <- apply(mom_bootstrap_results, 1, calculate_mse, true_value = true_theta)

# Calculate mean, median, and standard deviation for bias
mean_bias_mom <- mean(bias_mom)
median_bias_mom <- median(bias_mom)
sd_bias_mom <- sd(bias_mom)

mean_bias_bcmle <- mean(bias_bcmle)
median_bias_bcmle <- median(bias_bcmle)
sd_bias_bcmle <- sd(bias_bcmle)

mean_bias_momboot <- mean(bias_momboot)
median_bias_momboot <- median(bias_momboot)
sd_bias_momboot <- sd(bias_momboot)

# Calculate mean, median, and standard deviation for MSE
mean_mse_mom <- mean(mse_mom)
median_mse_mom <- median(mse_mom)
sd_mse_mom <- sd(mse_mom)

mean_mse_bcmle <- mean(mse_bcmle)
median_mse_bcmle <- median(mse_bcmle)
sd_mse_bcmle <- sd(mse_bcmle)

mean_mse_momboot <- mean(mse_momboot)
median_mse_momboot <- median(mse_momboot)
sd_mse_momboot <- sd(mse_momboot)

# Create a data frame for mean, median, and standard deviation
summary_data_separate <- data.frame(
  Estimator = rep(c("MoM", "BCMLE", "momboot"), each = 3),
  Category = rep(c("Mean", "Median", "SD"), times = 3),
  Bias = c(mean_bias_mom, median_bias_mom, sd_bias_mom,
           mean_bias_bcmle, median_bias_bcmle, sd_bias_bcmle,
           mean_bias_momboot, median_bias_momboot, sd_bias_momboot),
  MSE = c(mean_mse_mom, median_mse_mom, sd_mse_mom,
          mean_mse_bcmle, median_mse_bcmle, sd_mse_bcmle,
          mean_mse_momboot, median_mse_momboot, sd_mse_momboot)
)

# data frame for better ggplot aesthetics
summary_data_separate_melted <- melt(summary_data_separate, id.vars = c("Estimator", "Category"))

# Function to conduct a simulation study for efficiency
check_efficiency <- function(estimator_function, true_theta, sample_sizes, num_simulations) {
  efficiency_results <- matrix(NA, nrow = length(sample_sizes), ncol = num_simulations)
  
  for (i in 1:length(sample_sizes)) {
    n <- sample_sizes[i]
    for (j in 1:num_simulations) {
      data <- rmodel(n, true_theta)
      theta_hat <- estimator_function(data)
      efficiency_results[i, j] <- theta_hat
    }
  }
  
  return(efficiency_results)
}

# Run simulation for efficiency for each estimator
mom_efficiency_results <- check_efficiency(mom_estimator, true_theta, sample_sizes, num_simulations)
bcmle_efficiency_results <- check_efficiency(bcmle_estimator, true_theta, sample_sizes, num_simulations)
mom_bootstrap_efficiency_results <- check_efficiency(function(data) mom_bootstrap_estimator(data, num_bootstraps), true_theta, sample_sizes, num_simulations)

# Create a data frame for plotting
efficiency_data <- data.frame(
  Sample_Size = rep(sample_sizes, each = num_simulations),
  MOM_Efficiency = c(mom_efficiency_results),
  BCMLE_Efficiency = c(bcmle_efficiency_results),
  MOM_Bootstrap_Efficiency = c(mom_bootstrap_efficiency_results)
)

efficiency_data_long <- melt(efficiency_data, id.vars = "Sample_Size", variable.name = "Estimator", value.name = "Efficiency")

# Function to calculate Cramer-Rao Lower Bound
calculate_crlb <- function(variance, sample_size) {
  return(1 / (variance * sample_size))
}

# True parameters (assuming normal distribution for simplicity)
true_mean <- 0.2
true_variance <- 1.0

# CRLB for Method of Moments Estimator
crlb_mom <- calculate_crlb(true_variance, sample_sizes)

# CRLB for Bias-Corrected Maximum Likelihood Estimator
crlb_bcmle <- calculate_crlb(true_variance, sample_sizes)

# Display the results
cat("CRLB for MOM Estimator:", crlb_mom, "\n")
cat("CRLB for BCMLE Estimator:", crlb_bcmle, "\n")


#-------------------------------GRAPHS--------------------------------

# Plot the histogram of estimators
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (i in 1:length(sample_sizes)) {
  hist(mom_results[i, ], main = paste("MoM - n =", sample_sizes[i]), xlab = "Estimate", col = "lightblue", border = "black")
}
for (i in 1:length(sample_sizes)) {
  hist(bcmle_results[i, ], main = paste("BCMLE - n =", sample_sizes[i]), xlab = "Estimate", col = "lightgreen", border = "black")
}
for (i in 1:length(sample_sizes)) {
  hist(mom_bootstrap_results[i, ], main = paste("MoM Boot - n =", sample_sizes[i]), xlab = "Estimate", col = "orange", border = "black")
}

# Plot bias results
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
barplot(rbind(bias_mom, bias_bcmle, bias_momboot), beside = TRUE,
        col = c("lightblue", "lightgreen", "orange"),
        names.arg = sample_sizes,
        main = "Bias Comparison",
        xlab = "Sample Size", ylab = "Bias",
        #legend.text = c("MoM", "BCMLE", "MoM-Boot"),
        args.legend = list(title = "Estimator"))

# Plot mean squared error results
barplot(rbind(mse_mom, mse_bcmle, mse_momboot), beside = TRUE,
        col = c("lightblue", "lightgreen", "orange"),
        names.arg = sample_sizes,
        main = "MSE Comparison",
        xlab = "Sample Size", ylab = "MSE",
        legend.text = c("MoM", "BCMLE", "MoM-Boot"),
        args.legend = list(title = "Estimator"))

# Plot Bias and MSE results
ggplot(summary_data_separate_melted, aes(x = Category, y = value, fill = Estimator)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  labs(title = "Comparison of Estimators",
       x = "Category", y = "Value") +
  scale_fill_manual(values = c("lightblue", "lightgreen", "orange")) +
  theme_minimal()

# Plotting efficiency results using boxplots
ggplot(efficiency_data_long, aes(x = factor(Sample_Size), y = Efficiency, fill = Estimator)) +
  geom_boxplot() +
  labs(title = "Efficiency Comparison",
       x = "Sample Size",
       y = "Efficiency") +
  theme_minimal()

# QUESTION 2

rmodel = function(n, theta) {
  floor(log(runif(n)) / log(1 - theta))
}

# Wald Confidence Interval
wald_confidence_interval <- function(theta_hat, n, confidence_level = 0.90) {
  if (confidence_level < 0 || confidence_level > 1 || n <= 0) {
    stop("Invalid input values for confidence interval calculation.")
  }
  
  z <- qnorm(1 - (1 - confidence_level) / 2)
  margin_of_error <- z * sqrt(theta_hat * (1 - theta_hat) / n)
  
  lower_bound <- theta_hat - margin_of_error
  upper_bound <- theta_hat + margin_of_error
  
  result <- list(
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    margin_of_error = margin_of_error,
    confidence_level = confidence_level
  )
  
  return(result)
}

# Given Confidence Interval
variation_confidence_interval <- function(theta_hat, n, confidence_level = 0.90) {
  if (confidence_level < 0 || confidence_level > 1 || n <= 0) {
    stop("Invalid input values for confidence interval calculation.")
  }
  
  z_theta <- qnorm(1 - (1 - confidence_level) / 2)
  margin_of_error <- z_theta * theta_hat * sqrt((1 - theta_hat) / n)
  
  lower_bound <- theta_hat - margin_of_error
  upper_bound <- theta_hat + margin_of_error
  
  result <- list(
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    margin_of_error = margin_of_error,
    confidence_level = confidence_level
  )
  
  return(result)
}

# Simulation parameters
set.seed(123)
num_simulations <- 1000
sample_size <- 100
true_theta <- 0.2  # True value of the parameter

# Storage for results
wald_results <- matrix(NA, nrow = num_simulations, ncol = 2)
variation_results <- matrix(NA, nrow = num_simulations, ncol = 2)

# Simulation loop
for (i in 1:num_simulations) {
  #samples <- rmodel(sample_size, true_theta)
  samples <- rbinom(sample_size, size = 1, prob = true_theta)
  
  theta_hat <- mean(samples)
  #theta_hat <- pmax(0, pmin(1, theta_hat))
  
  # Calculate Wald Confidence Interval
  wald_ci <- wald_confidence_interval(theta_hat, sample_size)
  wald_results[i, ] <- c(wald_ci$lower_bound, wald_ci$upper_bound)
  
  # Calculate Variation of Wald Confidence Interval
  variation_ci <- variation_confidence_interval(theta_hat, sample_size)
  variation_results[i, ] <- c(variation_ci$lower_bound, variation_ci$upper_bound)
}

# Output results
cat("Wald Confidence Interval Formula: theta_hat +/- z * sqrt(theta_hat * (1 - theta_hat) / n)\n\n")
cat("Variation of Wald Confidence Interval Formula: theta_hat +/- z * theta_hat * sqrt((1 - theta_hat) / n)\n\n")

cat("Simulation Results (First 5 rows):\n")
cat("Wald Results:\n", wald_results[1:5, ], "\n\n")
cat("Variation Results:\n", variation_results[1:5, ], "\n")

# Calculate coverage probability
wald_coverage <- apply(wald_results, 1, function(row) {
  true_theta >= row[1] && true_theta <= row[2]
})

variation_coverage <- apply(variation_results, 1, function(row) {
  true_theta >= row[1] && true_theta <= row[2]
})

#-----------------DATAFRAMES------------------

# Combine coverage information into a data frame
coverage_df <- data.frame(
  Simulation = 1:num_simulations,
  Wald_Coverage = cumsum(wald_coverage) / seq_along(wald_coverage),
  Variation_Coverage = cumsum(variation_coverage) / seq_along(variation_coverage)
)


# Combine lower and upper bounds into a data frame
bounds_df <- data.frame(
  Lower_Bound_Wald = wald_results[, 1],
  Upper_Bound_Wald = wald_results[, 2],
  Lower_Bound_Variation = variation_results[, 1],
  Upper_Bound_Variation = variation_results[, 2]
)


# Combine results into data frames
wald_df <- data.frame(
  true_theta = rep(true_theta, num_simulations),
  point_estimate = colMeans(wald_results),
  lower_bound = wald_results[, 1],
  upper_bound = wald_results[, 2],
  method = rep("Wald", num_simulations)
)

variation_df <- data.frame(
  true_theta = rep(true_theta, num_simulations),
  point_estimate = colMeans(variation_results),
  lower_bound = variation_results[, 1],
  upper_bound = variation_results[, 2],
  method = rep("Variation", num_simulations)
)

# Combine data frames
combined_df <- rbind(wald_df, variation_df)

# Calculate interval widths
wald_widths <- wald_results[, 2] - wald_results[, 1]
variation_widths <- variation_results[, 2] - variation_results[, 1]

# Combine widths into a data frame
widths_df <- data.frame(
  Width = c(wald_widths, variation_widths),
  Method = rep(c("Wald", "Variation"), each = num_simulations)
)

# Calculate interval widths and coverage probability
wald_widths <- wald_results[, 2] - wald_results[, 1]
variation_widths <- variation_results[, 2] - variation_results[, 1]

wald_coverage <- apply(wald_results, 1, function(row) {
  true_theta >= row[1] && true_theta <= row[2]
})

variation_coverage <- apply(variation_results, 1, function(row) {
  true_theta >= row[1] && true_theta <= row[2]
})

# Combine information into a data frame
consistency_df <- data.frame(
  Simulation = 1:num_simulations,
  Wald_Width = wald_widths,
  Variation_Width = variation_widths,
  Wald_Coverage = cumsum(wald_coverage) / seq_along(wald_coverage),
  Variation_Coverage = cumsum(variation_coverage) / seq_along(variation_coverage)
)


#----------------------PLOTS--------------------------

library(ggplot2)

# Plotting
ggplot(widths_df, aes(x = Width, fill = Method)) +
  geom_histogram(binwidth = 0.02, position = "identity", alpha = 0.7, color = "black") +
  labs(title = "Histogram of Confidence Interval Widths",
       x = "Width",
       y = "Frequency") +
  theme_minimal()


# Histogram for lower bounds
ggplot(bounds_df, aes(x = Lower_Bound_Wald, fill = "Wald")) +
  geom_histogram(binwidth = 0.02, position = "identity", alpha = 0.7, color = "black") +
  geom_histogram(data = bounds_df, aes(x = Lower_Bound_Variation, fill = "Variation"), binwidth = 0.02, position = "identity", alpha = 0.7, color = "black") +
  labs(title = "Histogram of Confidence Intervals Lower Bounds",
       x = "Lower Bound",
       y = "Frequency") +
  scale_fill_manual(values = c("Wald" = "blue", "Variation" = "red")) +
  theme_minimal()

# Histogram for upper bounds
ggplot(bounds_df, aes(x = Upper_Bound_Wald, fill = "Wald")) +
  geom_histogram(binwidth = 0.02, position = "identity", alpha = 0.7, color = "black") +
  geom_histogram(data = bounds_df, aes(x = Upper_Bound_Variation, fill = "Variation"), binwidth = 0.02, position = "identity", alpha = 0.7, color = "black") +
  labs(title = "Histogram of Confidence Intervals Upper Bounds",
       x = "Upper Bound",
       y = "Frequency") +
  scale_fill_manual(values = c("Wald" = "blue", "Variation" = "red")) +
  theme_minimal()

# Plotting coverage probabilities
ggplot(coverage_df, aes(x = Simulation)) +
  geom_line(aes(y = Wald_Coverage), color = "blue", linetype = "solid", size = 1) +
  geom_line(aes(y = Variation_Coverage), color = "red", linetype = "solid", size = 1) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "black") +
  labs(title = "Coverage Probability of Confidence Intervals",
       x = "Simulation",
       y = "Coverage Probability") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(breaks = seq(0, num_simulations, by = 100)) +
  theme_minimal()


# Plotting interval widths
ggplot(consistency_df, aes(x = Simulation)) +
  geom_line(aes(y = Wald_Width), color = "blue", linetype = "solid", size = 1) +
  geom_line(aes(y = Variation_Width), color = "red", linetype = "solid", size = 1) +
  labs(title = "Consistency of Confidence Interval Widths",
       x = "Simulation",
       y = "Width") +
  theme_minimal()

