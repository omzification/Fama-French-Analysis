
factors <- read_excel("~/Downloads/factors.xlsx")
portfolios <- read_excel("~/Downloads/portfolios.xlsx")

factors <- factors[,-1]
portfolios <- portfolios[,-1]  # Remove the first column from the data frame

portfolios_excess <- portfolios
for(i in 2:ncol(portfolios)) {
  portfolios_excess[[i]] <- portfolios[[i]] - factors$RF
}
{
mean_returns <- apply(portfolios_excess[,-1], 2, mean)  # Assuming portfolios_excess contains the excess returns with the first column being 'Date'
std_dev_returns <- apply(portfolios_excess[,-1], 2, sd)

n <- 756  # Sample size
mu_0 <- 0  # Population mean under the null hypothesis

# Calculate t-statistics for each portfolio
t_statistics <- (mean_returns - mu_0) / (std_dev_returns / sqrt(n))

alpha <- 0.05  # Significance level
df <- n - 1  # Degrees of freedom

# Calculate critical t-value for two-tailed test
critical_t <- qt(1 - alpha / 2, df)

# Determine which portfolios have statistically significant mean returns
significant_means <- abs(t_statistics) > critical_t

# Print the t-statistics and whether they are significant
results <- data.frame(Portfolio=names(significant_means), T_Statistic=t_statistics, Significant=significant_means)
print(results)
}

{
# Calculate summary statistics
  
mean_returns <- apply(portfolios_excess,2, mean)
variance_returns <- apply(portfolios_excess, 2, var)
skewness_returns <- apply(portfolios_excess, 2, skewness)
kurtosis_returns <- apply(portfolios_excess, 2, kurtosis)

# Assuming mean_returns, variance_returns, skewness_returns, and kurtosis_returns
# are vectors containing the corresponding statistics for each portfolio.
# Assuming mean_returns, variance_returns, skewness_returns, and kurtosis_returns
# are vectors containing the corresponding statistics for each portfolio.


# Record statistics for specific portfolios
portfolio_names <- colnames(portfolios)
portfolio_stats <- data.frame(
  Portfolio = portfolio_names,
  Mean = mean_returns,
  Variance = variance_returns,
  Skewness = skewness_returns,
  Kurtosis = kurtosis_returns
)

# Filter specific portfolios
specific_portfolios <- c('SMALL LoBM', 'SMALL HiBM', 'BIG LoBM', 'BIG HiBM')
specific_stats <- portfolio_stats[portfolio_stats$Portfolio %in% specific_portfolios, ]
print(specific_stats)


portfolio_stats$Portfolio <- factor(portfolio_stats$Portfolio, levels = portfolio_stats$Portfolio)

stats_melted <- melt(portfolio_stats, id.vars = "Portfolio")

# If `Portfolio` needs to be numbered from 1 to 25
stats_melted$Portfolio <- as.factor(seq_len(nrow(stats_melted) / length(unique(stats_melted$variable))))


ggplot(stats_melted, aes(x = Portfolio, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Portfolio", y = "Value", fill = "Statistic") +
  facet_wrap(~variable, scales = "free_y") +
  ggtitle("Summary Statistics for Each Portfolio")

# Save plots as figures in an appendix or separately
# For example, you can save them using ggsave() or other methods
}

{
  
  # Step 2: Calculate excess returns for each portfolio

  
  data_for_regression <- left_join(portfolios_excess, factors, by = "Month")
  
  results <- list()
  coefficients_matrix <- matrix(NA, nrow = 4, ncol = 25)
  colnames(coefficients_matrix) <- names(portfolios_excess)[2:26]
  rownames(coefficients_matrix) <- c("alpha", "beta_m", "beta_s", "beta_v")
  
  # Adjust the loop to properly handle column names with spaces
  for(i in 2:26) {
    # Enclose the portfolio column name in backticks to handle spaces
    col_name <- names(portfolios_excess)[i]
    formula_str <- paste0("`", col_name, "` ~ Mkt.RF + SMB + HML")
    formula <- as.formula(formula_str)
    regression <- lm(formula, data = data_for_regression)
    results[[i-1]] <- tidy(regression)
    coefficients_matrix[1, i-1] <- regression$coefficients[1] # Intercept (alpha)
    coefficients_matrix[2:4, i-1] <- regression$coefficients[2:4] # Betas
    
  }
  
  

  
  # Assuming 'coefficients_matrix' contains the regression coefficients for each portfolio
  # and 'results' contains the list of regression summaries
  
  # Convert the coefficients matrix to a tidy dataframe for plotting
  coefficients_df <- as.data.frame(t(coefficients_matrix))
  colnames(coefficients_df) <- c("alpha", "beta_m", "beta_s", "beta_v")
  coefficients_df$Portfolio <- rownames(coefficients_df)
  

  coefficients_melted <- melt(coefficients_df, id.vars = "Portfolio")
  
  # Ensure 'coefficients_df' is created as before
  coefficients_df$Portfolio <- factor(coefficients_df$Portfolio, levels = coefficients_df$Portfolio)
  
  # Continue with the melted dataframe for ggplot2
  coefficients_melted <- melt(coefficients_df, id.vars = "Portfolio")
  
  # If `Portfolio` needs to be numbered from 1 to 25
  coefficients_melted$Portfolio <- as.factor(seq_len(nrow(coefficients_melted) / length(unique(coefficients_melted$variable))))
  
  
  # Plot using ggplot2 with ordered portfolios
  ggplot(coefficients_melted, aes(x = Portfolio, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Portfolio", y = "Coefficient Value", fill = "Coefficient") +
    facet_wrap(~variable, scales = "free_y") +
    ggtitle("Regression Coefficients for Each Portfolio")
  
  
  
  # Step 5: Compute and plot residuals and estimate variance for B(p)
  residuals_matrix <- matrix(NA, nrow = length(data_for_regression$Month), ncol = 25)
  variance_matrix <- matrix(NA, nrow = 4, ncol = 25)
  
  for(i in 1:25) {
    residuals_matrix[, i] <- residuals(lm(data_for_regression[[i+1]] ~ data_for_regression$Mkt.RF + data_for_regression$SMB + data_for_regression$HML))
    variance_matrix[, i] <- diag(vcov(lm(data_for_regression[[i+1]] ~ data_for_regression$Mkt.RF + data_for_regression$SMB + data_for_regression$HML)))
  
  }
  
}

{
# Extract unique dates to loop over for cross-sectional regressions
unique_dates <- unique(data_for_regression$Month)

# Initialize a matrix to store the estimated parameters and intercept for each time t
lambda_matrix <- matrix(NA, nrow = length(unique_dates), ncol = 4)
colnames(lambda_matrix) <- c("alpha_t", "lambda_m", "lambda_s", "lambda_v")
rownames(lambda_matrix) <- as.character(unique_dates)

betas_matrix <- as.data.frame(t(coefficients_matrix))

# Perform the cross-sectional regression for each time t
for (date in unique_dates) {
  # Extract the data for the current time t
  current_data <- filter(data_for_regression, Month == date)
  
  # Assuming the structure of current_data has portfolio returns starting from the 3rd column
  Rt <- as.data.frame(t(current_data[, 2:26]))
  
  # Prepare the X matrix using the estimated betas
  X <- betas_matrix[2:4] # Assuming betas_matrix is ordered correctly
  Rt[2:4] <- X
  # Cross-sectional regression
  cross_sectional_fit <- lm(Rt$V1 ~ Rt$beta_m + Rt$beta_s + Rt$beta_v)  # The '-1' removes the intercept from the model since it's included in X
  
  # Store the estimated parameters
  lambda_matrix[as.character(date), ] <- coef(cross_sectional_fit)
}


# Convert the results into a data frame for plotting
lambda_df <- as.data.frame(lambda_matrix) %>%
  rownames_to_column("Date") %>%
  pivot_longer(-Date, names_to = "Parameter", values_to = "Estimate")

lambda_df$Date <- as.Date(paste0(lambda_df$Date, "01"), format="%Y%m%d")

ggplot(lambda_df, aes(x = Date, y = Estimate, color = Parameter)) +
  geom_line() +  # Connect points over time
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Date", y = "Estimate", color = "Parameter") +
  facet_wrap(~Parameter, scales = "free_y") +
  ggtitle("Risk Premia Over Time")


mean_lambda <- apply(lambda_matrix,2, mean)
variance_lambda <- apply(lambda_matrix, 2, var)
skewness_lambda <- apply(lambda_matrix, 2, skewness)
kurtosis_lambda <- apply(lambda_matrix, 2, kurtosis)

lambda_stats <- data.frame(mean_lambda,
                           variance_lambda,
                           skewness_lambda,
                           kurtosis_lambda)

mean_lambda
variance_lambda
skewness_lambda
kurtosis_lambda
}


{
  alpha_t_estimates <- lambda_matrix[,1]
  # Assuming alpha_t_estimates is a vector of alpha_t estimates from each cross-sectional regression
  mean_alpha <- mean(alpha_t_estimates)
  std_alpha <- sd(alpha_t_estimates)
  n <- length(alpha_t_estimates)
  se_alpha <- std_alpha / sqrt(n)
  
  # Calculate the t-statistic for mean_alpha
  t_statistic_alpha <- mean_alpha / se_alpha
  
  # Critical value for a two-tailed test at the 5% significance level
  critical_value <- qnorm(1 - 0.05 / 2)
  
  # Determine if the average pricing error is significantly different from zero
  significant <- abs(t_statistic_alpha) > critical_value
  
  # Print the result
  cat("T-statistic for average pricing error (alpha):", t_statistic_alpha, "\n")
  cat("Is the average pricing error significantly different from zero at the 5% level? ", significant, "\n")
  
  
}

{
portfolios_excess <- portfolios_excess %>% rename(Date = Month)

portfolios_excess$Date <- as.Date(paste0(portfolios_excess$Date, "01"), format="%Y%m%d")

estimation_start <- as.Date("2010-01-01")
estimation_end <- as.Date("2017-12-01")
testing_start <- as.Date("2018-01-01")
testing_end <- as.Date("2022-12-01")

# Filter the dataset for the estimation period
estimation_data <- portfolios_excess %>%
  filter(Date >= estimation_start & Date <= estimation_end)

# Calculate the sample covariance matrix
sample_cov <- cov(estimation_data[,-1])  # Assuming the first column is Date

# Shrink the covariance matrix
avg_var <- mean(diag(sample_cov))
identity_matrix <- diag(nrow = ncol(sample_cov))
shrunk_cov <- 0.5 * avg_var * identity_matrix + 0.5 * sample_cov

# Define the objective function for the optimization (minimize variance)
solve_qp <- function(cov_matrix) {
  Dmat <- cov_matrix
  dvec <- rep(0, ncol(Dmat))
  Amat <- matrix(1, nrow = 1, ncol = ncol(Dmat))
  bvec <- 1
  print(solve.QP(Dmat, dvec, t(Amat), bvec, meq = 1)$solution)
  
}

# Calculate portfolio weights using the sample and shrunk covariance matrices
weights_sample <- solve_qp(sample_cov)
weights_shrunk <- solve_qp(shrunk_cov)
weights_equal <- rep(1/ncol(sample_cov), ncol(sample_cov))  # Equally weighted

# Filter the dataset for the testing period
testing_data <- portfolios_excess %>%
  filter(Date >= testing_start & Date <= testing_end)

numeric_testing_data <- as.matrix(testing_data[,-1])

aaa <- Return.portfolio(testing_data, weights = weights_sample, verbose = TRUE, geometric = FALSE)
# Calculate excess returns using the estimated weights
excess_returns_sample <- numeric_testing_data %*% weights_sample
excess_returns_shrunk <- numeric_testing_data %*% weights_shrunk
excess_returns_equal <- numeric_testing_data %*% weights_equal

factors <- factors %>% rename(Date = Month)

factors$Date <- as.Date(paste0(factors$Date, "01"), format="%Y%m%d")

portfolios[,-1] <- portfolios[,-1] / 100

# Market excess returns as benchmark (assuming it's the first factor in `factors`)
market_excess_returns <- factors %>%
  filter(Date >= testing_start & Date <= testing_end) %>%
  pull(Mkt.RF)

# Calculate performance metrics: mean, variance, skewness, and kurtosis
calculate_metrics <- function(returns) {
  c(mean = mean(returns), variance = var(returns), skewness = skewness(returns), kurtosis = kurtosis(returns))
}

metrics_sample <- calculate_metrics(excess_returns_sample)
metrics_shrunk <- calculate_metrics(excess_returns_shrunk)
metrics_equal <- calculate_metrics(excess_returns_equal)
metrics_market <- calculate_metrics(market_excess_returns)

# Print the results
metrics <- data.frame(
  w = metrics_sample,
  ws = metrics_shrunk,
  wEW = metrics_equal,
  benchmark = metrics_market
)

weigted_returns <- data.frame(
  Date = testing_data$Date,
  Sample = excess_returns_sample / 100,
  Shrunk = excess_returns_shrunk / 100,
  Equal = excess_returns_equal / 100,
  Market = market_excess_returns / 100
)

cumulative_returnss <- data.frame(
  Date = testing_data$Date,
  w = cumsum(weigted_returns$Sample), # Back to percentages for plotting
  ws = cumsum(weigted_returns$Shrunk),
  wEW = cumsum(weigted_returns$Equal),
  Benchmark = cumsum(weigted_returns$Market)
)


ggplot(melt(cumulative_returnss, id = "Date"), aes(x = Date, y = value, color = variable)) +
  geom_line() +
  labs(title = "Cumulative (monthly compounded) Excess Returns Over Time", x = "Date", y = "Cumulative Excess Return", color = "Portfolio") +
  theme_minimal()


}
print(sd(excess_returns_sample))
print(sqrt(var(excess_returns_sample)))
{
  compute_sharpe_ratio <- function(portfolio) {
    
    mean_excess_returns <- metrics['mean',portfolio]
    print(mean_excess_returns)
    sd_excess_returns <- metrics['variance',portfolio]
    print(sd_excess_returns)
    sharpe_ratio <- (sqrt(12) * mean_excess_returns / sqrt(sd_excess_returns))
    print(sharpe_ratio)
    return(sharpe_ratio)
  }
  
  sharpe_ratio_sample <- compute_sharpe_ratio('w')
  sharpe_ratio_shrunk <- compute_sharpe_ratio('ws')
  sharpe_ratio_equal <- compute_sharpe_ratio('wEW')
  sharpe_ratio_market <- compute_sharpe_ratio('benchmark')
  
  print(sharpe(excess_returns_sample, r=0))
  
  print(apply(excess_returns_shrunk, 2, var))
  
  print(sharpe_ratio_shrunk / sqrt( 1 + 0.5 * sharpe_ratio_shrunk * sharpe_ratio_shrunk / 60))
  
  significance_test <- function(sharpe_ratio) {
  n <- 60  # Sample size
  mu_0 <- 0  # Population mean under the null hypothesis
  sd <- sqrt( 1 + 0.5 * sharpe_ratio_shrunk * sharpe_ratio_shrunk / 60)
  
  # Calculate t-statistics for each portfolio
  t_statistic <- sharpe_ratio / sd
  
  alpha <- 0.1  # Significance level
  df <- n - 1  # Degrees of freedom
  
  # Calculate critical t-value for two-tailed test
  critical_t <- qt(1 - alpha, df)
  
  print(abs(t_statistic) > critical_t)
  
  }
  significance_test(sharpe_ratio_sample)
  significance_test(sharpe_ratio_shrunk)
  significance_test(sharpe_ratio_equal)
  significance_test(sharpe_ratio_market)
  
}
