# Load required libraries
library(tidyverse)
library(Matrix)

# bayesian VAR for feature extraction

dynamic_df <- read_csv('/Users/dwaipayanchakrabarti/Downloads/final_data.csv')

# variables
# gdp growth
# employment growth
# composite rate
# financial stress
# epu

gdp_growth <- dynamic_df$gdp_growth
employment_growth <- dynamic_df$employment_growth
EPU_value_Quarterly <- dynamic_df$EPU_value_Quarterly
composite_rates <- dynamic_df$composite_rates
growth_STLFSI4 <- dynamic_df$growth_STLFSI4

library(dplyr)

state_vars    <- c("gdp_growth", "employment_growth", "EPU_value_Quarterly")
national_vars <- c("composite_rates", "growth_STLFSI4")

# standardizing
df_standard <- dynamic_df %>%
  group_by(State) %>%
  mutate(across(all_of(state_vars),
                ~ (. - mean(., na.rm=TRUE)) / sd(., na.rm=TRUE),
                .names = "{.col}_std")) %>%
  ungroup() %>%
  
  mutate(across(all_of(national_vars),
                ~ (. - mean(., na.rm=TRUE)) / sd(., na.rm=TRUE),
                .names = "{.col}_std"))

# Read in the data
df_standard

# Define variables (already standardized in the dataset)
state_vars_std <- c("gdp_growth_std", "employment_growth_std", "EPU_value_Quarterly_std")
national_vars_std <- c("composite_rates_std", "growth_STLFSI4_std")
all_vars_std <- c(state_vars_std, national_vars_std)

# Original variable names (without _std)
state_vars <- c("gdp_growth", "employment_growth", "EPU_value_Quarterly")
national_vars <- c("composite_rates", "growth_STLFSI4")
all_vars <- c(state_vars, national_vars)

# Function to parse year and quarter from YearQuarter string
parse_year_quarter <- function(year_quarter) {
  year <- as.integer(substr(year_quarter, 1, 4))
  quarter <- as.integer(substr(year_quarter, 6, 6))
  return(list(year = year, quarter = quarter))
}

# Function to calculate overall fit measures for VAR models
calculate_overall_fit <- function(Y_dep, X_matrix, B_post, n_vars, k) {
  # Y_dep: Dependent variable matrix
  # X_matrix: Matrix of lagged variables and constant
  # B_post: Posterior coefficient matrix
  # n_vars: Number of variables in the system
  # k: Number of coefficients per equation
  
  # Calculate fitted values
  Y_fitted <- X_matrix %*% B_post
  
  # Calculate residuals
  residuals <- Y_dep - Y_fitted
  
  # Sample size
  n <- nrow(X_matrix)
  
  # Total number of parameters in the model
  total_params <- n_vars * k
  
  # Calculate multivariate R-squared (McElroy's R-squared)
  # This is based on the determinants of residual and total covariance matrices
  
  # Residual covariance matrix
  Sigma_res <- (t(residuals) %*% residuals) / n
  
  # Calculate total variation (around means)
  Y_demeaned <- Y_dep - matrix(rep(colMeans(Y_dep), each = n), nrow = n)
  Sigma_tot <- (t(Y_demeaned) %*% Y_demeaned) / n
  
  # To avoid singularity issues
  det_res <- det(Sigma_res)
  det_tot <- det(Sigma_tot)
  
  if (det_tot == 0 || is.na(det_tot) || is.nan(det_tot)) {
    R2_overall <- NA
    R2_adj <- NA
  } else {
    # McElroy's R-squared
    R2_overall <- 1 - det_res / det_tot
    
    # Adjusted R-squared for multivariate model
    R2_adj <- 1 - (1 - R2_overall) * ((n - 1) / (n - total_params - 1))
    
    # Cap R-squared at 1 (sometimes numerical issues can cause it to exceed 1)
    R2_overall <- min(R2_overall, 1)
    R2_adj <- min(R2_adj, 1)
  }
  
  # Return all fit measures
  return(list(
    R2_overall = R2_overall,
    R2_adjusted = R2_adj
  ))
}

# Function to select optimal lag based on BIC
select_optimal_lag <- function(data, max_lag = 4) {
  n_obs <- nrow(data)
  n_vars <- ncol(data)
  Y <- as.matrix(data)
  
  # If we don't have enough observations for even 1 lag, return NULL
  if (n_obs <= max_lag + 2) {
    return(NULL)
  }
  
  # Initialize BIC values
  bic_values <- numeric(max_lag)
  
  # Loop through possible lags
  for (p in 1:max_lag) {
    # Skip if we don't have enough observations
    if (n_obs <= p + 2) {
      bic_values[p] <- Inf
      next
    }
    
    # Create X matrix (lagged values and constant)
    X <- matrix(NA, nrow = n_obs - p, ncol = p * n_vars + 1)
    X[, 1] <- 1  # Constant
    
    # Fill in lagged values
    for (t in (p+1):n_obs) {
      for (lag in 1:p) {
        lag_idx <- (lag-1) * n_vars + 2
        X[t-p, lag_idx:(lag_idx+n_vars-1)] <- Y[t-lag, ]
      }
    }
    
    # Dependent variable
    Y_dep <- Y[(p+1):n_obs, ]
    
    # OLS estimates
    tryCatch({
      XtX <- t(X) %*% X
      XtY <- t(X) %*% Y_dep
      B_ols <- solve(XtX, XtY)
      
      # Residuals
      residuals <- Y_dep - X %*% B_ols
      SSE <- t(residuals) %*% residuals
      
      # Log likelihood
      n <- nrow(X)
      log_likelihood <- -0.5 * n * n_vars * log(2 * pi) - 0.5 * n * log(det(SSE / n)) - 0.5 * n * n_vars
      
      # BIC
      k <- ncol(X) * n_vars  # Number of parameters
      bic_values[p] <- -2 * log_likelihood + log(n) * k
    }, error = function(e) {
      bic_values[p] <- Inf
    })
  }
  
  # Find optimal lag (minimum BIC)
  optimal_lag <- which.min(bic_values)
  
  # If all BICs are Inf, return 1 as default
  if (all(is.infinite(bic_values))) {
    return(1)
  }
  
  return(optimal_lag)
}

# Simple Bayesian VAR with Minnesota prior and flexible lag selection
fit_state_bvar_with_optimal_lag <- function(state_data, max_lag = 4) {
  # Create a copy of the data
  df_state <- state_data
  
  # Define the variables to include in the VAR (already standardized)
  var_cols <- c(
    "gdp_growth_std", 
    "employment_growth_std", 
    "EPU_value_Quarterly_std",
    "composite_rates_std", 
    "growth_STLFSI4_std"
  )
  
  # Check if these columns exist
  if(!all(var_cols %in% colnames(df_state))) {
    missing_cols <- var_cols[!var_cols %in% colnames(df_state)]
    cat("  Missing columns for state", df_state$State[1], ":", paste(missing_cols, collapse=", "), "\n")
    return(NULL)
  }
  
  # Get the data
  data <- df_state[, var_cols]
  
  # Handle missing values
  missing_idx <- apply(data, 1, function(x) any(is.na(x)))
  if (all(missing_idx)) {
    cat("  All data rows contain missing values for state", df_state$State[1], "\n")
    return(NULL)
  }
  
  data <- data[!missing_idx, ]
  
  # Check if we have enough data
  if (nrow(data) <= 5) {  # Need at least a few observations
    cat("  Not enough data for state", df_state$State[1], "\n")
    return(NULL)
  }
  
  # Determine optimal lag length for this state
  optimal_lag <- select_optimal_lag(data, max_lag)
  
  # If we couldn't determine optimal lag, use 1 as default
  if (is.null(optimal_lag)) {
    cat("  Not enough observations for lag selection, using default lag=1 for state", df_state$State[1], "\n")
    optimal_lag <- 1
  } else {
    cat("  Selected optimal lag =", optimal_lag, "for state", df_state$State[1], "\n")
  }
  
  # Set lag length based on optimal selection
  p <- optimal_lag
  n_vars <- ncol(data)
  n_obs <- nrow(data)
  
  # Check if we have enough observations after lagging
  if (n_obs <= p + 2) {
    cat("  Not enough observations after using optimal lag, for state", df_state$State[1], "\n")
    return(NULL)
  }
  
  # Convert to matrix
  Y <- as.matrix(data)
  
  # Create lagged data matrix (X)
  X <- matrix(NA, nrow = n_obs - p, ncol = p * n_vars + 1) # +1 for constant
  X[, 1] <- 1  # Constant term
  
  # Fill in lagged values
  for (t in (p+1):n_obs) {
    for (lag in 1:p) {
      lag_idx <- (lag-1) * n_vars + 2  # +2 to account for constant
      X[t-p, lag_idx:(lag_idx+n_vars-1)] <- Y[t-lag, ]
    }
  }
  
  # Dependent variable matrix (current values)
  Y_dep <- Y[(p+1):n_obs, ]
  
  # Number of coefficients per equation
  k <- ncol(X)
  
  # MN prior hyperparameters
  lambda1 <- 0.1  # Overall tightness
  lambda2 <- 0.5  # Cross-variable tightness
  lambda3 <- 1.0  # Lag decay
  lambda4 <- 100  # Constant term
  
  # Prior mean matrix (k x n_vars)
  B_prior <- matrix(0, nrow = k, ncol = n_vars)
  # First own lag coefficients have prior mean 0.9, decaying for longer lags
  for (i in 1:n_vars) {
    for (lag in 1:p) {
      idx <- (lag-1) * n_vars + i + 1  # +1 for constant
      if (lag == 1) {
        B_prior[idx, i] <- 0.9
      } else {
        B_prior[idx, i] <- 0.9 / (lag^2)  # Decay for higher lags
      }
    }
  }
  
  # Prior variance
  # For each equation i and coefficient j
  V_prior <- matrix(0, nrow = k, ncol = n_vars)
  
  # Fill in prior variances
  for (i in 1:n_vars) {  # For each equation
    # Constant term
    V_prior[1, i] <- lambda4
    
    # For lagged terms
    for (lag in 1:p) {
      for (j in 1:n_vars) {
        # Index in the coefficient vector
        coef_idx <- (lag-1) * n_vars + j + 1  # +1 for constant
        
        if (j == i) {  # Own lag
          V_prior[coef_idx, i] <- lambda1 / (lag^lambda3)
        } else {  # Cross-variable lag
          # Use relative scale of variables - here simplified
          V_prior[coef_idx, i] <- (lambda1 * lambda2) / (lag^lambda3)
        }
      }
    }
  }
  
  # Convert to diagonal matrices for computational efficiency
  V_prior_diag <- lapply(1:n_vars, function(i) diag(V_prior[, i]))
  
  # OLS estimates
  XtX <- t(X) %*% X
  XtY <- t(X) %*% Y_dep
  B_ols <- solve(XtX, XtY)
  
  # Posterior calculation for each equation
  B_post <- matrix(0, nrow = k, ncol = n_vars)
  V_post_inv <- list()
  V_post <- list()
  
  for (i in 1:n_vars) {
    # Posterior precision
    V_post_inv[[i]] <- t(X) %*% X + solve(V_prior_diag[[i]])
    V_post[[i]] <- solve(V_post_inv[[i]])
    
    # Posterior mean
    B_post[, i] <- V_post[[i]] %*% (solve(V_prior_diag[[i]], B_prior[, i]) + 
                                      t(X) %*% Y_dep[, i])
  }
  
  # Standard errors from posterior variance
  SE_post <- matrix(0, nrow = k, ncol = n_vars)
  for (i in 1:n_vars) {
    SE_post[, i] <- sqrt(diag(V_post[[i]]))
  }
  
  # Create results dataframe
  coef_df <- data.frame()
  
  # Variable names (without _std suffix)
  var_names <- c("gdp_growth", "employment_growth", "EPU_value_Quarterly", 
                 "composite_rates", "growth_STLFSI4")
  
  # For each equation (dependent variable)
  for (i in 1:n_vars) {
    # For each right-hand side variable
    for (j in 2:k) {  # Skip constant (j=1)
      # Determine which variable and lag this coefficient represents
      # Calculate lag index and variable index
      relative_idx <- j - 1  # Adjust for constant term
      lag_idx <- ceiling(relative_idx / n_vars)
      var_idx <- ((relative_idx - 1) %% n_vars) + 1
      
      # Calculate Bayesian p-value (probability of opposite sign)
      coef_val <- B_post[j, i]
      coef_std <- SE_post[j, i]
      
      # Probability mass on other side of zero
      if (coef_val >= 0) {
        p_value <- pnorm(0, mean = coef_val, sd = coef_std)
      } else {
        p_value <- 1 - pnorm(0, mean = coef_val, sd = coef_std)
      }
      
      # Convert to two-tailed p-value
      p_value <- min(p_value * 2, 1.0)
      
      # Calculate 95% credible interval
      ci_lower <- coef_val - 1.96 * coef_std
      ci_upper <- coef_val + 1.96 * coef_std
      
      # Add to dataframe
      coef_row <- data.frame(
        State = df_state$State[1],
        Lag_Order = p,
        Equation = var_names[i],
        Variable = paste0(var_names[var_idx], "_lag", lag_idx),
        Coefficient = coef_val,
        Std_Error = coef_std,
        P_Value = p_value,
        CI_Lower = ci_lower,
        CI_Upper = ci_upper,
        Significant = p_value < 0.05
      )
      coef_df <- rbind(coef_df, coef_row)
    }
  }
  
  # Calculate model fit metrics
  
  # 1. Log likelihood
  residuals <- Y_dep - X %*% B_post
  SSE <- t(residuals) %*% residuals
  n <- nrow(X)
  log_likelihood <- -0.5 * n * n_vars * log(2 * pi) - 0.5 * n * log(det(SSE / n)) - 0.5 * n * n_vars
  
  # 2. AIC (approximation for Bayesian model)
  aic <- -2 * log_likelihood + 2 * (k * n_vars)
  
  # 3. BIC (approximation for Bayesian model)
  bic <- -2 * log_likelihood + log(n) * (k * n_vars)
  
  # 4. DIC (Deviance Information Criterion)
  # This is a rough approximation - a full implementation would require MCMC samples
  dic <- -2 * log_likelihood + 2 * (k * n_vars)
  
  # 5. R-squared for each equation
  r_squared <- numeric(n_vars)
  for (i in 1:n_vars) {
    SS_total <- sum((Y_dep[, i] - mean(Y_dep[, i]))^2)
    SS_residual <- sum(residuals[, i]^2)
    r_squared[i] <- 1 - SS_residual / SS_total
  }
  
  # 6. Calculate overall fit measures
  overall_fit <- calculate_overall_fit(Y_dep, X, B_post, n_vars, k)
  
  # Store model fit metrics
  fit_metrics <- list(
    lag_order = p,
    log_likelihood = as.numeric(log_likelihood),
    aic = as.numeric(aic),
    bic = as.numeric(bic),
    dic = as.numeric(dic),
    r_squared = r_squared,
    R2_overall = overall_fit$R2_overall,
    R2_adjusted = overall_fit$R2_adjusted
  )
  
  # Return both coefficients and fit metrics
  return(list(
    coefficients = coef_df,
    fit_metrics = fit_metrics
  ))
}

# Create standardized columns if they don't exist
if (!all(all_vars_std %in% colnames(df_standard))) {
  cat("Creating standardized columns that don't exist...\n")
  
  # For each state, standardize state-specific variables
  for (state in unique(df_standard$State)) {
    state_data <- df_standard[df_standard$State == state, ]
    
    for (var in state_vars) {
      std_var <- paste0(var, "_std")
      if (!std_var %in% colnames(df_standard)) {
        # Check if we have enough data
        if (sum(!is.na(state_data[[var]])) > 1) {
          # Calculate mean and sd by state
          mean_val <- mean(state_data[[var]], na.rm = TRUE)
          sd_val <- sd(state_data[[var]], na.rm = TRUE)
          
          if (sd_val > 0) {
            df_standard[df_standard$State == state, std_var] <- 
              (state_data[[var]] - mean_val) / sd_val
          } else {
            df_standard[df_standard$State == state, std_var] <- 0
          }
        } else {
          df_standard[df_standard$State == state, std_var] <- NA
        }
      }
    }
  }
  
  # Standardize national variables across all states
  for (var in national_vars) {
    std_var <- paste0(var, "_std")
    if (!std_var %in% colnames(df_standard)) {
      # Check if we have enough data
      if (sum(!is.na(df_standard[[var]])) > 1) {
        # Calculate mean and sd across all states
        mean_val <- mean(df_standard[[var]], na.rm = TRUE)
        sd_val <- sd(df_standard[[var]], na.rm = TRUE)
        
        if (sd_val > 0) {
          df_standard[[std_var]] <- (df_standard[[var]] - mean_val) / sd_val
        } else {
          df_standard[[std_var]] <- 0
        }
      } else {
        df_standard[[std_var]] <- NA
      }
    }
  }
}

# Process all states with optimal lag selection
all_states_var <- function(df, max_lag = 4) {
  # Get list of states
  states <- unique(df$State)
  
  # Initialize empty data frames for results
  all_coefs <- data.frame()
  all_metrics <- data.frame()
  
  for (state in states) {
    cat("Processing state:", state, "\n")
    
    # Get data for this state
    state_data <- df[df$State == state, ]
    
    # Sort by YearQuarter
    state_data <- state_data[order(state_data$YearQuarter), ]
    
    # Fit the model with optimal lag selection
    tryCatch({
      result <- fit_state_bvar_with_optimal_lag(state_data, max_lag)
      
      if (!is.null(result)) {
        # Add coefficient results
        all_coefs <- rbind(all_coefs, result$coefficients)
        
        # Add fit metrics
        metrics_df <- data.frame(
          State = state,
          Lag_Order = result$fit_metrics$lag_order,
          LogLikelihood = result$fit_metrics$log_likelihood,
          AIC = result$fit_metrics$aic,
          BIC = result$fit_metrics$bic,
          DIC = result$fit_metrics$dic,
          R2_overall = result$fit_metrics$R2_overall,
          R2_adjusted = result$fit_metrics$R2_adjusted
        )
        
        # Add R-squared for each equation
        for (i in 1:length(result$fit_metrics$r_squared)) {
          var_name <- c("gdp_growth", "employment_growth", "EPU_value_Quarterly", 
                        "composite_rates", "growth_STLFSI4")[i]
          metrics_df[[paste0("R2_", var_name)]] <- result$fit_metrics$r_squared[i]
        }
        
        all_metrics <- rbind(all_metrics, metrics_df)
        
        cat("  Successfully fit model for", state, "with lag order =", result$fit_metrics$lag_order, "\n")
        cat("  Overall R-squared =", round(result$fit_metrics$R2_overall, 4), 
            ", Adjusted R-squared =", round(result$fit_metrics$R2_adjusted, 4), "\n")
      }
    }, error = function(e) {
      cat("  Error processing", state, ":", conditionMessage(e), "\n")
    })
  }
  
  return(list(
    coefficients = all_coefs,
    fit_metrics = all_metrics
  ))
}

# Define maximum lag to consider for each state
MAX_LAG <- 4  # Can be adjusted based on data availability

# Run the analysis with optimal lag selection
cat("\nRunning Bayesian VAR analysis with optimal lag selection (max lag =", MAX_LAG, ")...\n")
results <- all_states_var(df_standard, MAX_LAG)
all_coefficients <- results$coefficients
all_metrics <- results$fit_metrics

# Save results to CSV
write.csv(all_coefficients, "bayesian_var_coefficients.csv", row.names = FALSE)
write.csv(all_metrics, "bayesian_var_fit_metrics.csv", row.names = FALSE)

# Print a summary
cat("\nAnalysis complete.\n")
cat("Processed", length(unique(all_coefficients$State)), "states.\n")
cat("Coefficient results saved to 'bayesian_var_coefficients.csv'\n")
cat("Model fit metrics saved to 'bayesian_var_fit_metrics.csv'\n")

# Function to generate a summary by state
summarize_by_state <- function(all_coefs, all_metrics) {
  states <- unique(all_coefs$State)
  summary_list <- list()
  
  for (state in states) {
    state_coefs <- all_coefs[all_coefs$State == state, ]
    state_metrics <- all_metrics[all_metrics$State == state, ]
    
    # Summarize
    significant_count <- sum(state_coefs$Significant)
    total_count <- nrow(state_coefs)
    lag_order <- state_coefs$Lag_Order[1]
    
    # Calculate average coefficient magnitude by equation
    eq_summary <- state_coefs %>%
      group_by(Equation) %>%
      summarize(
        Avg_Coef_Magnitude = mean(abs(Coefficient)),
        Significant_Count = sum(Significant),
        Total_Coefs = n(),
        Min_P_Value = min(P_Value),
        Max_P_Value = max(P_Value)
      )
    
    # Add state and lag information
    eq_summary$State <- state
    eq_summary$Lag_Order <- lag_order
    eq_summary$R2_overall <- state_metrics$R2_overall
    eq_summary$R2_adjusted <- state_metrics$R2_adjusted
    
    # Add to list
    summary_list[[state]] <- eq_summary
  }
  
  # Combine results
  all_summary <- do.call(rbind, summary_list)
  
  return(all_summary)
}

# Generate and save summary
if (nrow(all_coefficients) > 0) {
  state_summary <- summarize_by_state(all_coefficients, all_metrics)
  write.csv(state_summary, "bayesian_var_state_summary.csv", row.names = FALSE)
  cat("State summary saved to 'bayesian_var_state_summary.csv'\n")
} else {
  cat("No coefficients were generated. Check the error messages above.\n")
}

# Print the best and worst models based on BIC
if (nrow(all_metrics) > 0) {
  # Sort by BIC (lower is better)
  sorted_metrics <- all_metrics[order(all_metrics$BIC), ]
  
  cat("\nTop 5 states with best model fit (lowest BIC):\n")
  print(head(sorted_metrics[, c("State", "Lag_Order", "BIC", "R2_overall", "R2_adjusted")], 5))
  
  cat("\nBottom 5 states with worst model fit (highest BIC):\n")
  print(tail(sorted_metrics[, c("State", "Lag_Order", "BIC", "R2_overall", "R2_adjusted")], 5))
  
  # Rank states by overall R-squared
  r2_ranked <- all_metrics[order(-all_metrics$R2_overall), ]
  cat("\nStates ranked by overall R-squared (highest to lowest):\n")
  print(r2_ranked[, c("State", "Lag_Order", "R2_overall", "R2_adjusted", "BIC")])
  
  # Save the ranking to CSV
  write.csv(r2_ranked, "bayesian_var_state_ranking.csv", row.names = FALSE)
  cat("State ranking saved to 'bayesian_var_state_ranking.csv'\n")
  
  # Show distribution of optimal lag orders
  lag_distribution <- table(all_metrics$Lag_Order)
  cat("\nDistribution of optimal lag orders across states:\n")
  print(lag_distribution)
}




# Function to extract vectorized VAR coefficients and IRFs
extract_var_features <- function(coef_file = "bayesian_var_coefficients.csv", horizon = 8) {
  # Load the coefficient data
  coef_data <- read.csv(coef_file)
  
  # Get unique states
  states <- unique(coef_data$State)
  
  # Variables in the system
  vars <- c("gdp_growth", "employment_growth", "EPU_value_Quarterly", 
            "composite_rates", "growth_STLFSI4")
  n_vars <- length(vars)
  
  # Create data frames to store results
  coef_vectors <- data.frame(State = character(), Delta_coef = I(list()), stringsAsFactors = FALSE)
  irf_vectors <- data.frame(State = character(), Delta_irf = I(list()), stringsAsFactors = FALSE)
  
  # Process each state
  for (state in states) {
    # Get data for this state
    state_coefs <- coef_data[coef_data$State == state, ]
    lag_order <- unique(state_coefs$Lag_Order)
    
    # Create coefficient matrices for each lag
    A_matrices <- list()
    
    for (lag in 1:lag_order) {
      # Initialize matrix for this lag
      A <- matrix(0, nrow = n_vars, ncol = n_vars)
      rownames(A) <- vars
      colnames(A) <- vars
      
      # Fill in coefficients
      for (i in 1:n_vars) {
        eq_name <- vars[i]  # Equation (dependent variable)
        
        for (j in 1:n_vars) {
          var_name <- vars[j]  # Independent variable
          lag_var <- paste0(var_name, "_lag", lag)
          
          # Find the coefficient
          coef_row <- state_coefs[state_coefs$Equation == eq_name & 
                                    state_coefs$Variable == lag_var, ]
          
          if (nrow(coef_row) > 0) {
            A[i, j] <- coef_row$Coefficient
          }
        }
      }
      
      A_matrices[[lag]] <- A
    }
    
    # Vectorize the coefficients (Δᵢᶜᵒᵉᶠ) - equation (4)
    coef_vec <- c()
    for (lag in 1:lag_order) {
      coef_vec <- c(coef_vec, as.vector(A_matrices[[lag]]))
    }
    
    # Calculate IRFs
    irf_matrices <- calculate_irfs(A_matrices, horizon)
    
    # Vectorize IRFs (Δᵢᴵᴿᶠ) - equation (5)
    irf_vec <- c()
    for (h in 1:(horizon+1)) {  # +1 because we include the impact period
      irf_vec <- c(irf_vec, as.vector(irf_matrices[[h]]))
    }
    
    # Add to data frames
    coef_vectors <- rbind(coef_vectors, 
                          data.frame(State = state, 
                                     Delta_coef = I(list(coef_vec)), 
                                     stringsAsFactors = FALSE))
    
    irf_vectors <- rbind(irf_vectors, 
                         data.frame(State = state, 
                                    Delta_irf = I(list(irf_vec)), 
                                    stringsAsFactors = FALSE))
  }
  
  # Return the results
  return(list(coef_vectors = coef_vectors, 
              irf_vectors = irf_vectors))
}

# Function to calculate IRFs
calculate_irfs <- function(A_matrices, horizon = 8) {
  lag_order <- length(A_matrices)
  n_vars <- nrow(A_matrices[[1]])
  
  # Create companion matrix
  companion <- matrix(0, nrow = n_vars * lag_order, ncol = n_vars * lag_order)
  
  # Fill with coefficient matrices
  for (lag in 1:lag_order) {
    companion[1:n_vars, ((lag-1)*n_vars+1):(lag*n_vars)] <- A_matrices[[lag]]
  }
  
  # Fill identity matrices in lower part
  if (lag_order > 1) {
    for (i in 1:(lag_order-1)) {
      companion[(i*n_vars+1):((i+1)*n_vars), ((i-1)*n_vars+1):(i*n_vars)] <- diag(n_vars)
    }
  }
  
  # Calculate IRFs
  irf_matrices <- list()
  irf_matrices[[1]] <- diag(n_vars)  # Initial impact
  
  # Recursive calculation for horizon periods
  for (h in 1:horizon) {
    impact <- matrix(0, nrow = n_vars * lag_order, ncol = n_vars)
    impact[1:n_vars, ] <- irf_matrices[[h]]
    
    next_impact <- companion %*% impact
    irf_matrices[[h+1]] <- next_impact[1:n_vars, ]
  }
  
  return(irf_matrices)
}

# Function to save the vectors as CSV files for easy use
save_vectors_as_csv <- function(features) {
  # For coefficient vectors
  coef_data <- features$coef_vectors
  
  # Create a matrix with states as rows and coefficients as columns
  coef_mat <- matrix(NA, nrow = nrow(coef_data), ncol = length(coef_data$Delta_coef[[1]]))
  rownames(coef_mat) <- coef_data$State
  
  for (i in 1:nrow(coef_data)) {
    coef_mat[i, ] <- unlist(coef_data$Delta_coef[i])
  }
  
  # Convert to data frame and save
  coef_df <- as.data.frame(coef_mat)
  colnames(coef_df) <- paste0("coef_", 1:ncol(coef_df))
  coef_df$State <- rownames(coef_df)
  coef_df <- coef_df[, c("State", paste0("coef_", 1:ncol(coef_mat)))]
  
  write.csv(coef_df, "bayesian_var_vectorized_coefficients.csv", row.names = FALSE)
  
  # For IRF vectors
  irf_data <- features$irf_vectors
  
  # Create a matrix with states as rows and IRFs as columns
  irf_mat <- matrix(NA, nrow = nrow(irf_data), ncol = length(irf_data$Delta_irf[[1]]))
  rownames(irf_mat) <- irf_data$State
  
  for (i in 1:nrow(irf_data)) {
    irf_mat[i, ] <- unlist(irf_data$Delta_irf[i])
  }
  
  # Convert to data frame and save
  irf_df <- as.data.frame(irf_mat)
  colnames(irf_df) <- paste0("irf_", 1:ncol(irf_df))
  irf_df$State <- rownames(irf_df)
  irf_df <- irf_df[, c("State", paste0("irf_", 1:ncol(irf_mat)))]
  
  write.csv(irf_df, "bayesian_var_vectorized_irfs.csv", row.names = FALSE)
  
  cat("Saved vectorized coefficients to 'bayesian_var_vectorized_coefficients.csv'\n")
  cat("Saved vectorized IRFs to 'bayesian_var_vectorized_irfs.csv'\n")
}

# Main execution
# Extract features with a 10-period horizon for IRFs
features <- extract_var_features(coef_file = "bayesian_var_coefficients.csv", horizon = 8)

# Save as RDS for R use
saveRDS(features, "bayesian_var_features.rds")
cat("Saved features to 'bayesian_var_features.rds'\n")

# Also save as CSV files for easy use in other software
save_vectors_as_csv(features)

# Print dimensions of vectors for verification
cat("\nDimensions of vectorized features:\n")
cat("Coefficient vectors:\n")
for (i in 1:nrow(features$coef_vectors)) {
  state <- features$coef_vectors$State[i]
  coef_length <- length(features$coef_vectors$Delta_coef[[i]])
  cat(sprintf("  %s: %d elements\n", state, coef_length))
}

cat("\nIRF vectors:\n")
for (i in 1:nrow(features$irf_vectors)) {
  state <- features$irf_vectors$State[i]
  irf_length <- length(features$irf_vectors$Delta_irf[[i]])
  cat(sprintf("  %s: %d elements\n", state, irf_length))
}
