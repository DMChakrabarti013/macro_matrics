#!/usr/bin/env Rscript

# =============================================================================
# Full Synthetic Control + Block Bootstrap in R
# =============================================================================

options(error = recover)

# 0) Dependencies -------------------------------------------------------------
library(readr)    # read_csv
library(dplyr)    # data wrangling
library(tidyr)    # pivot_wider
library(purrr)    # map_dfr
library(CVXR)     # convex optimization
library(glmnet)   # penalized regression
library(caret)    # createFolds
library(quadprog) # fallback QP
library(zoo)      # as.yearqtr
library(ggplot2)  # plotting

options(error = recover)

# 1) Read CSVs (suppress col specs) ------------------------------------------
baseline_df       <- read_csv("/Users/dwaipayanchakrabarti/Downloads/baseline_data_synthetic.csv", show_col_types = FALSE)
coef_df           <- read_csv("/Users/dwaipayanchakrabarti/Downloads/bayesian_var_vectorized_coefficients.csv", show_col_types = FALSE)[,1:26]
irf_df            <- read_csv("/Users/dwaipayanchakrabarti/Downloads/selected_irfs.csv", show_col_types = FALSE)
dt_df             <- read_csv("/Users/dwaipayanchakrabarti/Downloads/decision_tree_scores.csv", show_col_types = FALSE)
gpr_df            <- read_csv("/Users/dwaipayanchakrabarti/Downloads/gpr_state_scores.csv", show_col_types = FALSE)
rf_df             <- read_csv("/Users/dwaipayanchakrabarti/Downloads/rf_proximity_scores.csv", show_col_types = FALSE)
xgb_df            <- read_csv("/Users/dwaipayanchakrabarti/Downloads/xgboost_scores.csv", show_col_types = FALSE)
gdp_panel         <- read_csv("/Users/dwaipayanchakrabarti/Downloads/gdp_growth_data.csv", show_col_types = FALSE)

# 2) Pivot GDP to state × quarter ---------------------------------------------
gdp_wide <-
  gdp_panel %>%
  mutate(YearQuarter = gsub("\\s+", "", YearQuarter)) %>%  # "2005 Q2" → "2005Q2"
  group_by(State, YearQuarter) %>%
  summarise(gdp_growth = mean(gdp_growth, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = State, names_from = YearQuarter, values_from = gdp_growth) %>%
  arrange(State) %>%
  as.data.frame()

rownames(gdp_wide) <- gdp_wide$State
gdp_wide$State    <- NULL

# 3) Aggregate covariates ----------------------------------------------------
baseline_cov <-
  baseline_df %>%
  group_by(State) %>%
  summarise(
    assetquality_diff  = mean(assetquality_diff,  na.rm = TRUE),
    low_diff           = mean(low_diff,           na.rm = TRUE),
    profitability_diff = mean(profitability_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()

# 4) Build ML‐scores table ---------------------------------------------------
ml_scores <-
  tibble(State = irf_df$State) %>%
  mutate(
    decision_tree_score = dt_df$G_score_tree[ match(State, dt_df$State) ],
    rf_score            = rf_df$G_score[       match(State, rf_df$State) ],
    gpr_score           = gpr_df$G_score[      match(State, gpr_df$State) ],
    xgb_score           = xgb_df$G_score[      match(State, xgb_df$State) ]
  ) %>%
  as.data.frame()

# 5) Filter to common states -------------------------------------------------
common_states <- Reduce(intersect, list(
  irf_df$State,
  coef_df$State,
  baseline_cov$State,
  ml_scores$State,
  rownames(gdp_wide)
))

irf_df       <- irf_df[ irf_df$State %in% common_states, ]
coef_df      <- coef_df[ coef_df$State %in% common_states, ]
baseline_cov <- baseline_cov[ baseline_cov$State %in% common_states, ]
ml_scores    <- ml_scores[ ml_scores$State %in% common_states, ]
gdp_wide     <- gdp_wide[common_states, , drop = FALSE]

# 6) Build numeric matrices --------------------------------------------------
irf_mat      <- scale( data.matrix(irf_df[match(common_states, irf_df$State), -1]) )
coef_mat     <- scale( data.matrix(coef_df[match(common_states, coef_df$State), -1]) )
baseline_mat <- data.matrix(baseline_cov[ match(common_states, baseline_cov$State), -1 ])

# ML scores → scaled numeric matrix
dfm       <- ml_scores %>% distinct(State, .keep_all = TRUE) %>% as.data.frame()
rownames(dfm) <- dfm$State; dfm$State <- NULL
dfm[]     <- lapply(dfm, function(col) as.numeric(as.character(col)))
ml_mat    <- scale( data.matrix(dfm[common_states, , drop = FALSE]) )

# 7) Pre-treatment GDP (2005Q2–2008Q2) ---------------------------------------
pre_cols <- which(colnames(gdp_wide) >= "2005Q2" & colnames(gdp_wide) <= "2008Q2")
y_pre     <- gdp_wide[, pre_cols, drop = FALSE]

# =============================================================================
# 8) SyntheticControlMatcher definition  (with corrected fit_predict)
# =============================================================================
SyntheticControlMatcher <- function(alpha = 0.5,
                                    lambda_grid = 10^seq(-1, 2, length.out = 10),
                                    gamma = 0.5) {
  state <- list(alpha = alpha, lambda_grid = lambda_grid, gamma = gamma)
  
  prepare_data <- function(irfs, coeffs, baseline_data, ml_scores) {
    state$states       <- irfs$State
    state$irf_mat      <- scale( data.matrix(irfs  %>% select(-State)) )
    state$coef_mat     <- scale( data.matrix(coeffs %>% select(-State)) )
    
    # IRF relevance weights
    m   <- min(ncol(state$irf_mat), ncol(state$coef_mat))
    cm  <- tryCatch(cor(state$irf_mat[,1:m], state$coef_mat[,1:m]), error = function(e) NULL)
    w   <- if (is.null(cm)) rep(1, m) else apply(abs(cm), 1, max)
    w[is.na(w)] <- 1
    state$irf_weights <- w
    
    # Baseline covariates → numeric matrix
    dfb <- baseline_data %>%
      select(State, assetquality_diff, low_diff, profitability_diff) %>%
      as.data.frame()
    rownames(dfb) <- dfb$State 
    dfb$State <- NULL
    state$baseline_static <- data.matrix(dfb)
    
    # ML scores → numeric matrix
    dfm <- ml_scores %>% distinct(State, .keep_all = TRUE) %>% as.data.frame()
    rownames(dfm) <- dfm$State; dfm$State <- NULL
    dfm[] <- lapply(dfm, function(col) as.numeric(as.character(col)))
    state$ml_scores <- data.matrix(dfm[state$states, , drop = FALSE])
    
    return(state)
  }
  
  optimize_ml_weights <- function(y_pre) {
    y_mean <- rowMeans(y_pre[state$states, ], na.rm = TRUE)
    X      <- scale(state$ml_scores)
    cv     <- tryCatch(cv.glmnet(X, y_mean, alpha = 0, lower.limits = 0), error = function(e) NULL)
    w      <- if (is.null(cv)) {
      rep(1, ncol(X))
    } else {
      co <- as.vector(coef(cv, s = "lambda.min"))[-1]
      if (sum(co) <= 0) rep(1, length(co)) else co
    }
    w <- w / sum(w)
    names(w) <- colnames(state$ml_scores)
    state$ml_weights <- w
    return(w)
  }
  
  combine_delta <- function() {
    m <- min(ncol(state$irf_mat), ncol(state$coef_mat))
    I <- state$irf_mat[,1:m]
    C <- state$coef_mat[,1:m]
    sweep(state$alpha * I + (1 - state$alpha) * C, 2, state$irf_weights[1:m], "*")
  }
  
  compute_cv_error <- function(G) {
    if (length(unique(G)) < 2) return(0)
    folds <- createFolds(G, k = min(5, length(G)-1))
    errs  <- c()
    for(tr in folds) {
      te <- setdiff(seq_along(G), tr)
      if (length(unique(G[tr])) < 2) next
      fit   <- glm(G[tr] ~ ., data = as.data.frame(state$ml_scores)[tr, ], family = "binomial")
      pred  <- predict(fit,  as.data.frame(state$ml_scores)[te, ], type = "response")
      errs  <- c(errs, mean((pred > 0.5) != G[te]))
    }
    if (length(errs) == 0) 0 else mean(errs)
  }
  
  synth_weights <- function(d_i, d_j, o_i, o_j, lam) {
    J  <- nrow(d_j)
    wv <- CVXR::Variable(J)
    
    obj <- CVXR::sum_squares(d_i - t(d_j) %*% wv) +
      lam * CVXR::sum_squares(o_i - t(o_j) %*% wv)
    
    pr  <- CVXR::Problem(CVXR::Minimize(obj),
                         list(wv >= 0, sum(wv) == 1))
    sol <- tryCatch(CVXR::solve(pr), error = function(e) NULL)
    if (!is.null(sol) &&
        sol$status %in% c("optimal","optimal_inaccurate")) {
      return(as.numeric(sol$getValue(wv)))
    }
    
    # fallback QP:
    D <- t(d_j) %*% d_j + lam * (t(o_j) %*% o_j)
    d <- as.vector(t(d_j) %*% d_i + lam * (t(o_j) %*% o_i))
    ev <- eigen(D, symmetric=TRUE, only.values=TRUE)$values
    if (min(ev) <= 0) D <- D + diag(nrow(D))*(abs(min(ev))+1e-6)
    
    solve.QP::solve.QP(D, d,
                       cbind(rep(1,J), diag(J)),
                       c(1, rep(0,J)),
                       meq = 1)$solution
  }
  
  optimize_lambda <- function(combined, y_pre, G) {
    cv_err <- compute_cv_error(G)
    best   <- list(loss = Inf, lam = state$lambda_grid[1])
    for (lam in state$lambda_grid) {
      errs <- c()
      for (i in which(G == 1)) {
        donors <- which(G == 0)
        d_i    <- combined[i, ]
        d_j    <- combined[donors, , drop = FALSE]
        o_i    <- state$baseline_static[i, ]
        o_j    <- state$baseline_static[donors, , drop = FALSE]
        w      <- synth_weights(d_i, d_j, o_i, o_j, lam)
        if (length(w) != length(donors)) next
        y_i <- as.numeric(y_pre[i, ])
        y_d <- as.matrix(y_pre[donors, , drop = FALSE])
        errs <- c(errs, mean((y_i - as.numeric(t(w) %*% y_d))^2, na.rm = TRUE))
      }
      pre_err <- if (length(errs) == 0) Inf else mean(errs)
      loss    <- state$gamma * pre_err + (1 - state$gamma) * cv_err
      if (loss < best$loss) {
        best$loss <- loss
        best$lam  <- lam
      }
    }
    return(best$lam)
  }
  
  fit_predict <- function(irfs, coeffs, baseline_data, ml_scores, y_pre, y_full, T0) {
    state <<- prepare_data(irfs, coeffs, baseline_data, ml_scores)
    ml_w   <- optimize_ml_weights(y_pre)
    state$ml_weights <- ml_w
    
    delta <- combine_delta()
    
    # composite ML score → grouping
    comp <- as.vector(scale(state$ml_scores) %*% ml_w)
    state$c_threshold <- median(comp, na.rm = TRUE)
    G <- as.integer(comp > state$c_threshold)
    
    lam <- optimize_lambda(delta, y_pre, G)
    
    cf <- list(); te <- list()
    for (i in which(G == 1)) {
      donors <- which(G == 0)
      
      # observed & synthetic series
      y_o      <- as.numeric(y_full[i, ])
      names(y_o) <- colnames(y_full)
      Y_d_mat  <- as.matrix(y_full[donors, , drop = FALSE])
      y_s_vals <- as.numeric(t( synth_weights(delta[i,], delta[donors,], 
                                              state$baseline_static[i,], 
                                              state$baseline_static[donors,], 
                                              lam) ) %*% Y_d_mat)
      names(y_s_vals) <- colnames(y_full)
      
      # full effect & split
      eff_full     <- y_o - y_s_vals
      t0_idx       <- which(colnames(y_full) == T0)
      if (length(t0_idx) == 0) t0_idx <- floor(length(eff_full)/2)
      active_eff     <- eff_full[ t0_idx:length(eff_full) ]
      cumulative_eff <- cumsum(active_eff)
      
      cf[[ common_states[i] ]] <- list(
        observed  = setNames(y_o,      colnames(y_full)),
        synthetic = setNames(y_s_vals, colnames(y_full)),
        weights   = setNames(
          synth_weights(delta[i,], delta[donors,], 
                        state$baseline_static[i,], 
                        state$baseline_static[donors,], lam),
          common_states[donors]
        )
      )
      te[[ common_states[i] ]] <- list(
        active     = setNames(active_eff,     colnames(y_full)[t0_idx:length(eff_full)]),
        cumulative = setNames(cumulative_eff, colnames(y_full)[t0_idx:length(eff_full)])
      )
    }
    
    list(
      alpha            = state$alpha,
      lambda           = lam,
      gamma            = state$gamma,
      c_threshold      = state$c_threshold,
      ml_weights       = state$ml_weights,
      groups           = G,
      counterfactuals  = cf,
      treatment_effects= te,
      treated_states   = names(cf),
      donor_states     = common_states[G == 0]
    )
  }
  
  list(fit_predict = fit_predict)
}

# 9) Run the SCM --------------------------------------------------------------
alpha           <- 0.5
gamma           <- 0.5
lambda_grid     <- 10^seq(-1, 2, length.out = 10)
treatment_start <- "2008Q3"

scm     <- SyntheticControlMatcher(alpha, lambda_grid, gamma)
results <- scm$fit_predict(irf_df, coef_df, baseline_cov, ml_scores,
                           y_pre, gdp_wide, treatment_start)

message("Optimal λ = ", results$lambda)
message("Treated states: ", paste(results$treated_states, collapse = ", "))
message("Donor states:   ", paste(results$donor_states, collapse = ", "))

# 10) stationary‐bootstrap 95% CIs on *active* effects -----------------------------
# … assume you’ve already run your SCM code through `results <- scm$fit_predict(...)`
# and have `gdp_wide`, `pre_cols`, `irf_df`, `coef_df`, `baseline_cov`, `ml_scores`,
# `treatment_start`, and the `scm` object in your workspace.

library(dplyr)
library(purrr)
library(zoo)
library(ggplot2)

# 1) Pre‐treatment dimensions & block parameter -------------------------------
pre_cols     <- which(colnames(gdp_wide) >= "2005Q2" &
                        colnames(gdp_wide) <= "2008Q2")
T_pre        <- length(pre_cols)
block_length <- 2   # smaller blocks → more variation (~T_pre^(1/3))

# 2) One stationary‐bootstrap replicate ---------------------------------------
# ------------- full‐series block bootstrap helper -------------
# ------------- full‐series block bootstrap helper -------------
bootstrap_one <- function(b) {
  if (b %% 10 == 0) cat("Bootstrap replicate:", b, "\n")
  
  # 1) draw blocks across the **entire** series
  full_T   <- ncol(gdp_wide)               # e.g. quarters 2005Q2–2020Q1
  starts   <- sample(1:(full_T - block_length + 1),
                     ceiling(full_T / block_length),
                     replace = TRUE)
  full_idx <- unlist(lapply(starts, function(s) s:(s + block_length - 1)))[1:full_T]
  
  # 2) assemble bootstrapped full series & rename quarters back to original
  y_full_boot <- gdp_wide[, full_idx, drop = FALSE]
  colnames(y_full_boot) <- colnames(gdp_wide)  # keep the Q labels
  
  # 3) extract the pre‐treatment portion by your original pre_cols
  y_pre_boot <- y_full_boot[, pre_cols, drop = FALSE]
  
  # 4) re-run SCM on these bootstrapped series
  scb <- tryCatch(
    scm$fit_predict(
      irf_df, coef_df, baseline_cov, ml_scores,
      y_pre_boot, y_full_boot, treatment_start
    ),
    error = function(e) NULL
  )
  if (is.null(scb) || length(scb$treatment_effects)==0) return(NULL)
  
  # 5) pull out the active effects
  map_dfr(names(scb$treatment_effects), function(st) {
    te <- scb$treatment_effects[[st]]$active
    tibble(
      State     = st,
      Quarter   = names(te),
      Effect    = as.numeric(te),
      replicate = b
    )
  })
}

# 3) Run B replicates ----------------------------------------------------------
set.seed(123)
B <- 500
block_length <- 2    # in quarters, say

boot_df <- map_dfr(1:B, bootstrap_one) %>% drop_na()

ci_df <- boot_df %>%
  group_by(State, Quarter) %>%
  summarize(
    lower = quantile(Effect, 0.025, na.rm = TRUE),
    upper = quantile(Effect, 0.975, na.rm = TRUE),
    mean  = mean(Effect,    na.rm = TRUE),
    .groups = "drop"
  )

# 5) Merge with your original point‐estimate effects --------------------------
obs_df <- map_dfr(names(results$treatment_effects), function(st) {
  te <- results$treatment_effects[[st]]$active
  tibble(State = st, Quarter = names(te), point = as.numeric(te))
})

plot_df <- obs_df %>%
  left_join(ci_df, by = c("State", "Quarter")) %>%
  mutate(Date = as.yearqtr(gsub("Q"," ", Quarter), "%Y %q"))

# 6) Plot all states, now without warnings ------------------------------------
vline <- as.yearqtr("2008 Q3", "%Y %q")
vline_df <- data.frame(x = as.numeric(vline))

ggplot(plot_df, aes(x = Date, y = point)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80") +
  geom_line(color = "blue") +
  facet_wrap(~ State, ncol = 5) +
  geom_vline(data = vline_df, aes(xintercept = x),
             linetype = "dashed", inherit.aes = FALSE) +
  scale_x_yearqtr(format = "%Y Q%q") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    strip.text  = element_text(size = 7)
  ) +
  labs(
    x     = "Quarter",
    y     = "Treatment Effect (GDP Growth)",
    title = "SCM Active Treatment Effects with 95% Stationary Bootstrap CIs"
  )
