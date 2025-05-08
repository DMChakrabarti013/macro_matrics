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
baseline_df <- baseline_df %>%
  mutate(
    YQ = as.yearqtr(YearQuarter, format = "%Y Q%q")
  )

# then filter to pre‐treatment
baseline_pre <- baseline_df %>%
  filter(YQ < as.yearqtr("2008 Q3", "%Y Q%q"))

baseline_cov <-
  baseline_pre %>%
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

######## results ##########
library(zoo)

obs <- results$counterfactuals$maryland$observed
syn <- results$counterfactuals$maryland$synthetic
qtr <- names(obs)

# convert "2005Q2" → yearqtr
obs <- results$counterfactuals$maryland$observed
syn <- results$counterfactuals$maryland$synthetic
qtr <- names(obs)

# numeric indices for x
ix <- seq_along(obs)

plot(ix, obs, type="l", col="black", pch=16,
     xaxt="n",                    # suppress default x‐axis
     xlab="Quarter", ylab="GDP Growth",
     main="Maryland: Observed vs Synthetic")
lines(ix, syn, type="l", col="blue", pch=17)

# add our own x‐axis with labels
axis(1, at=ix, labels=qtr, las=2, cex.axis=0.7)  

# vertical line at treatment
t0 <- which(qtr=="2008Q3")
t1 <- which(qtr=='2017Q1')
abline(v=t0, lty=2)
abline(v=t1, lty=2)

legend("topleft", legend=c("Observed","Synthetic"),
       col=c("black","blue"), pch=c(16,17), bty="n")

######### all results ######
library(ggplot2)
library(zoo)
library(reshape2)

# 1) Pull everything into a long data.frame
all_df <- do.call(rbind, lapply(names(results$counterfactuals), function(st) {
  obs <- results$counterfactuals[[st]]$observed
  syn <- results$counterfactuals[[st]]$synthetic
  qtr <- names(obs)
  data.frame(
    State     = st,
    Quarter   = qtr,
    Observed  = as.numeric(obs),
    Synthetic = as.numeric(syn),
    stringsAsFactors = FALSE
  )
}))

# 2) Melt so that Observed/Synthetic become a single “Type” column
all_long <- melt(all_df,
                 id.vars   = c("State","Quarter"),
                 measure.vars = c("Observed","Synthetic"),
                 variable.name = "Type",
                 value.name    = "Value")

# 3) Convert “2005Q2” strings into true yearqtr
all_long$Date <- as.yearqtr(gsub("Q", " ", all_long$Quarter), format = "%Y %q")

# 4) Plot 
end_qtr <- as.yearqtr("2017 Q2", "%Y Q%q")

ggplot(all_long, aes(x = Date, y = Value, color = Type)) +
  geom_line(size = 0.7) +
  facet_wrap(~ State, ncol = 5) +
  geom_vline(xintercept = as.yearqtr("2008 Q3", "%Y %q"),
             linetype = "dashed", alpha = 0.5) +
  scale_x_yearqtr(format = "%Y Q%q",
                  limits = c(min(all_long$Date), end_qtr)) +
  theme_minimal() +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1, size = 6),
    strip.text     = element_text(size = 8),
    legend.position = "bottom"
  ) +
  labs(
    x     = "Quarter",
    y     = "GDP Growth",
    color = NULL,
    title = "Observed vs Synthetic GDP Growth by State"
  )
