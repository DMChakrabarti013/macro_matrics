### stattionary 95% bootstrap

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
bootstrap_one <- function(b) {
  if (b %% 10 == 0) cat("Bootstrap replicate:", b, "\n")
  
  # 1) draw blocks across the **entire** series
  full_T   <- ncol(gdp_wide)               
  starts   <- sample(1:(full_T - block_length + 1),
                     ceiling(full_T / block_length),
                     replace = TRUE)
  full_idx <- unlist(lapply(starts, function(s) s:(s + block_length - 1)))[1:full_T]
  
 
  y_full_boot <- gdp_wide[, full_idx, drop = FALSE]
  colnames(y_full_boot) <- colnames(gdp_wide)  
  
 
  y_pre_boot <- y_full_boot[, pre_cols, drop = FALSE]
  
  # re-run SCM 
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

# 5) Merge with point‐estimate effects --------------------------
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
