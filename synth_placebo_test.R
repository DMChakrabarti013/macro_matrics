###### placebo test #######

library(dplyr)
library(purrr)
library(zoo)


irf_mat  <- scale( data.matrix( dplyr::select(irf_df,  -State) ) )
coef_mat <- scale( data.matrix( dplyr::select(coef_df, -State) ) )

m        <- min(ncol(irf_mat), ncol(coef_mat))
cor_mat  <- cor(irf_mat[,1:m], coef_mat[,1:m])
irf_weights <- apply(abs(cor_mat), 1, max)
irf_weights[is.na(irf_weights)] <- 1

delta_composite <- sweep(
  alpha * irf_mat[,1:m] + (1 - alpha) * coef_mat[,1:m],
  2, irf_weights[1:m], `*`
)


baseline_static <- data.matrix(
  baseline_cov[, c("assetquality_diff","low_diff","profitability_diff")]
)

rownames(delta_composite) <- common_states
rownames(baseline_static) <- common_states

donors <- match(results$donor_states, common_states)

delta_i      <- delta_composite[i, , drop = FALSE]
delta_donors <- delta_composite[donors, , drop = FALSE]
o_i          <- baseline_static[i, , drop = FALSE]
o_donors     <- baseline_static[donors, , drop = FALSE]

library(CVXR)

solve_w <- function(j, donor_idx){
  m <- ncol(delta_composite)
  p <- ncol(baseline_static)
  
  d_j   <- matrix(delta_composite[j,   ], nrow = m, ncol = 1)
  D_d   <- delta_composite[donor_idx, , drop = FALSE]  # K × m
  O_d   <- baseline_static[donor_idx, , drop = FALSE]  # K × p
  o_j   <- matrix(baseline_static[j,   ], nrow = p, ncol = 1)
  
  w     <- Variable(length(donor_idx))
  
  obj   <- sum_squares(d_j - t(D_d) %*% w) +
    results$lambda * sum_squares(o_j - t(O_d) %*% w)
  
  prob  <- Problem(Minimize(obj),
                   list(w >= 0, sum(w) == 1))
  sol   <- solve(prob, quiet = TRUE)
  as.numeric(sol$getValue(w))
}


t0        <- which(colnames(gdp_wide) == treatment_start)
t_end   <- which(colnames(gdp_wide) == "2017Q2")
active_ix <- t0:t_end

tau_perm <- map_dbl(donors, function(j){
  pool_j <- setdiff(donors, j)
  
  w_j    <- solve_w(j, pool_j)
  
  y_j    <- as.numeric(gdp_wide[j,     ])
  Y_pool <- as.matrix(   gdp_wide[pool_j, ])
  y_syn  <- as.numeric(t(w_j) %*% Y_pool)
  
  mean((y_j - y_syn)[active_ix])
})

i <- which(common_states == results$treated_states[15]) # change state number here to get different permutation test

w_i     <- solve_w(i, donors)
y_i     <- as.numeric(gdp_wide[i,])
y_syn_i <- as.numeric(t(w_i) %*% as.matrix(gdp_wide[donors,]))
tau_obs <- mean((y_i - y_syn_i)[active_ix])

p_val <- mean(abs(tau_perm) >= abs(tau_obs))
p_val

par(mfrow = c(1,1))

par(mar = c(5, 4, 4, 2) + 0.1)

hist(tau_perm, breaks = 30, main = "Permutation distribution of Treatment Effect",
     xlab = expression(tau[j]))
abline(v = tau_obs, col = "red", lwd = 2)
title(sub = sprintf("Obs. Treatment Effect = %.3f, p = %.3f", tau_obs, p_val))
