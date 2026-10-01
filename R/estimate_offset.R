#' Estimate Phenological Offset Shift Volume
#' 
#' @param model A fitted mgcv GAM/BAM object
#' @param newdata A data.frame with prediction grid (365 days) for the shifted scenario
#' @param ref_newdata A data.frame with prediction grid for the reference scenario
#' @param n_sims Integer, number of posterior draws
#' @param sim_mat Optional pre-computed simulation matrix for newdata
#' @param ref_sim_mat Optional pre-computed simulation matrix for ref_newdata
#' @param threshold_fraction Density threshold to filter winter noise (default 0.05)
#' @export
estimate_offset <- function(model, newdata, ref_newdata, n_sims = 1000, sim_mat = NULL, ref_sim_mat = NULL, threshold_fraction = 0.05) {
  if (is.null(sim_mat)) sim_mat <- simulate_posterior_density(model, newdata, n_sims)
  if (is.null(ref_sim_mat)) ref_sim_mat <- simulate_posterior_density(model, ref_newdata, n_sims)
  
  base_X <- predict(model, newdata = ref_newdata, type = "lpmatrix")
  re_cols <- grep("s\\(fYear\\)|s\\(SubSite\\)", colnames(base_X))
  if (length(re_cols) > 0) base_X[, re_cols] <- 0
  base_dens <- apply(exp(base_X %*% coef(model)), 1, mean)
  min_val <- min(base_dens)
  min_days <- which(base_dens <= min_val + 1e-9)
  theta <- min_days * 2 * pi / 365
  mean_theta <- atan2(mean(sin(theta)), mean(cos(theta)))
  mean_theta <- ifelse(mean_theta < 0, mean_theta + 2*pi, mean_theta)
  D_min <- round(mean_theta * 365 / (2 * pi))
  if (D_min == 0) D_min <- 365
  
  res <- numeric(ncol(sim_mat))
  for(i in seq_len(ncol(sim_mat))) {
    td <- sim_mat[, i]
    rd <- ref_sim_mat[, i]
    
    td_shift <- td[order(( (1:365 - D_min) %% 365 ) + 1)]
    rd_shift <- rd[order(( (1:365 - D_min) %% 365 ) + 1)]
    
    # 1. Filter out winter noise
    active_mask <- (td_shift > threshold_fraction * max(td_shift)) | (rd_shift > threshold_fraction * max(rd_shift))
    
    # 2. Find crossings
    diff_sign <- sign(td_shift - rd_shift)
    diff_sign[diff_sign == 0] <- 1 
    crossings <- which(diff(diff_sign) != 0)
    
    # 3. Last valid crossing
    valid_crossings <- crossings[active_mask[crossings] | active_mask[crossings + 1]]
    if (length(valid_crossings) == 0) {
      last_cross <- which.max(abs(cumsum(td_shift) - cumsum(rd_shift)))
    } else {
      last_cross <- valid_crossings[length(valid_crossings)]
    }
    
    # Positive means shifted curve persisted later than reference curve
    cdiff <- cumsum(rd_shift) - cumsum(td_shift)
    res[i] <- cdiff[last_cross]
  }
  return(list(Mean = mean(res), Lower_CI = quantile(res, 0.025), Upper_CI = quantile(res, 0.975), Posterior_Draws = res))
}
