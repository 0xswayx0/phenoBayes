#' Estimate Phenological Centroid (DOBY)
#' 
#' @param model A fitted mgcv GAM/BAM object
#' @param newdata A data.frame with prediction grid (365 days)
#' @param n_sims Integer, number of posterior draws
#' @param sim_mat Optional pre-computed simulation matrix
#' @export
estimate_centroid <- function(model, newdata, n_sims = 1000, sim_mat = NULL) {
  if (is.null(sim_mat)) sim_mat <- simulate_posterior_density(model, newdata, n_sims)
  
  # Base mean profile for DOBY centering
  base_X <- predict(model, newdata = newdata, type = "lpmatrix")
  re_cols <- grep("s\\(fYear\\)|s\\(SubSite\\)", colnames(base_X))
  if (length(re_cols) > 0) base_X[, re_cols] <- 0
  base_dens <- apply(exp(base_X %*% coef(model)), 1, mean)
  
  # Paper 1 robust D_min (center of minimum period using circular mean)
  min_val <- min(base_dens)
  min_days <- which(base_dens <= min_val + 1e-9)
  theta <- min_days * 2 * pi / 365
  mean_theta <- atan2(mean(sin(theta)), mean(cos(theta)))
  mean_theta <- ifelse(mean_theta < 0, mean_theta + 2*pi, mean_theta)
  D_min <- round(mean_theta * 365 / (2 * pi))
  if (D_min == 0) D_min <- 365
  
  res <- numeric(ncol(sim_mat))
  for (i in seq_len(ncol(sim_mat))) {
    td <- sim_mat[, i]
    shifted <- td[order(( (1:365 - D_min) %% 365 ) + 1)]
    doby_c <- sum((1:365) * shifted)
    res[i] <- ( (doby_c + D_min - 2) %% 365 ) + 1
  }
  
  return(list(
    Mean = mean(res),
    Lower_CI = quantile(res, 0.025),
    Upper_CI = quantile(res, 0.975),
    Posterior_Draws = res,
    D_min = D_min
  ))
}
