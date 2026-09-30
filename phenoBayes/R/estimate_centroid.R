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
  D_min <- which.min(base_dens)
  
  res <- numeric(ncol(sim_mat))
  for (i in seq_len(ncol(sim_mat))) {
    td <- sim_mat[, i]
    shifted <- td[order(( (1:365 - D_min) %% 365 ) + 1)]
    doby_c <- sum((1:365) * shifted)
    res[i] <- (doby_c + D_min - 1) %% 365
  }
  
  return(list(
    Mean = mean(res),
    Lower_CI = quantile(res, 0.025),
    Upper_CI = quantile(res, 0.975),
    Posterior_Draws = res
  ))
}
