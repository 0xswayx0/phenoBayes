#' Estimate Phenological Offset Shift Volume
#' 
#' @param model A fitted mgcv GAM/BAM object
#' @param newdata A data.frame with prediction grid (365 days) for the shifted scenario
#' @param ref_newdata A data.frame with prediction grid for the reference scenario
#' @param n_sims Integer, number of posterior draws
#' @param sim_mat Optional pre-computed simulation matrix for newdata
#' @param ref_sim_mat Optional pre-computed simulation matrix for ref_newdata
#' @export
estimate_offset <- function(model, newdata, ref_newdata, n_sims = 1000, sim_mat = NULL, ref_sim_mat = NULL) {
  if (is.null(sim_mat)) sim_mat <- simulate_posterior_density(model, newdata, n_sims)
  if (is.null(ref_sim_mat)) ref_sim_mat <- simulate_posterior_density(model, ref_newdata, n_sims)
  
  res <- numeric(ncol(sim_mat))
  for(i in seq_len(ncol(sim_mat))) {
    td <- sim_mat[, i]
    rd <- ref_sim_mat[, i]
    
    # Positive means td (shifted) persisted later than rd (reference)
    cdiff <- cumsum(rd) - cumsum(td)
    pr <- which.max(rd)
    
    # Find the day of maximum absolute divergence after the reference peak
    max_idx <- pr - 1 + which.max(abs(cdiff[pr:365]))
    res[i] <- cdiff[max_idx]
  }
  return(list(Mean = mean(res), Lower_CI = quantile(res, 0.025), Upper_CI = quantile(res, 0.975), Posterior_Draws = res))
}
