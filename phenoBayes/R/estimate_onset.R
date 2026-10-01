#' Estimate Phenological Onset Shift Volume
#' 
#' @param model A fitted mgcv GAM/BAM object
#' @param newdata A data.frame with prediction grid (365 days) for the shifted scenario
#' @param ref_newdata A data.frame with prediction grid for the reference scenario
#' @param n_sims Integer, number of posterior draws
#' @param sim_mat Optional pre-computed simulation matrix for newdata
#' @param ref_sim_mat Optional pre-computed simulation matrix for ref_newdata
#' @export
estimate_onset <- function(model, newdata, ref_newdata, n_sims = 1000, sim_mat = NULL, ref_sim_mat = NULL) {
  if (is.null(sim_mat)) sim_mat <- simulate_posterior_density(model, newdata, n_sims)
  if (is.null(ref_sim_mat)) ref_sim_mat <- simulate_posterior_density(model, ref_newdata, n_sims)
  
  res <- numeric(ncol(sim_mat))
  for(i in seq_len(ncol(sim_mat))) {
    td <- sim_mat[, i]
    rd <- ref_sim_mat[, i]
    
    # Positive means td (shifted) emerged earlier than rd (reference)
    cdiff <- cumsum(td) - cumsum(rd)
    pr <- which.max(rd)
    
    # Find the day of maximum absolute divergence before the reference peak
    max_idx <- which.max(abs(cdiff[1:pr]))
    res[i] <- cdiff[max_idx]
  }
  return(list(Mean = mean(res), Lower_CI = quantile(res, 0.025), Upper_CI = quantile(res, 0.975), Posterior_Draws = res))
}
