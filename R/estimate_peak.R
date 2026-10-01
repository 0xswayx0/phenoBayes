#' Estimate Phenological Peak (Mode)
#' 
#' @param model A fitted mgcv GAM/BAM object
#' @param newdata A data.frame with prediction grid (365 days)
#' @param n_sims Integer, number of posterior draws
#' @param sim_mat Optional pre-computed simulation matrix
#' @export
estimate_peak <- function(model, newdata, n_sims = 1000, sim_mat = NULL) {
  if (is.null(sim_mat)) sim_mat <- simulate_posterior_density(model, newdata, n_sims)
  res <- apply(sim_mat, 2, which.max)
  return(list(Mean = mean(res), Lower_CI = quantile(res, 0.025), Upper_CI = quantile(res, 0.975), Posterior_Draws = res))
}
