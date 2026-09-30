#' Estimate Phenological Active Duration
#' 
#' @param model A fitted mgcv GAM/BAM object
#' @param newdata A data.frame with prediction grid (365 days)
#' @param n_sims Integer, number of posterior draws
#' @param sim_mat Optional pre-computed simulation matrix
#' @param threshold_fraction The density threshold as a fraction of peak maximum (default 0.05)
#' @export
estimate_duration <- function(model, newdata, n_sims = 1000, sim_mat = NULL, threshold_fraction = 0.05) {
  if (is.null(sim_mat)) sim_mat <- simulate_posterior_density(model, newdata, n_sims)
  res <- apply(sim_mat, 2, function(x) sum(x > (threshold_fraction * max(x))))
  return(list(Mean = mean(res), Lower_CI = quantile(res, 0.025), Upper_CI = quantile(res, 0.975), Posterior_Draws = res))
}
