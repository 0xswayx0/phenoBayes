#' Estimate Phenological Modality (Number of Peaks)
#' 
#' @param model A fitted mgcv GAM/BAM object
#' @param newdata A data.frame with prediction grid (365 days)
#' @param n_sims Integer, number of posterior draws
#' @param sim_mat Optional pre-computed simulation matrix
#' @param rel_prominence Relative prominence threshold for valid peaks (default 0.25)
#' @export
estimate_modality <- function(model, newdata, n_sims = 1000, sim_mat = NULL, rel_prominence = 0.25) {
  if (is.null(sim_mat)) sim_mat <- simulate_posterior_density(model, newdata, n_sims)
  
  count_peaks <- function(dens) {
    diffs <- diff(dens); signs <- sign(diffs); signs[signs == 0] <- -1
    peaks <- which(diff(signs) == -2) + 1
    if (length(peaks) <= 1) return(max(1, length(peaks)))
    valid_peaks <- peaks[1]
    for (i in 2:length(peaks)) {
      trough_val <- min(dens[valid_peaks[length(valid_peaks)]:peaks[i]])
      peak_min_height <- min(dens[valid_peaks[length(valid_peaks)]], dens[peaks[i]])
      if ((peak_min_height - trough_val) / peak_min_height >= rel_prominence) {
        valid_peaks <- c(valid_peaks, peaks[i])
      } else {
        if (dens[peaks[i]] > dens[valid_peaks[length(valid_peaks)]]) valid_peaks[length(valid_peaks)] <- peaks[i]
      }
    }
    return(length(valid_peaks))
  }
  
  res <- apply(sim_mat, 2, count_peaks)
  return(list(Mean = mean(res), Lower_CI = quantile(res, 0.025), Upper_CI = quantile(res, 0.975), Posterior_Draws = res))
}
