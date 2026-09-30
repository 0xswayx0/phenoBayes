#' Plot Phenology Curve
#' @export
plot_phenology_curve <- function(model, dataset_for_quantiles, latitude = "mean", temperature = "mean", 
                                 n_sims = 1000, log = TRUE, add = FALSE, color = "#d73027") {
  
  lat_val <- if(is.numeric(latitude)) latitude else {
    if(latitude == "min") quantile(dataset_for_quantiles$Latitude, 0.1, na.rm=TRUE)
    else if(latitude == "max") quantile(dataset_for_quantiles$Latitude, 0.9, na.rm=TRUE)
    else mean(dataset_for_quantiles$Latitude, na.rm=TRUE)
  }
  
  temp_val <- if(is.numeric(temperature)) temperature else {
    if(temperature == "min") quantile(dataset_for_quantiles$TMAM, 0.1, na.rm=TRUE)
    else if(temperature == "max") quantile(dataset_for_quantiles$TMAM, 0.9, na.rm=TRUE)
    else mean(dataset_for_quantiles$TMAM, na.rm=TRUE)
  }
  
  nd <- data.frame(Yday = 1:365, LogEffort = 0, Latitude = lat_val, TMAM = temp_val, 
                   fYear = model$model$fYear[1], SubSite = model$model$SubSite[1])
                   
  sim_mat <- simulate_posterior_density(model, nd, n_sims)
  
  df_plot <- data.frame(
    Yday = 1:365,
    Density = rowMeans(sim_mat),
    Lower = apply(sim_mat, 1, quantile, probs = 0.025),
    Upper = apply(sim_mat, 1, quantile, probs = 0.975)
  )
  
  if(log) {
    df_plot$Density <- df_plot$Density + 1e-5
    df_plot$Lower <- df_plot$Lower + 1e-5
    df_plot$Upper <- df_plot$Upper + 1e-5
  }
  
  layer_ribbon <- ggplot2::geom_ribbon(data = df_plot, ggplot2::aes(x = Yday, ymin = Lower, ymax = Upper), fill = color, alpha = 0.2, inherit.aes = FALSE)
  layer_line <- ggplot2::geom_line(data = df_plot, ggplot2::aes(x = Yday, y = Density), color = color, linewidth = 1, inherit.aes = FALSE)
  
  if (add) {
    return(list(layer_ribbon, layer_line))
  } else {
    p <- ggplot2::ggplot() + layer_ribbon + layer_line + ggplot2::theme_minimal() +
      ggplot2::labs(x = "Day of Year", y = if(log) "Probability Density (log10)" else "Probability Density")
    if (log) p <- p + ggplot2::scale_y_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1), labels = c("0.0001", "0.001", "0.01", "0.1"))
    return(p)
  }
}
