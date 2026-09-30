#' Plot Metric Gradient
#' @export
plot_metric_gradient <- function(model, dataset_for_quantiles, metric = "centroid", 
                                 x_axis = "latitude", hold_constant = "mean", 
                                 n_sims = 1000, add = FALSE, color = "#4575b4") {
                                   
  if (x_axis == "latitude") {
    x_seq <- seq(min(dataset_for_quantiles$Latitude, na.rm=TRUE), max(dataset_for_quantiles$Latitude, na.rm=TRUE), length.out = 30)
    
    t_val <- if(is.numeric(hold_constant)) hold_constant else {
      if(hold_constant == "min") quantile(dataset_for_quantiles$TMAM, 0.1, na.rm=TRUE)
      else if(hold_constant == "max") quantile(dataset_for_quantiles$TMAM, 0.9, na.rm=TRUE)
      else mean(dataset_for_quantiles$TMAM, na.rm=TRUE)
    }
    nd <- data.frame(Yday = rep(1:365, 30), LogEffort = 0, Latitude = rep(x_seq, each = 365), 
                     TMAM = t_val, fYear = model$model$fYear[1], SubSite = model$model$SubSite[1])
                     
  } else {
    x_seq <- seq(min(dataset_for_quantiles$TMAM, na.rm=TRUE), max(dataset_for_quantiles$TMAM, na.rm=TRUE), length.out = 30)
    
    l_val <- if(is.numeric(hold_constant)) hold_constant else {
      if(hold_constant == "min") quantile(dataset_for_quantiles$Latitude, 0.1, na.rm=TRUE)
      else if(hold_constant == "max") quantile(dataset_for_quantiles$Latitude, 0.9, na.rm=TRUE)
      else mean(dataset_for_quantiles$Latitude, na.rm=TRUE)
    }
    nd <- data.frame(Yday = rep(1:365, 30), LogEffort = 0, Latitude = l_val, 
                     TMAM = rep(x_seq, each = 365), fYear = model$model$fYear[1], SubSite = model$model$SubSite[1])
  }
  
  # Base mean profile for DOBY centering
  base_nd <- data.frame(Yday=1:365, LogEffort=0, Latitude=mean(dataset_for_quantiles$Latitude, na.rm=TRUE), TMAM=mean(dataset_for_quantiles$TMAM, na.rm=TRUE), fYear=model$model$fYear[1], SubSite=model$model$SubSite[1])
  
  set.seed(42)
  beta <- coef(model)
  V <- vcov(model, unconditional = TRUE)
  br <- mgcv::rmvn(n_sims, beta, V)
  
  X <- predict(model, newdata = nd, type = "lpmatrix")
  base_X <- predict(model, newdata = base_nd, type = "lpmatrix")
  
  re_cols <- grep("s\\(fYear\\)|s\\(SubSite\\)", colnames(X))
  if(length(re_cols)>0) {
    X[, re_cols] <- 0
    base_X[, re_cols] <- 0
  }
  
  sim_link <- X %*% t(br)
  sim_val <- exp(sim_link)
  
  base_dens <- apply(exp(base_X %*% t(br)), 1, mean)
  D_min <- which.min(base_dens)
  
  metric_res <- data.frame(X_Val = x_seq, Mean = 0, Lower = 0, Upper = 0)
  
  for (i in 1:30) {
    idx <- ((i-1)*365 + 1):(i*365)
    point_sims <- sim_val[idx, ]
    point_sims <- apply(point_sims, 2, function(x) if(sum(x)>0) x/sum(x) else rep(1/365, 365))
    
    if (metric == "centroid") {
      res_sims <- numeric(n_sims)
      for (s in 1:n_sims) {
        td <- point_sims[, s]
        shifted <- td[order(( (1:365 - D_min) %% 365 ) + 1)]
        doby_c <- sum((1:365)*shifted)
        res_sims[s] <- (doby_c + D_min - 1) %% 365
      }
      metric_res$Mean[i] <- mean(res_sims)
      metric_res$Lower[i] <- quantile(res_sims, 0.025)
      metric_res$Upper[i] <- quantile(res_sims, 0.975)
    }
  }
  
  layer_ribbon <- ggplot2::geom_ribbon(data = metric_res, ggplot2::aes(x = X_Val, ymin = Lower, ymax = Upper), fill = color, alpha = 0.2, inherit.aes = FALSE)
  layer_line <- ggplot2::geom_line(data = metric_res, ggplot2::aes(x = X_Val, y = Mean), color = color, linewidth = 1, inherit.aes = FALSE)
  
  if (add) {
    return(list(layer_ribbon, layer_line))
  } else {
    p <- ggplot2::ggplot() + layer_ribbon + layer_line + ggplot2::theme_minimal() +
      ggplot2::labs(x = tools::toTitleCase(x_axis), y = paste(tools::toTitleCase(metric), "DOY"))
    return(p)
  }
}
