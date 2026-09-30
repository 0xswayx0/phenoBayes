# Internal helper to generate posterior density curves
simulate_posterior_density <- function(model, newdata, n_sims = 1000, seed = 42) {
  set.seed(seed)
  beta <- coef(model)
  V <- vcov(model, unconditional = TRUE)
  br <- mgcv::rmvn(n_sims, beta, V)
  
  X <- predict(model, newdata = newdata, type = "lpmatrix")
  re_cols <- grep("s\\(fYear\\)|s\\(SubSite\\)", colnames(X))
  if (length(re_cols) > 0) X[, re_cols] <- 0
  
  sim_link <- X %*% t(br)
  sim_val <- exp(sim_link)
  
  # Normalize
  apply(sim_val, 2, function(x) {
    if (sum(x) > 0) x / sum(x) else rep(1/nrow(newdata), nrow(newdata))
  })
}
