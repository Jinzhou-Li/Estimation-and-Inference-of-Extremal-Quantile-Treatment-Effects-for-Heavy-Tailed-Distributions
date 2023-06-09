calculate_true_qtes <- function(sampling.fun, pn, Nsim=1e8, ran.seed=1){
  ## Approximates the true QTE by sampling from the Model and calculating the quantiles
  set.seed(ran.seed)
  
  samples <- sampling.fun(Nsim)
  q1 = quantile(samples$Y1, probs = c(1-pn))
  q0 = quantile(samples$Y0, probs = c(1-pn))
  return(matrix(q1 - q0, nrow=3))  # here nrow=3 only for our simulations.
}
