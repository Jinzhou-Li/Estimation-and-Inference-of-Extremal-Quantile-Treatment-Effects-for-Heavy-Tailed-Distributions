# Bootstrap CIs based on extrapolated-qte estimators.

qte_extrapolation_bootstrap <- function(extrapolator, Y, X, D, pn, ks, CI_level=0.9, N_bootstrap){
  
  if(extrapolator=="Hill" | extrapolator=="Pickands"){
    if(extrapolator=="Hill"){qte_extrapolation_func <- qte_extrapolation_hill}
    if(extrapolator=="Pickands"){qte_extrapolation_func <- qte_extrapolation_pickands}
  }else{
    stop("Input 'extrapolator' must be 'Hill' or 'Pickands'!")
  }
  
  # qte estimator based on full sample
  qtes <- qte_extrapolation_func(Y, X, D, pn, ks)$qtes
  
  bootstrapped_qtes <- matrix(0, nrow=N_bootstrap, ncol=length(ks)) 
  for (i in 1:N_bootstrap) {
    ## Sample bootstrap dataset
    selected <- sample(1:length(Y), replace=TRUE)
    Yboot <- Y[selected]
    Dboot <- D[selected]
    Xboot <- X[selected]
    
    bootstrapped_qtes[i,] <- qte_extrapolation_func(Yboot, Xboot, Dboot, pn, ks)$qtes
  }
  
  # to check the bias and variance of the bootstrap CIs.
  sd_boot <- apply(bootstrapped_qtes, 2, sd)
  mean_boot <- apply(bootstrapped_qtes, 2, mean)
  
  ## Bootstrap CIs
  Gaussian_q <- qnorm(1-(1-CI_level)/2)
  CIupper_Gaussian <- qtes + Gaussian_q * sd_boot
  CIlower_Gaussian <- qtes - Gaussian_q * sd_boot
  
  quantile_up <- 1-(1-CI_level)/2
  lower_up <- (1-CI_level)/2
  CIupper_quantile <- apply(bootstrapped_qtes, 2, quantile, quantile_up)
  CIlower_quantile <- apply(bootstrapped_qtes, 2, quantile, lower_up)
  
  res <- list(qtes=qtes, mean_boot= mean_boot, sd_boot=sd_boot,
              CIupper_Gaussian=CIupper_Gaussian, CIlower_Gaussian=CIlower_Gaussian,
              CIupper_quantile=CIupper_quantile, CIlower_quantile=CIlower_quantile)
  
  return(res)
}
