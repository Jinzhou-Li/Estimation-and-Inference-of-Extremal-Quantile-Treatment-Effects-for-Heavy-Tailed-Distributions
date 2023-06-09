############## Setting H1 in the paper
h1.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  S <- rt(n, df=3)
  Y1 <- 5*S*(1+X)
  Y0 <- S*(1+X)
  
  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}

############## Setting H2 in the paper
h2.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  C1 <- rfrechet(n, shape=2)
  C2 <- rfrechet(n, shape=3)
  Y1 <- C1*exp(X)
  Y0 <- C2*exp(X)
  
  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}

############## Setting H3 in the paper
h3.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  Y1 <- rpareto(n, shape=1.75+X, scale=2)
  Y0 <- rpareto(n, shape=1.75+5*X, scale=1) 
  
  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}

###################### 
Simu_func <- function(simu_index, setting, n, pn, k, ks, hn, CI_level, N_bootstrap){
  
  ## Sample data
  if(setting=="H1"){sim_sample <- h1.fun(n)}
  if(setting=="H2"){sim_sample <- h2.fun(n)}
  if(setting=="H3"){sim_sample <- h3.fun(n)}
  Y1 <- sim_sample$Y1
  Y0 <- sim_sample$Y0
  D <- sim_sample$D
  X <- sim_sample$X
  Y <- D*Y1 + (1-D)*Y0
  
  ## use the same propensity scores for all methods
  # hn <- floor(2*n^(1/11))
  prop.fit <- glm(D ~ poly(X, hn), family=binomial)
  prop_scores <- fitted(prop.fit)
  
  ################## Simulations in the main paper
  ## run different methods with one value 'k'
  result_hill <- qte_extrapolation_hill(Y, X, D, pn, k, CI_level, prop_scores=prop_scores)
  result_pickands <- qte_extrapolation_pickands(Y, X, D, pn, k, prop_scores=prop_scores)
  result_firpo_zhang <- qte_firpo_zhang(Y, X, D, pn, CI_level, N_bootstrap, prop_scores=prop_scores, replacement=TRUE)
  
  result_bootstrap_hill <- qte_extrapolation_bootstrap("Hill", Y, X, D, pn, k, CI_level, N_bootstrap)
  result_bootstrap_pickands <- qte_extrapolation_bootstrap("Pickands", Y, X, D, pn, k, CI_level, N_bootstrap)
  
  ## Save results
  qte_est_vec <- c(result_hill$qtes, result_pickands$qtes, result_firpo_zhang$qte)
  CI_vec <- c(result_hill$CIupper, result_hill$CIlower,
              result_firpo_zhang$CIupper, result_firpo_zhang$CIlower,
              result_bootstrap_hill$CIupper_Gaussian, result_bootstrap_hill$CIlower_Gaussian,
              result_bootstrap_hill$CIupper_quantile, result_bootstrap_hill$CIlower_quantile,
              result_bootstrap_pickands$CIupper_Gaussian, result_bootstrap_pickands$CIlower_Gaussian,
              result_bootstrap_pickands$CIupper_quantile, result_bootstrap_pickands$CIlower_quantile)
  sd_est_vec <- c(result_hill$sd, result_bootstrap_hill$sd_boot, result_bootstrap_pickands$sd_boot)
  boot_mean_vec <- c(result_bootstrap_hill$mean_boot, result_bootstrap_pickands$mean_boot)
  
  ################## Simulations in Supplement
  ## run different methods with vector 'ks'
  result_hill_supp <- qte_extrapolation_hill(Y, X, D, pn, ks, CI_level, prop_scores=prop_scores)
  result_pickands_supp <- qte_extrapolation_pickands(Y, X, D, pn, ks, prop_scores=prop_scores)
  
  ## Save results
  qte_est_hill_supp <- result_hill_supp$qtes
  qte_est_pickands_supp <- result_pickands_supp$qtes
  qte_est_firpo_zhang_supp <- result_firpo_zhang$qte
  CI_upper_hill_supp <- result_hill_supp$CIupper
  CI_lower_hill_supp <- result_hill_supp$CIlower
  CI_upper_firpo_zhang_supp <- result_firpo_zhang$CIupper
  CI_lower_firpo_zhang_supp <- result_firpo_zhang$CIlower
  
  # put all results in one list
  res_list <- list(qte_est_vec, CI_vec, sd_est_vec, boot_mean_vec, # below are for supplements
                   qte_est_hill_supp, qte_est_pickands_supp, qte_est_firpo_zhang_supp,
                   CI_upper_hill_supp, CI_lower_hill_supp,
                   CI_upper_firpo_zhang_supp, CI_lower_firpo_zhang_supp)
  return(res_list)
}