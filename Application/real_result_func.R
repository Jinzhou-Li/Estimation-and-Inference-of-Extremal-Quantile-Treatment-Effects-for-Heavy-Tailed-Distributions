# Comparing results based on different estimated propensity scores:
# give the model for propensity score, df, and QTEs, return the results of our method and Zhang's.

real_data_results <- function(formula.prop, df, pn){
  
  fit <- glm(formula.prop, data=df, family="binomial")
  prop.pred <- fitted(fit)
  
  ## Extreme QTE Estimation
  Y = df$wage
  D = df$college
  ks <- 1:300
  
  set.seed(2021)    # zhang's result is random: for reproducibility we set a seed here.
  hill_list <- list()
  firpo_zhang_list <- list()
  for (i in 1:length(pn)) {
    hill_temp <- qte_extrapolation_hill(Y, X=NULL, D, pn[i], ks, CI_level=0.95, prop_scores=prop.pred)
    hill_df <- data.frame(qtes=hill_temp$qtes, CIupper=hill_temp$CIupper, CIlower=hill_temp$CIlower, gamma0=hill_temp$p0s, gamma1=hill_temp$p1s)
    hill_list[[i]] <- hill_df
    
    firpo_zhang_temp <- qte_firpo_zhang(Y, X=NULL, D, pn[i], CI_level=0.95, N_bootstrap=1000, prop_scores=prop.pred, replacement=TRUE)
    firpo_zhang_df <- data.frame(qtes=firpo_zhang_temp$qte, CIupper=firpo_zhang_temp$CIupper, CIlower=firpo_zhang_temp$CIlower)
    firpo_zhang_list[[i]] <- firpo_zhang_df 
  }
  return(list(hill_list, firpo_zhang_list))
}