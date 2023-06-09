## Directly estimates the extreme QTE using the inverse propensity score weights w1, w0 for extreme 1-pn quantiles as in Firpo(2006) and Zhang (2018).
## For construction of the confidence intervals, the b out of n bootstrap explained in Zhang (2018) is used. 
## In particular, we use the same hyperparameters b and m as suggested in Zhang (2018).

qte_firpo_zhang <- function(Y, X, D, pn, CI_level=0.9, N_bootstrap, prop_scores=NULL, replacement=TRUE){
  
  ## Choice of hyperparameters
  n <- length(Y)
  b <- floor(0.4*n - max(0,n-300)/7 - 2.3/28*max(n-1000,0) - 7/40*(1-log(5000)/log(n))*max(n-5000,0))
  k <- pn*n
  tau_b1 <- k/b
  tau_n1 <- pn
  
  tau_n0 <- min(10/n, 0.1*b/n)
  k0 <- tau_n0*n
  m <- 1 + 10/k0
  tau_b0 <- n*tau_n0/b
  
  if(is.null(prop_scores)){ 
    ## Fit polynomial sieve
    hn <- 2*n^(1/11)
    prop.fit <- glm(D ~ poly(X, hn), family=binomial)
    prop_scores <- fitted(prop.fit)
  }
  
  ## Inverse propensity score weights
  w1 <- D / prop_scores
  w0 <- (1-D) / (1-prop_scores) 
  
  ## Quantile estimation
  q1 <- BMisc::getWeightedQuantiles(cvec=Y, weights=w1, tau=c(1-tau_n1, 1-tau_b1, 1-tau_n0, 1-m*tau_n0))
  q0 <- BMisc::getWeightedQuantiles(cvec=Y, weights=w0, tau=c(1-tau_n1, 1-tau_b1, 1-tau_n0, 1-m*tau_n0))
  
  ## QTE estimates
  qte_n <- q1[1] - q0[1]      
  qte_b <- q1[2] - q0[2]      
  
  ## Normalizing factor
  an <- sqrt(tau_n0*n)/max(q1[3]-q1[4], q0[3]-q0[4])      # (4.3) in zhang's paper
  
  Vb <- numeric(N_bootstrap)
  
  ## Loop over all b out of n bootstrap datasets
  for (i in 1:N_bootstrap){
    ## Sample b out of n samples
    selected <- sample(1:length(Y), size=b, replace=replacement)
    Yboot <- Y[selected]
    Dboot <- D[selected]
    w1boot <- w1[selected]
    w0boot <- w0[selected]
    
    ## Quantile estimation for bootstrap data
    q1b <- BMisc::getWeightedQuantiles(cvec=Yboot, weights=w1boot, tau=c(1-tau_b1, 1-tau_b0, 1-m*tau_b0))
    q0b <- BMisc::getWeightedQuantiles(cvec=Yboot, weights=w0boot, tau=c(1-tau_b1, 1-tau_b0, 1-m*tau_b0))
    
    ## Bootstrap QTE estimate and normalization factor 
    qte_boot <- q1b[1] - q0b[1]
    
    ab <- sqrt(tau_b0*b)/max(q1b[2]-q1b[3], q0b[2]-q0b[3])
    
    Vb[i] <- ab*(qte_boot - qte_b)
  }
  
  
  ## Calculate upper and lower bounds of confidence interval
  C_low <- quantile(Vb, 1-(1-CI_level)/2)
  C_up <- quantile(Vb, (1-CI_level)/2)
  CIupper <- qte_n - C_up/an
  CIlower <- qte_n - C_low/an
  
  res <- list(qte=qte_n, q1=q1[1], q0=q0[1], CIupper=CIupper, CIlower=CIlower)
  
  return(res)
}