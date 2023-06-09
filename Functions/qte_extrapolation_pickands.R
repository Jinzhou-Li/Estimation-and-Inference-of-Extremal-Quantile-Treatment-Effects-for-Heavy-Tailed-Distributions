## Estimating the extreme (1-pn)-QTE and its CI using the Pickands type EVI estimator. 

# input: 'Y': response vector
#        'X': covariates
#        'D': treatment
#        'pn': target is (1-pn)-QTE
#        'ks': tuning parameter related to intermediate quantile (\tau_n = ks/n)
#
# output: the list 'res', which includes:
#        'qtes': point estimation of extremal QTE
#        'q1ext': estimated quantile for Y(1)
#        'q0ext': estimated quantile for Y(0)
#        'p1s': estimated extreme value index for Y(1)
#        'p0s': estimated extreme value index for Y(0)

###################
qte_extrapolation_pickands <- function(Y, X, D, pn, ks, prop_scores=NULL){
  
  n <- length(Y)
  
  qtes <- numeric(length(ks))           # For storing resulting qte
  q1ext <- numeric(length(ks))          # For storing extrapolated quantile q1(1-pn)
  q0ext <- numeric(length(ks))          # For storing extrapolated quantile q0(1-pn)
  p1s <- numeric(length(ks))            # For storing estimated extreme value index gamma1
  p0s <- numeric(length(ks))            # For storing estimated extreme value index gamma0

  if(is.null(prop_scores)){ 
    ## Fit polynomial sieve
    hn <- 2*n^(1/11)
    prop.fit <- glm(D ~ poly(X, hn), family=binomial)
    prop_scores <- fitted(prop.fit)
  }
  
  ## Inverse propensity score weights
  w1 <- D / prop_scores
  w0 <- (1-D) / (1-prop_scores) 
  
  ## Estimate intermediate counterfactual quantiles
  q1s <- BMisc::getWeightedQuantiles(cvec=Y, weights=w1, tau=c(1-ks/n, 1-2*ks/n, 1-4*ks/n))
  q0s <- BMisc::getWeightedQuantiles(cvec=Y, weights=w0, tau=c(1-ks/n, 1-2*ks/n, 1-4*ks/n)) 
  
  ## Loop over different choices of k
  for (j in 1:length(ks)){
    k <- ks[j] 
    
    ## Extract needed intermediate quantiles
    q1 <- c(q1s[j], q1s[j+length(ks)], q1s[j+2*length(ks)])
    q0 <- c(q0s[j], q0s[j+length(ks)], q0s[j+2*length(ks)])
    
    # Pickands type estimates of the extreme value indices
    p1 <- log((q1[1]-q1[2])/(q1[2]-q1[3]))/log(2)  
    p0 <- log((q0[1]-q0[2])/(q0[2]-q0[3]))/log(2) 
    
    ## Quantile extrapolation
    dn <- k/n/pn
    qe1 <- q1[1]*dn^p1
    qe0 <- q0[1]*dn^p0
    
    p1s[j] <- p1
    p0s[j] <- p0
    qtes[j] <- qe1-qe0
    q1ext[j] <- qe1
    q0ext[j] <- qe0
  }
  
  res <- list(qtes=qtes, q1ext=q1ext, q0ext=q0ext, p1s=p1s, p0s=p0s)
  
  return(res)
}
