## Estimating the extreme (1-pn)-QTE and its CI using the Hill type EVI estimator. 

# input: 'Y': response vector
#        'X': covariates
#        'D': treatment
#        'pn': target is (1-pn)-QTE
#        'ks': tuning parameter related to intermediate quantile (\tau_n = ks/n)
#        'CI_level': level of confidence interval
#
# output: the list 'res', which includes:
#        'qtes': point estimation of extremal QTE
#        'q1ext': estimated quantile for Y(1)
#        'q0ext': estimated quantile for Y(0)
#        'p1s': estimated extreme value index for Y(1)
#        'p0s': estimated extreme value index for Y(0)
#        'sd': estimated sd used to construct CI
#        'CIupper': upper bound of the CI
#        'CIlower': lower bound of the CI

##########
qte_extrapolation_hill <- function(Y, X, D, pn, ks, CI_level=0.9, prop_scores=NULL){
  
  n <- length(Y)
  
  qtes <- numeric(length(ks))           # For storing resulting qte
  q1ext <- numeric(length(ks))          # For storing extrapolated quantile q1(1-pn)
  q0ext <- numeric(length(ks))          # For storing extrapolated quantile q0(1-pn)
  p1s <- numeric(length(ks))            # For storing estimated extreme value index gamma1
  p0s <- numeric(length(ks))            # For storing estimated extreme value index gamma0
  variances <- numeric(length(ks))      # For storing estimated variance
  CIupper <- numeric(length(ks))        # For storing the upper bound of CI
  CIlower <- numeric(length(ks))        # For storing the lower bound of CI
  
  if(is.null(prop_scores)){ 
    ## Fit polynomial sieve
    hn <- floor(2*n^(1/11))
    prop.fit <- glm(D ~ poly(X, hn), family=binomial)
    prop_scores <- fitted(prop.fit)
  }
  
  ## Inverse propensity score weights
  w1 <- D / prop_scores
  w0 <- (1-D) / (1-prop_scores) 
  
  vw1 <- D/prop_scores^2                ## Weights for variance estimaton
  vw0 <- (1-D)/(1-prop_scores)^2
  
  ## Estimate intermediate counterfactual quantiles
  q1s <- BMisc::getWeightedQuantiles(cvec=Y, weights=w1, tau=1-ks/n)
  q0s <- BMisc::getWeightedQuantiles(cvec=Y, weights=w0, tau=1-ks/n) 
  
  ## Loop over different choices of k
  for (j in 1:length(ks)){
    k <- ks[j] 
    
    ## Extract needed intermediate quantiles
    q1 <- q1s[j]
    q0 <- q0s[j]
    
    ## The Hill type estimator requires the intermediate quantiles to be positive.
    if(q1 <= 0 || q0 <= 0 ){
      print("Intermediate quantiles are non-positive: Hill estimator is not defined!")
    }
    
    ## Hill type estimates for extreme value indices
    p1 <- sum(log(Y[which(Y>q1)]/q1)*w1[which(Y>q1)])/k
    p0 <- sum(log(Y[which(Y>q0)]/q0)*w0[which(Y>q0)])/k
    
    ## Quantile extrapolation
    dn <- k/n/pn
    qe1 <- q1*(dn)^p1
    qe0 <- q0*(dn)^p0
    
    p1s[j] <- p1
    p0s[j] <- p0
    qtes[j] <- qe1-qe0
    q1ext[j] <- qe1
    q0ext[j] <- qe0
    
    ########## Variance estimation #########
    G1 <- sum(log(Y[which(Y>q1)]/q1)^2*vw1[which(Y>q1)])/k
    G0 <- sum(log(Y[which(Y>q0)]/q0)^2*vw0[which(Y>q0)])/k
    J11 <- sum(log(Y[which(Y>q1)]/q1)*vw1[which(Y>q1)])/k
    J00 <- sum(log(Y[which(Y>q0)]/q0)*vw0[which(Y>q0)])/k
    
    Sig <- matrix(0,4,4)
    Sig[1,1] <- G1
    Sig[2,2] <- G0
    Sig[3,3] <- sum(vw1*(Y>q1))/k 
    Sig[4,4] <- sum(vw0*(Y>q0))/k
    Sig[1,3] <- Sig[3,1] <- J11
    Sig[2,4] <- Sig[4,2] <- J00
    
    B <- matrix(0, 2, 4)
    B[1,1] <- B[2,2] <- 1  
    B[1,3] <- -p1 
    B[2,4] <- -p0
    
    alpha <- qe1/qe0
    v <- c(min(1,alpha), -min(1,1/alpha))
    
    sigmasq <- t(v) %*% B %*% Sig %*% t(B) %*% v
    norm_factor <- sqrt(k)/(max(qe1, qe0))/log(dn)
    variances[j] <- sigmasq/norm_factor^2
    
    Gaussian_q <- qnorm(1-(1-CI_level)/2)
    CIupper[j] <- (qe1-qe0) + Gaussian_q * sqrt(sigmasq/norm_factor^2)
    CIlower[j] <- (qe1-qe0) - Gaussian_q * sqrt(sigmasq/norm_factor^2)
  }
  
  res <- list(qtes=qtes, q1ext=q1ext, q0ext=q0ext, p1s=p1s, p0s=p0s, sd=sqrt(variances), CIupper=CIupper, CIlower=CIlower)
  
  return(res)
}