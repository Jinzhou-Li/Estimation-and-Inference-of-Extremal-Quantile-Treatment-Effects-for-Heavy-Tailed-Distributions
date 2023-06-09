library(tidyverse)
library(grf)
library(here)
library("ggpubr")
library(evd)
library(ismev)
library(future)
library(doFuture)

################
n <- 1000
df <- 3
n_sim <- 1000
threshold <- 0.9

seedR <- 1232
set.seed(seedR)

levels <- c(0.99,0.995, 0.999,0.9995, 0.9999, 0.99995)
levels <- c(0.9,0.999, 0.9999, 0.99999)
n_levels <- length(levels)

q_GPD <- function(p, p0, t_x0, sigma, xi){
  ## produce the estimated extreme quantiles of GPD
  (((1-p)/(1-p0))^{-xi} - 1) * (sigma / xi) + t_x0
}

predict_Hill_quantiles <- function(Y, threshold, p = c(0.95)){
  p0 <- threshold
  t0 <- quantile(Y, p0)
  k <- sum(Y > t0)
  n <- length(Y)
  
  ## The Hill type estimator requires the intermediate quantiles to be positive.
  if(t0 <= 0){
    print("Intermediate quantiles are non-positive: Hill estimator is not defined!")
  }
  
  ## Hill type estimates for extreme value indices
  hill <- sum(log(Y[which(Y>t0)]/t0))/k
  
  ## Quantile extrapolation
  dn <- k/n/(1-p)
  qe <- t0*(dn)^hill
  
  return(qe)
}

predict_unconditional_quantiles <- function(threshold, quantiles = c(0.95),
                                            Y, ntest){
  p0 <- threshold
  t0 <- quantile(Y, p0)
  pars <- ismev::gpd.fit(Y, t0, show = FALSE)$mle
  sigma <- pars[1]
  xi <- pars[2]
  
  q_hat <- q_GPD(quantiles, p0, t0, sigma, xi)
  
  predictions <- matrix(q_hat, nrow = ntest, ncol = length(quantiles), byrow = T)
  pars <- cbind(rep(sigma, ntest) , rep(xi, ntest))
  return(list(predictions = predictions, pars = pars))
}


results <- foreach(i=1:n_sim, .errorhandling="stop", .combine=bind_rows) %do% {
  
  Y <- rfrechet(n, loc = 0, scale = 1, shape = df) 
  
  true_Q <-  qfrechet(p = levels, loc = 0, scale = 1, shape = df) 
  names(true_Q) <- paste0("True_",levels*100)
  
  empirical_Q <- quantile(Y,probs=levels)
  names(empirical_Q) <- paste0("Empirical_",levels*100)
  
  pot_Q <- c(predict_unconditional_quantiles(threshold=threshold, quantiles = levels, Y = Y, ntest = 1)$predictions)
  hill_Q <- predict_Hill_quantiles(Y, threshold = threshold, p = levels)
  names(pot_Q) <- paste0("PoT_",levels*100)
  names(hill_Q) <- paste0("hill_",levels*100)
  
  table_i <- bind_cols(as_tibble_row(true_Q),
                       as_tibble_row(empirical_Q),
                       as_tibble_row(hill_Q))
  table_i
}

means <- sapply(results, function(x){mean(x)})
sds <- sapply(results, function(x){sd(x)})
sems <- sds/sqrt(n_sim)

meanq_methods_plot <- data.frame(q=factor(rep(paste0(100*levels,"% (", round(1000*(1-levels),2), ")"),3)), 
                                 mean_quantile=(means),
                                  Method=rep(c("True","Empirical", "Hill"), each=n_levels), err=sems) %>%
  ggplot( aes(x=q, y=mean_quantile, group=Method, color=Method, linetype=Method, fill=Method)) +
  geom_line(size=1) + geom_point() + 
  scale_linetype_manual(values=c("dotted", "dashed", "solid", "dotted")) +
  scale_color_manual(values = c("darkred","darkgreen","black")) +
  labs(title=NULL, x= expression(paste("Probability level ", 1-p[n]," (effective sample size ", np[n],")")), 
       y= "Averaged quantile estimate", color=NULL, fill=NULL, linetype=NULL) + 
  theme_bw(base_size = 15) + theme(legend.position="bottom")
plot(meanq_methods_plot)

ggsave("Figure1/Figure1.pdf", plot=meanq_methods_plot, device="pdf", width=150, height=100, units="mm", dpi=300)