# Load functions and simulation settings
source("LoadFunctions.R")

## Simulation parameters
setting="H1"
nlist <- c(1000, 2000, 5000)
CI_level=0.9
N_bootstrap = 1000
nsimu <- 1000
para_cores <- 50
set.seed(2023)

## Loop over different choices of the sample size n
for(j in 1:length(nlist)){
  n <- nlist[j]
  print(paste("----------------", n, "------------------"))
  hn <- floor(2*n^(1/11))
  k <- n^(0.65)
  
  if(n==1000) {ks <- seq(10,200,10)}
  if(n==2000) {ks <- seq(10,300,15)}
  if(n==5000) {ks <- seq(10,500,20)}
  
  ## Parallel simulations for pn=5/n 
  pn <- 5/n
  Simu_results <- t(mcmapply(Simu_func, 1:nsimu, 
                             MoreArgs=list(setting, n, pn, k, ks, hn, CI_level, N_bootstrap), 
                             mc.cores=para_cores))
  save(Simu_results, file=paste("Simulations/",setting,"/sim_pn5_n_", n, ".RData", sep="")) ## Save simulation data
  
  ## Parallel simulations for pn=1/n 
  pn <- 1/n
  Simu_results <- t(mcmapply(Simu_func, 1:nsimu, 
                             MoreArgs=list(setting, n, pn, k, ks, hn, CI_level, N_bootstrap), 
                             mc.cores=para_cores))
  save(Simu_results, file=paste("Simulations/",setting,"/sim_pn1_n_", n, ".RData", sep="")) ## Save simulation data
  
  # Parallel simulations for pn=5/nlogn
  pn <- 5/n/log(n)
  Simu_results <- t(mcmapply(Simu_func, 1:nsimu, 
                             MoreArgs=list(setting, n, pn, k, ks, hn, CI_level, N_bootstrap), 
                             mc.cores=para_cores))
  save(Simu_results, file=paste("Simulations/",setting,"/sim_pn5nlogn_n_", n, ".RData", sep="")) ## Save simulation data
}
