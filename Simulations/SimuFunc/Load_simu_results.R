### extract simulation results and create dataframes for plots

Load_simu_results <- function(folderstring, pnstrings, pntext, nlist, true_qte_mat){
  
  ## Iterate over all n and pn
  df.qte.mse <- data.frame(method=character(), qte=numeric(), mse=numeric(), n=character(), pn=character())
  df.qte.coverage <- data.frame(method=character(), cover=numeric(), n=character(), pn=character(),
                                    sd_est=numeric(), boot_mean= numeric())
  for (i in 1:length(nlist)){
    for (j in 1:length(pnstrings)){
      load(paste("Simulations/",folderstring, "sim_",  pnstrings[j], "_n_", nlist[i], ".RData", sep=""))
      
      Simu_results <- Simu_results[,1:4] # extract simulation results related to Main
      
      nsimu <- nrow(Simu_results)
      nstring <- paste("n==", nlist[i], sep="")
      pnstring <- paste("p[n]==", pntext[j], sep="")
      true_qte <- true_qte_mat[i,j]
      qte_vec <- unlist(Simu_results[,1])

      ### data frame for mse plots
      method.mse <- rep(c("Extremal QTE estimator", "Extrapolation with causal Pickands estimator", "Firpo-Zhang estimator"), 
                        nsimu)
      mse_vec <- (qte_vec - true_qte)^2
      df.qte.mse.new <- data.frame(method=method.mse, qte=qte_vec, mse=mse_vec, n=nstring, pn=pnstring)
      
      df.qte.mse <- rbind(df.qte.mse, df.qte.mse.new)
      
      ### data frame for coverage plots
      method.coverage <- rep(c("Extremal QTE CI", "Zhang", 
                               "BS Hill", "BS Hill-Per",
                               "BS Pickands", "BS Pickands-Per"), 
                             nsimu)
      
      CI_list <- Simu_results[,2]
      sd_est_list <- Simu_results[,3]
      boot_mean_list <- Simu_results[,4]
      cover_vec <- sd_vec <- boot_mean_vec <- c()
      for (simu_index in 1:nsimu) {
        
        CI_temp <- CI_list[[simu_index]]
        cover_vec_temp <- as.numeric(c(CI_temp[1]>=true_qte & CI_temp[2]<=true_qte,
                       CI_temp[3]>=true_qte & CI_temp[4]<=true_qte,
                       CI_temp[5]>=true_qte & CI_temp[6]<=true_qte,
                       CI_temp[7]>=true_qte & CI_temp[8]<=true_qte,
                       CI_temp[9]>=true_qte & CI_temp[10]<=true_qte,
                       CI_temp[11]>=true_qte & CI_temp[12]<=true_qte))
        sd_vec_temp <- c(sd_est_list[[simu_index]][1], NA,
                    sd_est_list[[simu_index]][2], NA,
                    sd_est_list[[simu_index]][3], NA)
        boot_mean_vec_temp <- c(NA,NA,boot_mean_list[[simu_index]][1],NA,boot_mean_list[[simu_index]][2],NA)
      
        cover_vec <- c(cover_vec,cover_vec_temp)
        sd_vec <- c(sd_vec, sd_vec_temp)
        boot_mean_vec <- c(boot_mean_vec, boot_mean_vec_temp)
      }
      
      df.qte.coverage.new <- data.frame(method=method.coverage, cover=cover_vec, 
                                        n=nstring, pn=pnstring,
                                        sd_est=sd_vec, boot_mean= boot_mean_vec)
      df.qte.coverage <- rbind(df.qte.coverage, df.qte.coverage.new)
      
    } # end loop pn
  } # end loop n
  
  names(df.qte.mse) <- c("Method", "QTE", "SQE", "n", "pn")
  return(list(df.qte.mse, df.qte.coverage))
}


### extract simulation results and create dataframes for plots in supplements

Load_simu_results_supp <- function(folderstring, pnstrings, pntext, nlist, true_qte_mat){
  
  ## Iterate over all n and pn
  df.qte.mse <- data.frame(method=character(), qte=numeric(), mse=numeric(), n=character(), pn=character(), k=numeric())
  df.qte.coverage <- data.frame(method=character(), cover=numeric(), n=character(), pn=character(), k=numeric())
  
  for (i in 1:length(nlist)){
    
    if(nlist[i]==1000) {ks <- seq(10,200,10)}
    if(nlist[i]==2000) {ks <- seq(10,300,15)}
    if(nlist[i]==5000) {ks <- seq(10,500,20)}
    
    for (j in 1:length(pnstrings)){
      load(paste("Simulations/",folderstring, "sim_",  pnstrings[j], "_n_", nlist[i], ".RData", sep=""))
      
      Simu_results <- Simu_results[,5:11] # extract simulation results related to Supp
      
      nsimu <- nrow(Simu_results)
      nstring <- paste("n==", nlist[i], sep="")
      pnstring <- paste("p[n]==", pntext[j], sep="")
      true_qte <- true_qte_mat[i,j]
      
      qte_vec <- c(unlist(Simu_results[,1]), unlist(Simu_results[,2]), unlist(Simu_results[,3]))
      k_vec <- c(rep(ks, nsimu), rep(ks, nsimu), rep(NA, nsimu))
      
      ### data frame for mse plots
      method.mse <- c(rep("Extremal QTE estimator",length(ks)*nsimu), 
                      rep("Extrapolation with causal Pickands estimator",length(ks)*nsimu), 
                      rep("Firpo-Zhang estimator",nsimu))
      mse_vec <- (qte_vec - true_qte)^2
      df.qte.mse.new <- data.frame(method=method.mse, qte=qte_vec, mse=mse_vec, n=nstring, pn=pnstring, k=k_vec)
      
      df.qte.mse <- rbind(df.qte.mse, df.qte.mse.new)
      
      ### data frame for coverage plots
      cover_vec <- k_vec_cover <-  c()
      for (simu_index in 1:nsimu) {
        
        df_temp <- Simu_results[simu_index,]
        cover_vec_temp <- as.numeric(c(df_temp[[4]]>=true_qte & df_temp[[5]]<=true_qte,
                                       df_temp[[6]]>=true_qte & df_temp[[7]]<=true_qte))
        
        cover_vec <- c(cover_vec,cover_vec_temp)
        k_vec_cover <- c(k_vec_cover, c(ks, NA))
      }
      
      method.coverage <- rep(c(rep("Extremal QTE CI", length(ks)), "Zhang's b out of n bootstrap"),
                             nsimu)
      
      df.qte.coverage.new <- data.frame(method=method.coverage, cover=cover_vec, 
                                        n=nstring, pn=pnstring, k=k_vec_cover)
      df.qte.coverage <- rbind(df.qte.coverage, df.qte.coverage.new)
      
    } # end loop pn
  } # end loop n
  
  names(df.qte.mse) <- c("Method", "QTE", "SQE", "n", "pn", "k")
  return(list(df.qte.mse, df.qte.coverage))
}

