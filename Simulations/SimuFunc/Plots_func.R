# Functions to get plots
library(ggplot2)

######### Create MSE boxplots
squared_error_boxplot <- function(df.qte.mse, colors){
  ## This function allows us to change the whiskers of the boxplot to the 0.1 and the 0.9 quantile.
  bp.pctiles = function (x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    r <- quantile(x, probs = probs, na.rm = TRUE)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  # reorder methods in plot
  df.qte.mse$Method <- factor(df.qte.mse$Method, levels=c("Firpo-Zhang estimator", "Extrapolation with causal Pickands estimator", "Extremal QTE estimator"))
  df.qte.mse$pn <- factor(df.qte.mse$pn, levels=c("p[n]==5/n", "p[n]==1/n", "p[n]==5/nlogn"))
  
  sqe.boxplot <- ggplot(df.qte.mse, aes(x=Method, y=SQE, fill=Method)) +
    stat_summary(fun.data=bp.pctiles, geom="boxplot", width=0.5, position=position_dodge(width=0.5)) +
    stat_summary(fun=mean, geom="point", shape=15,  position=position_dodge(width=0.5)) +
    scale_fill_manual(values=colors) +
    facet_grid(pn ~ n, scales="free", labeller=label_parsed) + 
    theme_bw() +
    coord_trans(y="log10") +
    ylab("Squared Error") +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
    scale_y_continuous(breaks=c(1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11)) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(text = element_text(size = 22))
  
  return(sqe.boxplot)
}

########## Create Coverage plots
coverage_plot_func <- function(df.qte.coverage, colors, 
                               coverage.target, CI.level.simu=0.95){
  
  df.coverage <- data.frame(Method=character(), Coverage=numeric(), 
                            Upper=numeric(), Lower=numeric(), 
                            n=character(), pn=character(),
                            sd_est_mean=numeric(), boot_mean= numeric())
  
  Gaussian_q <- qnorm(1-(1-CI.level.simu)/2)
  method_vec <- unique(df.qte.coverage$method)
  n_vec <- unique(df.qte.coverage$n)
  pn_vec <- unique(df.qte.coverage$pn)
  for (i in 1:length(method_vec)) {
    for (j in 1:length(n_vec)) {
      for(k in 1:length(pn_vec)){
        index_temp <- which(df.qte.coverage$method == method_vec[i]
                            & df.qte.coverage$n == n_vec[j]
                            & df.qte.coverage$pn == pn_vec[k])
        df_temp <- df.qte.coverage[index_temp,]
        coverage_temp <- mean(df_temp$cover)
        sd_est_temp <- mean(df_temp$sd_est)
        boot_mean_temp <- mean(df_temp$boot_mean)
        
        nsimu <- length(df_temp$cover)
        ci.upper <- coverage_temp + Gaussian_q*sqrt(coverage_temp*(1-coverage_temp)/nsimu)
        ci.lower <- coverage_temp - Gaussian_q*sqrt(coverage_temp*(1-coverage_temp)/nsimu)
      
        df.coverage.temp <- data.frame(Method=method_vec[i], Coverage=coverage_temp,
                                       Upper=ci.upper, Lower= ci.lower,
                                       n=n_vec[j], pn=pn_vec[k],
                                       sd_est_mean=sd_est_temp, boot_mean=boot_mean_temp)
        
        df.coverage <- rbind(df.coverage, df.coverage.temp)
      }
    }
  } # end for-loops
  
  # reorder methods in plot
  df.coverage$pn <- factor(df.coverage$pn, levels=c("p[n]==5/n", "p[n]==1/n", "p[n]==5/nlogn"))
  df.coverage$Method <- factor(df.coverage$Method, 
                               levels=c("Zhang", 
                                        "BS Pickands", "BS Pickands-Per",
                                        "BS Hill", "BS Hill-Per",
                                        "Extremal QTE CI"))
  
  #### coverage plot: main
  df_main <- df.coverage[df.coverage$Method %in% c("Zhang", "BS Pickands", "BS Hill", "Extremal QTE CI"),]
  coverage.plot.main <- ggplot(df_main, aes(x=Method, y=Coverage, col=Method)) +
    geom_point() + geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.5, position=position_dodge(0.05)) +
    geom_hline(yintercept = coverage.target) + 
    facet_grid(pn ~ n, scales="free", labeller=label_parsed) +
    scale_color_manual(values=colors[c(1,2,4,6)]) +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="None") +
    guides(colour = guide_legend(override.aes = list(shape = NA), title="")) +
    scale_y_continuous(limits=c(0.6, 1)) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme(text = element_text(size = 22))
  
  #### coverage plot: comparing two bootstrap CIs
  df_boot <- df.coverage[df.coverage$Method %in% c("BS Pickands", "BS Pickands-Per",
                                                   "BS Hill", "BS Hill-Per"),]
  
  coverage.plot.boot <- ggplot(df_boot, aes(x=Method, y=Coverage, col=Method)) +
    geom_point() + geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.5, position=position_dodge(0.05)) +
    geom_hline(yintercept = coverage.target) + 
    facet_grid(pn ~ n, scales="free", labeller=label_parsed) +
    scale_color_manual(values=colors[c(2,3,4,5)]) +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="None") +
    guides(colour = guide_legend(override.aes = list(shape = NA), title="")) +
    scale_y_continuous(limits=c(0.6, 1)) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme(text = element_text(size = 18))
  
  return(list(coverage.plot.main, coverage.plot.boot, df.coverage))
}


####### Create QQ plots
qqplot_func <- function(df.qte.mse, nlist, pnstrings){
  qqplot_list <- list()
  count <- 0
  
  # reorder methods in plot
  df.qte.mse$Method <- factor(df.qte.mse$Method, levels=c("Firpo-Zhang estimator", "Extrapolation with causal Pickands estimator", "Extremal QTE estimator"))
  
  for (i in 1:length(nlist)){
    for (j in 1:length(pnstrings)){
      nstring <- paste("n==", nlist[i], sep="")
      pnstring <- paste("p[n]==", pntext[j], sep="")
      target_index <- which(df.qte.mse$n == nstring & df.qte.mse$pn ==pnstring)
      df.qte.new <- df.qte.mse[target_index,]
      
      qqplot <- ggplot(df.qte.new, aes(sample=QTE)) + stat_qq() + stat_qq_line() + 
        facet_wrap(. ~ Method, scale="free", ncol=3) + theme_bw() + ylim(-1000,3000) +
        theme(plot.title = element_text(hjust = 0.5)) + labs(x="Theoretical Normal Quantiles", y="Sample QTE Quantiles")
      
      count <- count+1
      qqplot_list[[count]] <- qqplot
    }
  } # end for loops
  return(qqplot_list)
}

######### Create MSE boxplots for supplement
squared_error_plot_supp <- function(df.qte.mse, colors){
  
  df.mse.final <- data.frame(Method=character(), qte=numeric(), mse=numeric(), n=character(), pn=character(), k=numeric())
  method_vec <- unique(df.qte.mse$Method)
  n_vec <- unique(df.qte.mse$n)
  pn_vec <- unique(df.qte.mse$pn)
  
  for (i in 1:length(method_vec)) {
    for (j in 1:length(n_vec)) {
      if(n_vec[j]=="n==1000") {k_vec <- seq(10,200,10); n_num_temp <- 1000}
      if(n_vec[j]=="n==2000") {k_vec <- seq(10,300,15); n_num_temp <- 2000}
      if(n_vec[j]=="n==5000") {k_vec <- seq(10,500,20); n_num_temp <- 5000}
      for (k_temp in 1:length(k_vec)) {
        for(l in 1:length(pn_vec)){
          
          # as Firpo-Zhang has k value NA
          if(method_vec[i]=="Firpo-Zhang estimator"){
            index_temp <- which(df.qte.mse$Method == method_vec[i]
                                & df.qte.mse$n == n_vec[j]
                                & df.qte.mse$pn == pn_vec[l])
          }else{
            index_temp <- which(df.qte.mse$Method == method_vec[i]
                                & df.qte.mse$n == n_vec[j]
                                & df.qte.mse$pn == pn_vec[l]
                                & df.qte.mse$k == k_vec[k_temp])
          }
          
          df_temp <- df.qte.mse[index_temp,]
          mse_temp <- mean(df_temp$SQE)
          
          df.mse.temp <- data.frame(Method=method_vec[i], MSE=mse_temp,
                                         n=n_vec[j], pn=pn_vec[l], k=k_vec[k_temp],
                                    n_numeric=n_num_temp)
          df.mse.final <- rbind(df.mse.final, df.mse.temp)
        }
      }
    }
  } # end for-loops
  
  # reorder methods in plot
  df.mse.final$Method <- factor(df.mse.final$Method, levels=c("Firpo-Zhang estimator", "Extrapolation with causal Pickands estimator", "Extremal QTE estimator"))
  df.mse.final$pn <- factor(df.mse.final$pn, levels=c("p[n]==5/n", "p[n]==1/n", "p[n]==5/nlogn"))
  
  mse.plot <- ggplot(df.mse.final, aes(x=k, y=MSE, group=Method, col=Method, linetype=Method)) + 
    geom_line(lwd=1) +
    facet_grid(pn ~ n, scales="free", labeller=label_parsed) + 
    scale_color_manual(values=colors) +
    scale_linetype_manual(values=c("dashed", "dotted", "solid")) +
    labs(color  = "", linetype = "") +
    scale_y_log10() +
    geom_vline(data = df.mse.final,
               aes(xintercept = n_numeric^(0.65)), linetype="dashed", size=0.3) +
    theme_bw() + 
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom") +
    guides(colour = guide_legend(override.aes = list(shape = NA), title="")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 16))
  
  return(mse.plot)
}


######### Create coverage boxplots for supplement
coverage_plot_supp <- function(df.qte.coverage, coverage.colors, 
                               coverage.target){
  
  df.coverage.final <- data.frame(method=character(), coverage=numeric(), n=character(), pn=character(), k=numeric())
  method_vec <- unique(df.qte.coverage$method)
  n_vec <- unique(df.qte.coverage$n)
  pn_vec <- unique(df.qte.coverage$pn)
  
  for (i in 1:length(method_vec)) {
    for (j in 1:length(n_vec)) {
      if(n_vec[j]=="n==1000") {k_vec <- seq(10,200,10); n_num_temp <- 1000}
      if(n_vec[j]=="n==2000") {k_vec <- seq(10,300,15); n_num_temp <- 2000}
      if(n_vec[j]=="n==5000") {k_vec <- seq(10,500,20); n_num_temp <- 5000}
      for (k_temp in 1:length(k_vec)) {
        for(l in 1:length(pn_vec)){
          
          # as 'Zhang's b out of n bootstrap' has k value NA
          if(method_vec[i]=="Zhang's b out of n bootstrap"){
            index_temp <- which(df.qte.coverage$method == method_vec[i]
                                & df.qte.coverage$n == n_vec[j]
                                & df.qte.coverage$pn == pn_vec[l])
          }else{
            index_temp <- which(df.qte.coverage$method == method_vec[i]
                                & df.qte.coverage$n == n_vec[j]
                                & df.qte.coverage$pn == pn_vec[l]
                                & df.qte.coverage$k == k_vec[k_temp])
          }
          
          df_temp <- df.qte.coverage[index_temp,]
          coverage_temp <- mean(df_temp$cover)
          
          df.coverage.temp <- data.frame(method=method_vec[i], coverage=coverage_temp,
                                    n=n_vec[j], pn=pn_vec[l], k=k_vec[k_temp],
                                    n_numeric=n_num_temp)
          df.coverage.final <- rbind(df.coverage.final, df.coverage.temp)
        }
      }
    }
  } # end for-loops
  
  # reorder methods in plot
  df.coverage.final$Method <- factor(df.coverage.final$method, levels=c("Extremal QTE CI", "Zhang's b out of n bootstrap"))
  df.coverage.final$pn <- factor(df.coverage.final$pn, levels=c("p[n]==5/n", "p[n]==1/n", "p[n]==5/nlogn"))
  
  coverage.plot <- ggplot(df.coverage.final, aes(x=k, y=coverage, group=method, col=method, linetype=Method)) + 
    geom_line(lwd=1) +
    facet_grid(pn ~ n, scales="free", labeller=label_parsed) + 
    scale_color_manual(values=coverage.colors) +
    scale_linetype_manual(values=c("solid", "dashed")) +
    geom_vline(data = df.coverage.final,
               aes(xintercept = n_numeric^(0.65)), linetype="dashed", size=0.3) +
    geom_hline(yintercept= coverage.target, linetype="solid") +
    theme_bw() + 
    ylab("Coverage") +
    scale_y_continuous(limits=c(0.4, 1)) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom") +
    guides(colour = guide_legend(override.aes = list(shape = NA), title="")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 16))
  
  return(coverage.plot)
}

