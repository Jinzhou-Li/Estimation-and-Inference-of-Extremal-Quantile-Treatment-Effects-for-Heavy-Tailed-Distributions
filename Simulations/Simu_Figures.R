source("LoadFunctions.R")

###
pnstrings <- c("pn5", "pn1", "pn5nlogn")
pntext <- c("5/n", "1/n", "5/nlogn")
nlist <- c(1000, 2000, 5000)

### Simulate true qtes for all n and pn: take a few minutes
pnmatrix <- cbind(5/nlist, 1/nlist, 5/nlist/log(nlist))
true_qte_mat1 <- calculate_true_qtes(h1.fun, pnmatrix, Nsim = 1e8)
true_qte_mat2 <- calculate_true_qtes(h2.fun, pnmatrix, Nsim = 1e8)
true_qte_mat3 <- calculate_true_qtes(h3.fun, pnmatrix, Nsim = 1e8)

### Get data.frame for different models
dataframe_list_h1 <- Load_simu_results(folderstring="H1/", pnstrings, pntext, nlist, true_qte_mat1)
df.qte.mse.h1 <- dataframe_list_h1[[1]]
df.qte.coverage.h1 <- dataframe_list_h1[[2]]

dataframe_list_h2 <- Load_simu_results(folderstring="H2/", pnstrings, pntext, nlist, true_qte_mat2)
df.qte.mse.h2 <- dataframe_list_h2[[1]]
df.qte.coverage.h2 <- dataframe_list_h2[[2]]

dataframe_list_h3 <- Load_simu_results(folderstring="H3/", pnstrings, pntext, nlist, true_qte_mat3)
df.qte.mse.h3 <- dataframe_list_h3[[1]]
df.qte.coverage.h3 <- dataframe_list_h3[[2]]

############ MSE plots ############
mse.colors <- c("#52854C", "#D55E00", "#4E84C4")
h1_plot_mse <- squared_error_boxplot(df.qte.mse.h1, mse.colors)+ ggtitle(expression("H"[1])) 
h2_plot_mse <- squared_error_boxplot(df.qte.mse.h2, mse.colors)+ ggtitle(expression("H"[2])) 
h3_plot_mse <- squared_error_boxplot(df.qte.mse.h3, mse.colors)+ ggtitle(expression("H"[3])) 

h2_plot_mse <- h2_plot_mse + labs(y = NULL)
h3_plot_mse <- h3_plot_mse + labs(y = NULL)
plot_comb_mse <- ggarrange(h1_plot_mse, h2_plot_mse, h3_plot_mse,
                            nrow = 1, ncol = 3,
                           common.legend=TRUE, legend="bottom") #+ bgcolor("White")
ggsave("Simulations/Figures/Figure2.png", plot=plot_comb_mse, width = 18, height = 8)

############ Coverage plots: main ############
coverage.colors <- c("#52854C", "#D55E00", "#D55E00", "#4E84C4","#4E84C4","#4E84C4")
h1_plot_cover <- coverage_plot_func(df.qte.coverage.h1, coverage.colors, 
                                    coverage.target=0.9, CI.level.simu=0.95)[[1]] + ggtitle(expression("H"[1])) 
h2_plot_cover <- coverage_plot_func(df.qte.coverage.h2, coverage.colors, 
                                    coverage.target=0.9, CI.level.simu=0.95)[[1]] + ggtitle(expression("H"[2])) 
h3_plot_cover <- coverage_plot_func(df.qte.coverage.h3, coverage.colors, 
                                    coverage.target=0.9, CI.level.simu=0.95)[[1]] + ggtitle(expression("H"[3])) 

plot_comb_cover <- ggarrange(h1_plot_cover, h2_plot_cover, h3_plot_cover,
                           nrow = 1, ncol = 3,
                           common.legend=TRUE, legend="none") #+ bgcolor("White")
ggsave("Simulations/Figures/Figure3.png", plot=plot_comb_cover, width = 19, height = 8.8)

############ QQ plots ############
h2_qqplot_list <- qqplot_func(df.qte.mse.h2, nlist, pnstrings)

# only present the following one in the paper to save space
h2.qqplot <- h2_qqplot_list[[9]] + ggtitle(expression("Q-Q Plots of Estimated QTEs for Model H"[2],", n=5000, p"["n"]*"=5/nlog(n)"))
ggsave("Simulations/Figures/Figure5.pdf", h2.qqplot, width=11, height=4)


######################################################################
############################ Plots in appendix #######################
######################################################################

### Get data.frame for different models
dataframe_list_h1 <- Load_simu_results_supp(folderstring="H1/", pnstrings, pntext, nlist, true_qte_mat1)
df.qte.mse.h1 <- dataframe_list_h1[[1]]
df.qte.coverage.h1 <- dataframe_list_h1[[2]]

dataframe_list_h2 <- Load_simu_results_supp(folderstring="H2/", pnstrings, pntext, nlist, true_qte_mat2)
df.qte.mse.h2 <- dataframe_list_h2[[1]]
df.qte.coverage.h2 <- dataframe_list_h2[[2]]

dataframe_list_h3 <- Load_simu_results_supp(folderstring="H3/", pnstrings, pntext, nlist, true_qte_mat3)
df.qte.mse.h3 <- dataframe_list_h3[[1]]
df.qte.coverage.h3 <- dataframe_list_h3[[2]]

############ MSE plots ############
mse.colors <- c("#52854C", "#D55E00", "#4E84C4")
h1_plot_mse_supp <- squared_error_plot_supp(df.qte.mse.h1, mse.colors)+ ggtitle(expression("H"[1])) 
h2_plot_mse_supp <- squared_error_plot_supp(df.qte.mse.h2, mse.colors)+ ggtitle(expression("H"[2])) 
h3_plot_mse_supp <- squared_error_plot_supp(df.qte.mse.h3, mse.colors)+ ggtitle(expression("H"[3])) 

h2_plot_mse_supp <- h2_plot_mse_supp + labs(y = NULL)
h3_plot_mse_supp <- h3_plot_mse_supp + labs(y = NULL)
plot_comb_mse_supp <- ggarrange(h1_plot_mse_supp, h2_plot_mse_supp, h3_plot_mse_supp,
                                nrow = 1, ncol = 3,
                                common.legend=TRUE, legend="bottom") # + bgcolor("White")
ggsave("Simulations/Figures/Figure6.png", plot=plot_comb_mse_supp, width = 15, height = 6)

############ Coverage plot ############
coverage.colors <- c("#4E84C4","#52854C")
h1_plot_cover_supp <- coverage_plot_supp(df.qte.coverage.h1, coverage.colors, 
                                         coverage.target=0.9) + ggtitle(expression("H"[1])) 
h2_plot_cover_supp <- coverage_plot_supp(df.qte.coverage.h2, coverage.colors, 
                                         coverage.target=0.9) + ggtitle(expression("H"[2])) 
h3_plot_cover_supp <- coverage_plot_supp(df.qte.coverage.h3, coverage.colors, 
                                         coverage.target=0.9) + ggtitle(expression("H"[3])) 

h2_plot_cover_supp <- h2_plot_cover_supp + labs(y = NULL)
h3_plot_cover_supp <- h3_plot_cover_supp + labs(y = NULL)
plot_comb_cover_plot_cover_supp <- ggarrange(h1_plot_cover_supp, h2_plot_cover_supp, h3_plot_cover_supp,
                                             nrow = 1, ncol = 3,
                                             common.legend=TRUE, legend="none") #+ bgcolor("White")
ggsave("Simulations/Figures/Figure7.png", plot=plot_comb_cover_plot_cover_supp, width = 15, height = 6)
