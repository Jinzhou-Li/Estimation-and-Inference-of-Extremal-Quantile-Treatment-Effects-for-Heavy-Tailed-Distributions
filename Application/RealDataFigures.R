# The code in this file is used to estimate (extremal) QTEs of college education on wage.
library(haven)
library(gridExtra)
library(boot)
library(qte)
library(kableExtra)
library(rJava)
library(glmulti)
library(ggpubr)
source("Functions/qte_extrapolation_hill.R")
source("Functions/qte_firpo_zhang.R")
source("Application/real_result_func.R")

###########################################################################
#####                       Preprocess data                          ######
###########################################################################
df.raw <- read_dta("Application/hhv-neo_v3.dta")
df.raw <- df.raw[!is.na(df.raw$wage) & (df.raw$collegedeg | df.raw$somecollege | df.raw$hsgrad),]       # select only people graduated from high-school
df.nan <- df.raw[,c("wage", "black", "hisp", "south", "west", "northeast", "curban", "broken", 
                                    "age80", "mhgc_mi", "fhgc_mi", "faminc79_th", "numsibs", 
                                    "sasvab1", "sasvab2", "sasvab3", "sasvab4", "sasvab5", "sasvab6" ,
                                    "sgr9_lang_gpa", "sgr9_math_gpa", "sgr9_scosci_gpa", "sgr9_sci_gpa"
                                    )]  
df.nan <- cbind(college = (df.raw$collegedeg | df.raw$somecollege), df.nan)                             # college is 1 if participant went to college for some time
df <- df.nan[complete.cases(df.nan),]                                                                   # remove missing values

### Generate categorical variables 'race' and 'region of residence' based on 
# ('black', 'hisp') and ('south', 'west', 'northeast'), respectively.

race <- rep("OtherRace", nrow(df))
race[df[,3]==1] <- "black"
race[df[,4]==1] <- "hisp"
df <- data.frame(df, race=as.factor(race))

region <- rep("OtherRegion", nrow(df))
region[df[,5]==1] <- "south"
region[df[,6]==1] <- "west"
region[df[,7]==1] <- "northeast"
df <- data.frame(df, region=as.factor(region))
df <- df[,-c(3:7)] # final data frame with 19 covariates

###########################################################################
#####                 Extremal QTE estimation                        ######
###########################################################################

## Target QTES: 0.99, 0.995, 0.997, 0.999
targets <- c(0.99, 0.995, 0.997, 0.999)
pn <- 1 - targets

### We try two versions of the GLM to estimate propensity scores as described in the paper
# 1. Fit a logistic regression model with all linear terms
# 2. Using AIC, first do a model selection based on linear terms, then selecting second order terms based on selected variables

formula.linear <- as.formula(paste0("college~",paste(names(df)[-c(1,2)], collapse = '+')))

# glm.select.first.order <- glmulti(formula.linear, data=df,intercept = TRUE,
#                       level=1, method="g", crit="aic")
#formula.first.order <- as.formula("college~1+race+region+age80+mhgc_mi+fhgc_mi+sasvab2+sasvab5+sasvab6+sgr9_lang_gpa+sgr9_scosci_gpa")
# glm.select.second.order <- glmulti(formula.first.order, data=df,intercept = TRUE,
#                       level=2, method="g", crit="aic")

# the following formula results from the above codes: we saved it to save time
formula.second.order <- as.formula("college~1+race+region+age80+mhgc_mi+fhgc_mi+sasvab5+sasvab6+sgr9_scosci_gpa+sasvab2:mhgc_mi+sasvab5:age80+sasvab5:sasvab2+sasvab6:sasvab2+sgr9_lang_gpa:age80+sgr9_lang_gpa:mhgc_mi+sgr9_scosci_gpa:age80+race:age80+race:mhgc_mi+region:age80+region:fhgc_mi+region:sgr9_scosci_gpa")

############ implement methods
res_PROP1 <- real_data_results(formula.linear, df, pn)
res_PROP2 <- real_data_results(formula.second.order, df, pn)

###################### Plots
Plot_list <- list()
x_axis <- rep(paste(targets, sep=""), each=4)
method <- rep(c("Firpo-Zhang (PROP1)", "Extremal QTE (PROP1)",
                "Firpo-Zhang (PROP2)", "Extremal QTE (PROP2)"), 
              length.out=length(x_axis))
point_estimate <- CI_low <- CI_up <- gamma0 <- gamma1 <- c()
for (pn_index in 1:length(pn)) {
  res_hill1 <- res_PROP1[[1]][[pn_index]]
  res_hill2 <- res_PROP2[[1]][[pn_index]]
  res_zhang1 <- res_PROP1[[2]][[pn_index]]
  res_zhang2 <- res_PROP2[[2]][[pn_index]]
  res_hill1$PROP <- "PROP1"
  res_hill2$PROP <- "PROP2"
  res_hill1$k <- 1:300
  res_hill2$k <- 1:300
  resdf.combo <- rbind(res_hill1, res_hill2)
  
  # For Gamma_Hill plot
  df.gamma <- data.frame(gamma=c(resdf.combo$gamma0, resdf.combo$gamma1), 
                         k=c(1:300,1:300,1:300,1:300), 
                         PROP=c(rep("PROP1",300), rep("PROP2",300), rep("PROP1",300), rep("PROP2",300)), 
                         Gamma=c(rep("gamma0",600), rep("gamma1",600)))
  
  ### To get results for the main plot
  res_k_PROP1 <- resdf.combo[which(resdf.combo$PROP=="PROP1" & resdf.combo$k==85),]
  res_k_PROP2 <- resdf.combo[which(resdf.combo$PROP=="PROP2" & resdf.combo$k==85),]
  
  point_estimate <- c(point_estimate, res_zhang1$qtes, res_k_PROP1$qtes,
                      res_zhang2$qtes, res_k_PROP2$qtes)
  CI_low <- c(CI_low, res_zhang1$CIlower, res_k_PROP1$CIlower,
                      res_zhang2$CIlower, res_k_PROP2$CIlower)
  CI_up <- c(CI_up, res_zhang1$CIupper, res_k_PROP1$CIupper,
                      res_zhang2$CIupper, res_k_PROP2$CIupper)
  
  ### The plots in appendix
  Plot_list[[pn_index]] <-  ggplot(resdf.combo, aes(x=k, y=qtes, colour=PROP)) + geom_line(linewidth=1) +
    geom_ribbon(aes(ymin=CIlower, ymax=CIupper, fill = PROP), linetype=2, alpha=0.3) +
    geom_hline(yintercept = 0) +
    ylab("QTE [$ / hour]") + 
    ggtitle(paste0("Estimates of ",1-pn[pn_index],"-QTE")) +
    theme_bw() +
    guides(color=guide_legend(title=""), fill=guide_legend(title="")) +
    scale_x_continuous(breaks=seq(0,300,50))
  
  # same for different pn_index
  Plot_gamma <-  ggplot(df.gamma, aes(x=k, y=gamma, 
                                      group = interaction(PROP, Gamma), 
                                      colour=PROP, shape=Gamma)) +
    geom_point(size=2) +
    geom_line(linewidth=0.5) +
    geom_hline(yintercept = 0) +
    ylab("") + 
    ggtitle("Estimated EVIs") +
    theme_bw() +
    guides(color=guide_legend(title=""), fill=guide_legend(title=""), shape=guide_legend(title="")) +
    scale_x_continuous(breaks=seq(0,300,50)) +
    theme(text = element_text(size = 15))
  
}

df_gg <- data.frame(x_axis=x_axis, method=method, point_estimate=point_estimate, CI_low=CI_low, CI_up=CI_up)
df_gg$method <- factor(df_gg$method, levels=c("Firpo-Zhang (PROP1)", "Extremal QTE (PROP1)",
                                                "Firpo-Zhang (PROP2)", "Extremal QTE (PROP2)"))
Method_color <- c("Firpo-Zhang (PROP1)"="#52854C", "Extremal QTE (PROP1)"="#4E84C4",
                  "Firpo-Zhang (PROP2)"="#52854C", "Extremal QTE (PROP2)"="#4E84C4")
Method_line <- c("Firpo-Zhang (PROP1)"="solid", "Extremal QTE (PROP1)"="solid",
                 "Firpo-Zhang (PROP2)"="dotted", "Extremal QTE (PROP2)"="dotted")

### Main plot
plot_real <- ggplot(df_gg, aes(x=x_axis, y=point_estimate, color=method)) + 
  geom_point(size=2, position=position_dodge(width=0.4)) +
  scale_color_manual(values=Method_color) +
  scale_linetype_manual(values=Method_line) +
  labs(x = 'Probability level', y = 'Estimated QTE value') +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up, width = 0.35, linetype=method), size=0.7, 
                position=position_dodge(width=0.4)) +
  ylim(-50,420) +
  theme_bw(base_size = 20)+
  theme(legend.title= element_blank(),
        legend.justification = c(0, 1), legend.position = c(0, 1))

ggsave("Application/Figures/Figure4.png", plot = plot_real, width = 15, height = 15, 
       units = "cm", bg="white")

### Plots in the appendix
ggsave(plot=Plot_gamma, width=8, height=5, filename="Application/Figures/Figure8.pdf")

comboplot <- ggarrange( Plot_list[[1]], Plot_list[[2]],
                        Plot_list[[3]],Plot_list[[4]],
                        ncol = 2, nrow = 2, 
                        common.legend = TRUE, legend="bottom")
ggsave(plot=comboplot, width=8, height=5, filename="Application/Figures/Figure9.pdf")