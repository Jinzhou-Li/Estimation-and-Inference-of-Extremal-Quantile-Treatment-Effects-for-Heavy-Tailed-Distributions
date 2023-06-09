library(evd) # for simulation
library(actuar) # for simulation
library(parallel)
library(ggpubr)

###########
source("Functions/calculate_true_qtes.R")
source("Functions/qte_extrapolation_hill.R")
source("Functions/qte_extrapolation_pickands.R")
source("Functions/qte_firpo_zhang.R")
source("Functions/qte_extrapolation_bootstrap.R")

source("Simulations/SimuFunc/Simu_func.R")
source("Simulations/SimuFunc/Plots_func.R")
source("Simulations/SimuFunc/Load_simu_results.R")