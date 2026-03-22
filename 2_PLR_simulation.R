rm(list=ls()) 

######Utils pacakges
library(foreach)
library(doParallel)
library(openxlsx)
library(doSNOW)
library(progress)
library(clusterGeneration)
library(mvtnorm)
library(tidyverse)
library(purrr)
library(openxlsx)

######ML packages
library(glmnet)
library(randomForest)
library(xgboost)
library(ranger)
library(gbm)
library(nnet)


dir_path <- "~/Causal/multi_causal_proj/codes"
set.seed(16) # = Seed for Replication = #


#souce helper functions
source("1_PLR_functions.R")


######################################
################run
######################################

#######multiple concurrent treatments

all_dat <- readRDS(paste0(dir_path, "/res/data500_PLR_202507_2As_inter_N1000.rds"))  



#set parameters
A_names <- c("A_1", "A_2", "A_inter")
folds_num = 5   #number of cross-fitting folds
split_num = 50  #number of splits


# machine learning parameters if needed
current_params_A = list(
  regression = list(),
  binary = list()  
)


all_res_list <- list()

# how many datasets to run
data_num_range = 1:500


#actuall run the DML
for(i in data_num_range){
  #for(i in 1:length(all_dat)){
  if(i %% 10 == 0)
    print(i)
  current_dat <- all_dat[[i]]
  current_dat$A_1 <- as.factor(current_dat$A_1)
  single_res <- get_all_method_res(A_names, current_dat, folds_num, split_num, params = list(), params_A = current_params_A)
  saveRDS(single_res, file = paste0(dir_path, "/res/PLR_sims/single_res_", i, ".rds"))
}

data_num = 500
all_res_list <- lapply(1:data_num, function(i) readRDS(paste0(dir_path, "/res/PLR_sims/single_res_", i, ".rds")))
all_res_list



######evaluate results: bias
true_values <- c(A_1 = 4, A_2 = 6, A_inter = 4)

lm_estimates <- extract_coefs(all_res_list, "lm")
lm_ses <- extract_ses(all_res_list, "lm")

rf_estimates <- extract_coefs(all_res_list, "rf")
rf_ses <- extract_ses(all_res_list, "rf")

gbm_estimates <- extract_coefs(all_res_list, "gbm")
gbm_ses <- extract_ses(all_res_list, "gbm")

nnet_estimates <- extract_coefs(all_res_list, "nnet")
nnet_ses <- extract_ses(all_res_list, "nnet")


# Compute metrics
lm_metrics <- calc_bias_rmse(lm_estimates, true_values)
rf_metrics <- calc_bias_rmse(rf_estimates, true_values)
gbm_metrics <- calc_bias_rmse(gbm_estimates, true_values)
nnet_metrics <- calc_bias_rmse(nnet_estimates, true_values)


temp_res <- data.frame(lm_metrics, rf_metrics, gbm_metrics, nnet_metrics )
temp_res


date_suffix <- format(Sys.Date(), "%m%d")
write.csv(temp_res, file = paste0(dir_path, "/res/2As_inter_results_folds", folds_num, "_N", N_num, "_", date_suffix, ".csv"))



#####evaluate results: variance and CP
estimate_list <- list(
  lm_estimates,
  rf_estimates,
  gbm_estimates,
  nnet_estimates
)
# Assign names
names(estimate_list) <- c("regression", "RF", "GBM", "nnet")

se_list <- list()
for (method in names(estimate_list)) {
  se_list[[method]] <- apply(estimate_list[[method]], 2, sd)
}
empirical_sd_df <- as.data.frame(do.call(cbind, se_list))
empirical_sd_df


lm_cis <- calc_confint(lm_estimates, lm_ses, 0.95, true_values)
rf_cis <- calc_confint(rf_estimates, rf_ses, 0.95, true_values)
gbm_cis <- calc_confint(gbm_estimates, gbm_ses, 0.95, true_values)
nnet_cis <- calc_confint(nnet_estimates, nnet_ses, 0.95, true_values)


analytical_sd_df <- as.data.frame(cbind(apply(lm_cis$ses, 2, mean),
apply(rf_cis$ses, 2, mean),
apply(gbm_cis$ses, 2, mean),
apply(nnet_cis$ses, 2, mean)))
names(analytical_sd_df) <- c("regression", "RF", "GBM", "nnet")


temp_res_CP <- cbind(lm_cis$coverage, rf_cis$coverage, gbm_cis$coverage, nnet_cis$coverage)
colnames(temp_res_CP) <- c("regression",  "RF", "GBM", "nnet")

temp_res_CP

#####save CIs output

output_list <- list(
  Empirical_SD = empirical_sd_df,
  Analytical_SD = analytical_sd_df,
  Coverage_Probability = as.data.frame(temp_res_CP)
)

date_suffix <- format(Sys.Date(), "%m%d")
write.xlsx(output_list, file =  paste0(dir_path, "/res/PLR_sims/2As_inter_CP_folds", folds_num, "_N", N_num, "_", date_suffix,".xlsx"), rowNames = TRUE, overwrite = TRUE)

