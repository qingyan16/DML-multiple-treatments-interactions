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
library(PSweight)
library(reshape2)


######ML packages
library(DoubleML)
library(mlr3)
library(glmnet)
library(randomForest)
library(xgboost)
library(nnet)
library(ranger)


library(reticulate)
use_virtualenv("/Users/qxiang/Library/CloudStorage/OneDrive-VUMC/Vandy/Research/Method/Single_Cell/DML_single_cell/codes/learning/python/py3.11/.venv")
use_python("/Users/qxiang/Library/CloudStorage/OneDrive-VUMC/Vandy/Research/Method/Single_Cell/DML_single_cell/codes/learning/python/py3.11/.venv/bin/python", required = TRUE)
dir_path <- "~/Library/CloudStorage/OneDrive-VUMC/Vandy/Research/Method/Causal/multi_causal_proj/codes/simulations"

set.seed(123) # Seed for Replication


####################################
######helper functions
###################################
sklearn_linear_model <- import("sklearn.linear_model")
sklearn_ensemble <- import("sklearn.ensemble")
dml <- import("doubleml")

# Access specific classes or functions
#LogisticRegression <- sklearn_linear_model$LogisticRegression
#LinearRegression <- sklearn_linear_model$LinearRegression
#RandomForestClassifier <- sklearn_ensemble$RandomForestClassifier
#RandomForestRegressor <- sklearn_ensemble$RandomForestRegressor
GradientBoostingClassifier <- sklearn_ensemble$GradientBoostingClassifier
GradientBoostingRegressor <- sklearn_ensemble$GradientBoostingRegressor

##function to run DML for multi-valued regimens
get_IRM_res <- function(ml_g, ml_m, data_num, treatment_levels, n_folds){
  all_coefs <- c()
  
  for (i in 1:data_nums){
    if(i %% 10 == 0)
      print(i)
    current_dat <- all_data_sets_multi_binary[[i]]
    dml_data <- dml$DoubleMLData(current_dat, 'Y', 'regimen')
    
    dml_obj <- dml$DoubleMLAPOS(
      dml_data,
      ml_g,
      ml_m,
      treatment_levels=treatment_levels,
      n_folds=n_folds,
      trimming_threshold=0.02
    )
    
    dml_fit <- dml_obj$fit()
    
    #get the resulted coefficients
    dml_fit$summary[['coef']]
    all_coefs <- rbind(all_coefs, as.vector(dml_fit$summary[['coef']]))
  }
  
  return (all_coefs)
}


##function to traditional methods (PS-based)
tradition_PSmethod <- function(all_data_sets, data_nums, weight, ps.method = "glm", trt_num){
  
  # Identify the columns with "X"
  x_columns <- colnames(all_data_sets[[1]])[grepl("^X", colnames(all_data_sets[[1]]))]
  
  # Create the formula
  ps.formula <- as.formula(paste("regimen ~", paste(x_columns, collapse = " + ")))
  
  all_res <- c()
  
  for (i in 1:data_nums){
    if(i %% 20 == 0)
      print(i)
    
    temp_dat <- all_data_sets[[i]]
    drop_cols <- grep("A", colnames(temp_dat), value = TRUE)
    temp_dat <- temp_dat %>% select(-all_of(drop_cols))
    temp_dat$regimen <- as.factor(temp_dat$regimen)
    
    ato1<-PSweight(ps.formula = ps.formula, yname = 'Y', data = temp_dat, weight = weight, ps.method = ps.method)
    a <- summary(ato1)
    
    if (trt_num == 8){
      all_res <- rbind(all_res, c(a$estimates[4, 1], a$estimates[7, 1]))
    }
    if (trt_num == 4){
      all_res <- rbind(all_res, a$estimates[1:3, 1])
    }
    if(trt_num == 3){
      all_res <- rbind(all_res, a$estimates[1:3, 1])
    }
    if(trt_num == 2){
      all_res <- rbind(all_res, a$estimates[1])
    }
  }
  
  return(all_res)
}


###results summary function for DML
ATE_est <- function(sim_coef, truth, trt_num){
  if(trt_num == 8){
    all_ATE1 <- sim_coef[, 5] - sim_coef[, 1]
    all_ATE2 <- sim_coef[, 8] - sim_coef[, 1]
    all_ATE_dat <- data.frame(all_ATE1, all_ATE2)
    
    mat_subtract_truth <- sweep(as.matrix(all_ATE_dat), 2, truth, '-')
    bias <- colMeans(mat_subtract_truth) 
    rMSE <- sqrt(colMeans(mat_subtract_truth**2))
  }
  if(trt_num == 4){
    all_ATE1 <- sim_coef[, 2] - sim_coef[, 1]
    all_ATE2 <- sim_coef[, 3] - sim_coef[, 1]
    all_ATE3 <- sim_coef[, 4] - sim_coef[, 1]
    all_ATE_dat <- data.frame(all_ATE1, all_ATE2, all_ATE3)
    
    mat_subtract_truth <- sweep(as.matrix(all_ATE_dat), 2, truth, '-')
    bias <- colMeans(mat_subtract_truth) 
    rMSE <- sqrt(colMeans(mat_subtract_truth**2))
  }
  if(trt_num == 3){
    all_ATE1 <- sim_coef[, 2] - sim_coef[, 1]
    all_ATE2 <- sim_coef[, 3] - sim_coef[, 1]
    all_ATE3 <- sim_coef[, 3] - sim_coef[, 2]
    all_ATE_dat <- data.frame(all_ATE1, all_ATE2, all_ATE3)
    
    mat_subtract_truth <- sweep(as.matrix(all_ATE_dat), 2, truth, '-')
    bias <- colMeans(mat_subtract_truth) 
    rMSE <- sqrt(colMeans(mat_subtract_truth**2))
  }
  if(trt_num == 2){
    all_ATE_dat <- sim_coef[, 2] - sim_coef[, 1]
    
    bias <- mean(all_ATE_dat) - truth
    rMSE <- sqrt(mean((all_ATE_dat - truth)**2))
  }
  
  return(list(bias = bias, rMSE = rMSE, relative_bias = sweep(mat_subtract_truth, 2, truth, `/`)))
}


###results summary for traditional methods
ATE_summaries_tradition <- function(all_ATE_dat, truth){
  mat_subtract_truth <- sweep(as.matrix(all_ATE_dat), 2, truth, '-')
  bias <- colMeans(mat_subtract_truth) 
  rMSE <- sqrt(colMeans(mat_subtract_truth**2))
  
  return (list(bias = bias, rMSE = rMSE, relative_bias = sweep(mat_subtract_truth, 2, truth, `/`)))
}




####################################
######DoubleML results
###################################

#######read data with 3 trts
all_data_sets_multi_binary <- readRDS(paste0(dir_path, "/res/all_data_sets_multi_regimen_2025_3Rs.rds"))
data_nums = length(all_data_sets_multi_binary)

#specify learners
ml_g <- GradientBoostingRegressor()
ml_m <- GradientBoostingClassifier()
data_nums = 100
treatment_levels <- c(1,2,3)

#
all_coefs_IRM_gbm <- get_IRM_res(ml_g, ml_m, data_num, treatment_levels = c(1,2,3), n_folds = 2L)
res_DML_gbm_contrast_2folds <- ATE_est(all_coefs_IRM_gbm, truth = c(5, 10.5, 5.5), trt_num = 3)$relative_bias
colMeans(res_DML_gbm_contrast_2folds)


all_coefs_IRM_gbm <- get_IRM_res(ml_g, ml_m, data_num, treatment_levels = c(1,2,3), n_folds = 5L)
res_DML_gbm_contrast_5folds <- ATE_est(all_coefs_IRM_gbm, truth = c(5, 10.5, 5.5), trt_num = 3)$relative_bias
colMeans(res_DML_gbm_contrast_5folds)



####################################
######other methods (PS) results
###################################

all_data_sets_multi_binary <- readRDS(paste0(dir_path, "/res/all_data_sets_multi_regimen_2025_3Rs.rds"))
data_nums = length(all_data_sets_multi_binary)

#########glm propensity score
res_IPW_glm <- tradition_PSmethod(all_data_sets_multi_binary, data_nums,  "IPW", ps.method = "glm", trt_num = 3)
res_ATM_glm <- tradition_PSmethod(all_data_sets_multi_binary, data_nums,  "matching", ps.method = "glm", trt_num = 3)

res_IPW_glm_contrast <- ATE_summaries_tradition(res_IPW_glm, c(5, 10.5, 5.5))$relative_bias
res_ATM_glm_contrast <- ATE_summaries_tradition(res_ATM_glm, c(5, 10.5, 5.5))$relative_bias

colMeans(res_IPW_glm_contrast)
colMeans(res_ATM_glm_contrast)


########gbm propensity score
###gbm propensity score
res_IPW_gbm <- tradition_PSmethod(all_data_sets_multi_binary, data_nums,  "IPW", ps.method = "gbm", trt_num = 3)
res_ATM_gbm <- tradition_PSmethod(all_data_sets_multi_binary, data_nums,  "matching", ps.method = "gbm", trt_num = 3)

res_IPW_gbm_contrast <- ATE_summaries_tradition(res_IPW_gbm, c(5, 10.5, 5.5))$relative_bias
res_ATM_gbm_contrast <- ATE_summaries_tradition(res_ATM_gbm, c(5, 10.5, 5.5))$relative_bias

colMeans(res_IPW_gbm_contrast)
colMeans(res_ATM_gbm_contrast)




###########################
###Results reproting
############################

data <- data.frame(
  Bias = c(res_DML_gbm_contrast_2folds, res_DML_gbm_contrast_5folds, res_IPW_glm_contrast, res_IPW_gbm_contrast, res_ATM_glm_contrast, res_ATM_gbm_contrast),
  Method = factor(rep(c("DML GBM 2-fold", "DML GBM 5-fold", "IPW GLM", "IPW GBM", "PSM GLM", "PSM GBM"), each = nrow(res_DML_gbm_contrast_2folds) * ncol(res_DML_gbm_contrast_2folds))),
  Column = factor(rep(rep(1:3, each = nrow(res_DML_gbm_contrast_2folds)), times = 6))
)



# Highlight the first boxplot ("DML GBM")
p <- ggplot(data, aes(x = Method, y = Bias, fill = (Method == "DML GBM 2-fold" | Method ==  "DML GBM 5-fold"))) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Dashed line at y = 0
  facet_wrap(~ Column, scales = "fixed", nrow = 1,
             labeller = labeller(Column = c("1" = "ATE 2-1", "2" = "ATE 3-1", "3" = "ATE 3-2")),
             strip.position = "top") +
  theme_minimal() +
  labs(x = "", y = "Relative bias") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    strip.text = element_text(size = 12, face = "bold"), # Bold facet titles
    plot.title = element_text(size = 16, face = "bold"), # Bold plot title
    panel.spacing = unit(2, "lines") # Add more space between facets
  ) +
  scale_x_discrete(limits = c("DML GBM 2-fold", "DML GBM 5-fold", "IPW GLM", "IPW GBM", "PSM GLM", "PSM GBM")) + # Preserve order
  scale_y_continuous(breaks = seq(-0.8, 0.8, 0.1),
                     labels = c("-80%", "-70%", "-60%", "-50%", "-40%", "-30%", "-20%", "-10%", "0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%"),
                     limits = c(-0.8, 0.8)) +
  scale_fill_manual(values = c("TRUE" = "grey80", "FALSE" = "white"), guide = "none") # Highlight "DML GBM" in red

print(p)
ggsave(filename = paste0(dir_path, "/res/plot_IRM_all_covariates.pdf"), plot = p, width = 8, height = 6)














