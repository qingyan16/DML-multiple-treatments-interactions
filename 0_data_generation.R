rm(list=ls()) 

######Utils pacakges
library(foreach)
library(doParallel)
library(openxlsx)
library(doSNOW)
library(progress)
library(clusterGeneration)
library(mvtnorm)
library(caret)
library(tidyverse)

######ML packages
#library(DoubleML)
#library(mlr3)
library(glmnet)
library(randomForest)
library(xgboost)
library(ranger)

dir_path <- "~/Library/CloudStorage/OneDrive-VUMC/Vandy/Research/Method/Causal/multi_causal_proj/codes/simulations"
set.seed(123)



gen_single_data <- function(N, all_dim, sigma, b, trt_type = "linear"){


  if (trt_type == "multi_regimen"){
    X_1 <- rnorm(N, mean = 0, sd = 1)
    X_2 <- rnorm(N, mean = 0, sd = 1)
    X_3 <- rnorm(N, mean = 0, sd = 1)
    X_4 <- rnorm(N, mean = 0, sd = 1)
    X_5 <- rnorm(N, mean = 0, sd = 1)
    X_6 <- rbinom(N, size = 1, prob = 0.1)
    X_7 <- rbinom(N, size = 1, prob = 0.3)
    X_8 <- rbinom(N, size = 1, prob = 0.5)
    X_9 <- rbinom(N, size = 1, prob = 0.7)
    X_10 <- rbinom(N, size = 1, prob = 0.9)
    
    X_all <- cbind(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10)
    
    ####################################################################
    ############################ Old codes##############################
    #g =  ifelse(X_1 < 0.5, -5, 5) + 6 * X_2 * X_3 + 4*X_4 + 4 * X_5 + 2 * X_6
    
    # eta <- cbind(
    #   1,   #R = 1
    #   1.5 * (X_1^2) - 0.5 * X_2 + 0.4 * (X_3 > 0) + 0.8 * X_4 * X_5 + 0.3 * X_6,   #R = 2
    #   5 * (X_1 * X_3) + 0.4 * (X_2 > -1) + 1.2 * X_4 - 0.7 * X_5,                  #R = 3
    #   2 * X_1^2  + 0.5 * X_2 * X_3 - 0.3 * X_4 + 1.5 * X_5 +  0.3 * X_6            #R = 4
    # )
    
    
    outcome_formula <- "ifelse(X_1 < 0, -5, 5) + ifelse(X_2 < 1, -8, 8) + 2*X_3 + 4*X_5 + 
                          1*X_6 + 2*X_7 + 4*X_9 + 5*X_10 +
                          4*X_3*X_4 + 6*X_5*X_10 +
                          6*X_5^2 + 4*X_9^2"
    
    g <- eval(parse(text = outcome_formula)) + rnorm(N)
    
    
    # eta_trt_assign_formula <- cbind(
    #   "1",
    #   "0.8*X_1^2 + 0.4*X_2^2 - 0.4*X_3 + 0.5*X_4      + 0.3*X_6 + 0.9*X_7*X_8 - 1.3* X_9 ", #R = 2
    #   "-1.2*X_1^2            + 1.8*X_3 + 2.5*(X_4>0)  + 0.3*X_6*X_7 - 1.2*X_8 + 0.5 * X_5 * X_10"   #R = 3                                                  #R = 4
    # )
    
    
    eta_trt_assign_formula <- cbind(
      "1",
      "0.8*X_1^2 + 0.4*X_2^2 - 0.4*X_3 + 0.7*X_4      + 0.3*X_6 + 0.9*X_7*X_8 - 1.3* X_9 ", #R = 2
      "-1.2*X_1^2            + 1.8*X_3 + 2.5*(X_4>0)  + 0.3*X_6*X_7 - 1.2*X_8 + 0.5 * X_5 * X_10"   #R = 3                                                  #R = 4
    )
    
    
    trt_nums = length(eta_trt_assign_formula)
    
    eta <- matrix(NA, nrow = N, ncol=trt_nums)
    for(i in 1:trt_nums){
      eta[,i] <- eval(parse(text = eta_trt_assign_formula[i]))
    }
    
    colMeans(eta)
    exp_eta <- exp(eta)
    probabilities <- exp_eta / (rowSums(exp_eta))  # Add reference class probability
    
    #probabilities <- matrix(rep(single_probabilities, N), nrow = N, byrow = TRUE)
    
    
    regimen <- apply(probabilities, 1, function(p) which(rmultinom(1, size = 1, prob = p) == 1))
    epsilon <- rnorm(N)
    
    
    # Y_1 = g
    # Y_2 = g + 5
    # Y_3 = g + 10
    # Y_4 = g + 15
    # 
    # Y = case_when(regimen == 2 ~ Y_2,
    #               regimen == 3 ~ Y_3,
    #               regimen == 4 ~ Y_4,
    #               TRUE ~ Y_1)
    # 
    # 
    # temp_dat <- cbind(Y, Y_1, Y_2, Y_3, Y_4, regimen)
    
    
    
    Y_1 = g
    Y_2 = g + 5
    Y_3 = g + 15 * X_9
    
    Y = case_when(regimen == 2 ~ Y_2,
                  regimen == 3 ~ Y_3,
                  TRUE ~ Y_1)
    
    
    temp_dat <- cbind(Y, Y_1, Y_2, Y_3, regimen)
    Y = Y + epsilon
    
    dat <- data.frame(Y, regimen, X_all)
  }
  
  if (trt_type == "2As_inter"){
    X_1 <- rnorm(N, mean = 0, sd = 1)
    X_2 <- rnorm(N, mean = 0, sd = 1)
    X_3 <- rnorm(N, mean = 0, sd = 1)
    X_4 <- rnorm(N, mean = 0, sd = 1)
    X_5 <- rnorm(N, mean = 0, sd = 1)
    X_6 <- rbinom(N, size = 1, prob = 0.1)
    X_7 <- rbinom(N, size = 1, prob = 0.3)
    X_8 <- rbinom(N, size = 1, prob = 0.5)
    X_9 <- rbinom(N, size = 1, prob = 0.7)
    X_10 <- rbinom(N, size = 1, prob = 0.9)
    
    X <- cbind(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10)
    
    
    # outcome_formula <- "ifelse(X_1 < 0, -2, 2) + ifelse(X_2 < 1, -1, 1) + 2*X_3 + 2*X_5 + 
    #                       1*X_6 + 1*X_7 - 2*X_9 - 0.5*X_10 +
    #                       2*X_3*X_4 + 2*X_5*X_10 +
    #                       2*X_5^2 + 2*X_9^2"
    
    outcome_formula <- "ifelse(X_1 < 0, -2, 2) + ifelse(X_2 < 1, -1, 1) + 2*X_3 + 2*X_5 + 
                          1*X_6 + X_7 - 2*X_9 - 0.5*X_10 +
                          2*X_3*X_4 + 2*X_5*X_10 +
                          2*X_5^2 + 2*X_9^2"
    
    g <- eval(parse(text = outcome_formula)) + rnorm(N)
    
    pA_1 <- plogis(1.3*X_1*X_2 + 0.7*X_2^2 - 0.4*X_3 +exp(X_4) + 1.5*X_7*X_9 - 1.5*X_10) 
    A_1 = rbinom(N, 1, pA_1)
    
    A_2 <- 1/(1+exp(X[,1])) - 1/(1+exp(X[,2])) + 0.5*X[, 3] + 0.25*((X[,5]>0) - X[,6]) + 0.1*(X[,7]+X[,9]*X[,10])
    #A_2 <- 1/(1+exp(X[,1])) - 1/(1+exp(X[,2])) + 0.5*X[, 3] + 1.5*((X[,5]>0) - (X[,6]>0)) + 0.8*(X[,7]+X[,9]*X[,10])
    A_2 <- A_2 + rnorm(N)
    
    
    Y = 4 * A_1 + 6 * A_2  + 4 * A_1 * A_2 + g + rnorm(N) # = Generate y = #
    dat <- data.frame(Y=Y, A_1=A_1, A_2=A_2, A_inter=A_1 * A_2, X)
  }
  
  
  ######################################
  ###### settings that are not used
  #######################################

  if (trt_type == "single_linear"){
    X_1 <- rnorm(N, mean = 0, sd = 1)
    X_2 <- rnorm(N, mean = 0, sd = 1)
    X_3 <- rnorm(N, mean = 0, sd = 1)
    X_4 <- rnorm(N, mean = 0, sd = 1)
    X_5 <- rnorm(N, mean = 0, sd = 1)
    X_6 <- rbinom(N, size = 1, prob = 0.1)
    X_7 <- rbinom(N, size = 1, prob = 0.3)
    X_8 <- rbinom(N, size = 1, prob = 0.5)
    X_9 <- rbinom(N, size = 1, prob = 0.7)
    X_10 <- rbinom(N, size = 1, prob = 0.9)
    
    X <- cbind(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10)
    
    
    outcome_formula <- "ifelse(X_1 < 0, -2, 2) + ifelse(X_2 < 1, -1, 1) + 2*X_3 + 2*X_5 + 
                          1*X_6 + 0.5*X_7 - 0.5*X_9 - 0.5*X_10 +
                          2*X_3*X_4 + 1*X_5*X_10 +
                          2*X_5^2 + 2*X_9^2"
    
    g <- eval(parse(text = outcome_formula)) + rnorm(N)
    
    #pA_1 <- plogis(0.8*X_1^2 + 0.4*X_2^2 - 0.4*X_3 - 0.8*X_4 + 0.2 * X_5  + 0.5*X_6*X_7 - 1.5* X_9 + 0.5*X_10)
    #A_1 = rbinom(N, 1, pA_1)
    
    A_2 <- 1/(1+exp(X[,1])) - 1/(1+exp(X[,2])) + 0.5*X[, 3] + 0.25*((X[,5]>0) - (X[,6]>0)) + 0.1*(X[,7]+X[,9]*X[,10])
    A_2 <- A_2 + rnorm(N)
    
    Y =  15 * A_2 +  g + rnorm(N) # = Generate y = #
    
    #Y = 1 * A_1 + 2 * A_2 + 1 * A_1 * A_2 +  g + rnorm(N) # = Generate y = #
    dat <- data.frame(Y=Y, A_2=A_2, X)
    #dat <- data.frame(Y=Y, A_1=A_1, A_2=A_2, A_inter = A_1*A_2, X)
  }
  if (trt_type == "single_binary"){
    X_1 <- rnorm(N, mean = 0, sd = 1)
    X_2 <- rnorm(N, mean = 0, sd = 1)
    X_3 <- rnorm(N, mean = 0, sd = 1)
    X_4 <- rnorm(N, mean = 0, sd = 1)
    X_5 <- rnorm(N, mean = 0, sd = 1)
    X_6 <- rbinom(N, size = 1, prob = 0.1)
    X_7 <- rbinom(N, size = 1, prob = 0.3)
    X_8 <- rbinom(N, size = 1, prob = 0.5)
    X_9 <- rbinom(N, size = 1, prob = 0.7)
    X_10 <- rbinom(N, size = 1, prob = 0.9)
    
    X <- cbind(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10)
    
    
    
    outcome_formula <- "ifelse(X_1 < 0, -2, 2) + ifelse(X_2 < 1, -1, 1) + 2*X_3 + 2*X_5 + 
                          1*X_6 + X_7 - 2*X_9 - 0.5*X_10 +
                          2*X_3*X_4 + 2*X_5*X_10 +
                          2*X_5^2 + 2*X_9^2"
    
    g <- eval(parse(text = outcome_formula)) + rnorm(N)
    
    pA_1 <- plogis(1.3*X_1*X_2 + 0.7*X_2^2 - 0.4*X_3 +exp(X_4) + 1.5*X_7*X_9 - 1.5*X_10) 
    A_1 = rbinom(N, 1, pA_1)
    
    Y =  5 * A_1 +  g + rnorm(N) # = Generate y = #
    
    #Y = 1 * A_1 + 2 * A_2 + 1 * A_1 * A_2 +  g + rnorm(N) # = Generate y = #
    dat <- data.frame(Y=Y, A_1=A_1, X)
    #dat <- data.frame(Y=Y, A_1=A_1, A_2=A_2, A_inter = A_1*A_2, X)
  }
  
  
  
  return (dat)
}




#generate data of multiple treatments for PLR model (with interactions)
data_nums = 500
N = 1000

all_data_sets_PLR <- list()
for (i in 1:data_nums){
  if(i %%10 == 0)
    print(i)
  dat <- gen_single_data(N, all_dim, sigma, b, trt_type = "2As_inter")
  all_data_sets_PLR[[i]] <- dat
}


saveRDS(all_data_sets_PLR, file = paste0(dir_path, "/res/data", data_nums, "_PLR_202507_2As_inter_N", N, ".rds"))





###generate data of multi-valued regimens for IRM model 
data_nums = 100
N = 1000

all_data_sets_multi_binary <- list()
for (i in 1:data_nums){
  if(i %%10 == 0)
    print(i)
  dat <- gen_single_data(N, all_dim, sigma, b, trt_type = "multi_regimen")
  all_data_sets_multi_binary[[i]] <- dat
}

saveRDS(all_data_sets_multi_binary, file = paste0(dir_path, "/res/all_data_sets_multi_regimen_2025_3Rs.rds"))




