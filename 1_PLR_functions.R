
########################################
###########DML functions
###########################################[]

#function to generate y_residuals
gen_y_residuals <- function(A_names = A_names, train_data, test_data, method_name, params = list()){
  #model_response <- caret::train(Y ~ . - A_1 - A_2 - A_inter, data = train_data, method = method_name, verbosity = 0)
  formula <- as.formula(paste("Y ~ . - ", paste(A_names, collapse = " - ")))
  
  if (method_name == "lasso"){
    y <- train_data$Y
    x <- model.matrix(formula, train_data)[, -1]  
    
    model_response <-  cv.glmnet(x, y, alpha = 1)
    predictions <- predict(model_response, newx = model.matrix(formula, test_data)[, -1], s = "lambda.min")
    residuals <- test_data$Y - predictions
  }
  if (method_name == "regression"){     #In our simulations we only consider continuous outcomes Y
    model_response <- lm(formula, data = train_data)
    predictions <- predict(model_response, newdata = test_data)
    residuals <- test_data$Y - predictions
  }
  if (method_name == "rf"){
    model_response <- ranger(formula, data = train_data)
    predictions <- predict(model_response, data = test_data)$predictions
    residuals <- test_data$Y - predictions
  }
  if (method_name == "gbm"){
    
    final_params <- modifyList(
      list(
        formula = formula,
        distribution = "gaussian",
        data = train_data,
        
        shrinkage = 0.01, 
        n.trees = 300, 
        interaction.depth = 5, 
        n.minobsinnode = 10
      ),
      params
    )
    
    # Execute GBM with parameters
    model_response <- do.call(gbm::gbm, final_params)
    
    predictions <- predict(model_response, newdata = test_data)
    residuals <- test_data$Y - predictions
  }
  if(method_name == "nnet"){
    x_train <- model.matrix(formula, train_data)[, -1]
    y_train <- train_data$Y
    
    # Train the neural network with scaled features and linear output
    final_params <- modifyList(
      list(
        x = x_train,
        y = y_train,
        size = 16,
        decay = 0.01,
        maxit = 200,
        linout = TRUE,
        trace = FALSE
      ),
      params  # override any from tuning grid
    )
    
    model_response <- do.call(nnet::nnet, final_params)
    
    # Scale the test data inputs
    x_test <- model.matrix(formula, test_data)[, -1]
    #x_test[, cols_to_scale] <- x_test_scaled_cols <- scale(x_test[,cols_to_scale], center = attr(x_train_scaled_cols, "scaled:center"), scale = attr(x_train_scaled_cols, "scaled:scale"))
    #x_test_scaled <- x_test
    
    # Predict using the neural network
    predictions <- predict(model_response, newdata = x_test)
    residuals <- test_data$Y - predictions
    
  }
  return (residuals)
}




#function to generate treatment_residuals
gen_treatment_residuals <- function(A_names, train_data, test_data, method_name, params = list()){
  # Loop through each treatment
  residuals_list <- list()
  
  #A = "A_1"
  
  for (A in A_names){
    other_terms <- paste(setdiff(A_names, A), collapse = " - ")
    formula_str <- if (nchar(other_terms) > 0) {
      paste(A, "~ . - Y -", other_terms)
    } else {
      paste(A, "~ . - Y")
    }
    
    #if(length(unique(train_data$A)) < 4){
    #  train_data[[A]] <- factor(train_data[[A]])
    #  test_data[[A]] <- factor(test_data[[A]])
    #}      
    
    # Convert to formula
    formula <- as.formula(formula_str)
    
    if (method_name == "lasso"){
      y <- train_data[[A]]
      x <- model.matrix(formula, train_data)[, -1]  
      if(is.factor(test_data[[A]])){
        model_A <-  cv.glmnet(x, y, alpha = 1, family = "binomial")
        predictions <- predict(model_A, newx = model.matrix(formula, test_data)[, -1], s = "lambda.min", type = "response")
        residuals <- as.numeric(test_data[[A]]) - 1 - predictions
      }else{
        model_A <-  cv.glmnet(x, y, alpha = 1)
        predictions <- predict(model_A, newx = model.matrix(formula, test_data)[, -1], s = "lambda.min")
        residuals <- test_data[[A]] - predictions
      }
      residuals_list[[A]] <- residuals
    }
    if (method_name == "regression"){
      if(is.factor(test_data[[A]])){
        model_A <- glm(formula, data = train_data, family = "binomial")
        predictions <- predict(model_A, newdata = test_data, type = "response")
        residuals <- as.numeric(test_data[[A]]) - 1 - predictions
      }else{
        model_A <- lm(formula, data = train_data)
        predictions <- predict(model_A, newdata = test_data)
        residuals <- test_data[[A]] - predictions
      }
      residuals_list[[A]] <- residuals
    }
    if (method_name == "rf"){
      if(is.factor(test_data[[A]])){
        #model_A <- ranger(formula, data = train_data, probability = TRUE, mtry = 2, num.trees = 1000)
        model_A <- ranger(formula, data = train_data, probability = TRUE)
        residuals <- as.numeric(test_data[[A]]) - 1 - predict(model_A, data = test_data)$predictions[,2]
      }
      else{
        final_params <- modifyList(
          #For n = 1000
          list(
            formula,
            data = train_data,
            mtry = 4,
            min.node.size = 10
          ),
          #list(
          #  formula,
          #  data = train_data
          #),
          params  # override any from tuning grid
        )
        model_A <- do.call(ranger::ranger, final_params)
        residuals <- test_data[[A]] - predict(model_A, data = test_data)$predictions
      }
      residuals_list[[A]] <- residuals
    }
    if (method_name == "gbm"){
      if(is.factor(test_data[[A]]) & length(levels(test_data[[A]])) == 2 ){
        train_data[[A]] <- as.numeric(train_data[[A]]) - 1
        
        default_params <- list(
          formula = formula,
          distribution = "bernoulli",
          data = train_data,
          n.trees = 100,
          shrinkage = 0.05,
          n.minobsinnode = 10
        )
        
        custom_params <- params[["binary"]]
        if (is.null(custom_params)) 
          custom_params <- list()
        final_params <- modifyList(default_params, custom_params)
        
        model_A <- do.call(gbm, final_params)
        residuals <- as.numeric(test_data[[A]]) - 1 - predict(model_A, newdata = test_data, type = "response")
        
        
        residuals_list[[A]] <- residuals
      }
      else if(is.factor(test_data[[A]]) & length(levels(test_data[[A]])) > 2){
        model_A <- gbm(formula, distribution = "multinomial", data = train_data, n.minobsinnode = 10)
        predictions <- predict(model_A, newdata = test_data, type = "response")
        
        one_hot_formula <- as.formula(paste("~", A, "-1"))
        A_one_hot_coding <- model.matrix(one_hot_formula, data = test_data)
        residuals_mat <- A_one_hot_coding - predictions[,,1]
        
        if (grepl("inter", A)){
          if (ncol(A_one_hot_coding) == 4)
            residuals_mat <- residuals_mat[,4]
          if (ncol(A_one_hot_coding) == 6)
            residuals_mat <- residuals_mat[,c(5,6)]
        } else{
          residuals_mat <- residuals_mat[,-1]
        }
        
        new_residuals <- setNames(as.list(as.data.frame(residuals_mat)), colnames(residuals_mat))
        new_residuals
        
        # Append to the existing list without overwriting
        residuals_list <- c(residuals_list, new_residuals)
        # residuals <- as.numeric(test_data[[A]]) - 1 - which.max(predictions, arr.ind = TRUE)[,2]
      }
      else{
        
        default_params <- list(
          formula = formula,
          distribution = "gaussian",
          data = train_data,
          shrinkage = 0.01, 
          n.trees = 500, 
          #interaction.depth = 3, 
          interaction.depth = 5, 
          n.minobsinnode = 1
        )
        
        custom_params <- params[["regression"]]
        if (is.null(custom_params)) 
          custom_params <- list()
        final_params <- modifyList(default_params, custom_params)
        
        model_A <- do.call(gbm, final_params)
        residuals <- test_data[[A]] - predict(model_A, newdata = test_data)
        residuals_list[[A]] <- residuals
      }
    }
    if (method_name == "nnet"){
      if (is.factor(test_data[[A]])) {
        model_A <- nnet(
          formula,
          data = train_data,
          size = 16,               # Number of hidden neurons
          decay = 0.1,           # Regularization parameter
          maxit = 500,           # Maximum number of iterations
          linout = FALSE,         # Logistic output for classification
          trace = FALSE           # Suppress training output
        )
        predictions <- predict(model_A, newdata = test_data, type = "raw")
        residuals <- as.numeric(test_data[[A]]) - 1 - as.numeric(predictions)
      } else {
        # Regression scenario
        final_params <- modifyList(
          list(
            formula,
            data = train_data,
            size = 16,
            decay = 0.1,
            maxit = 500,
            linout = TRUE,
            trace = FALSE
          ),
          params  # override any from tuning grid
        )
        model_A <- do.call(nnet::nnet, final_params)
        
        print(final_params)
        
        predictions <- predict(model_A, newdata = test_data)
        residuals <- test_data[[A]] - predictions
      }
      residuals_list[[A]] <- residuals
    }
    
  }
  #names(residuals_list) <- A_names
  
  return(residuals_list)
}



#function to perform residual regression and variance calculation
DML_residual_reg <- function(input_dat, folds_num, method_name_A, method_name_Y, A_names, params = list(), params_A = list()){
  ###preprocessing for neural network
  
  if (method_name_A == "nnet"){
    #####z-score normalization
    # f_cols <- sapply(dat, is.factor)
    # columns_to_scale <- !f_cols
    # dat_scaled <- input_dat
    # dat_scaled[, columns_to_scale] <- scale(input_dat[, columns_to_scale])
    # dat_scaled_cols <- scale(input_dat[,columns_to_scale])
    # dat = dat_scaled
    
    #####min-max normalization
    f_cols <- sapply(input_dat, is.factor)
    columns_to_scale <- !f_cols
    
    dat_scaled <- input_dat
    
    maxs <- apply(dat_scaled[, columns_to_scale ], 2, max) 
    mins <- apply(dat_scaled[, columns_to_scale ], 2, min)
    
    dat_scaled <- input_dat
    dat_scaled[, columns_to_scale] <- scale(input_dat[, columns_to_scale], center = mins, scale = maxs - mins)
    dat = dat_scaled
    
  } else{
    dat = input_dat
  }
  ################
  #########run estimates
  
  folds <- caret::createFolds(dat$Y, k = folds_num, list = TRUE)
  
  all_residuals <- list()
  residuals_y_all <- c()
  residuals_A_all <- c()
  
  
  # Loop over each fold, treating it as the test set once
  for(i in 1:folds_num) {
    
    # Split the data into training and testing sets
    test_indices <- folds[[i]]
    train_data <- dat[-test_indices, ]
    test_data <- dat[test_indices, ]
    
    residuals_y_temp <- gen_y_residuals(A_names = A_names, train_data, test_data, method_name_Y, params)
    residuals_A_temp <- gen_treatment_residuals(A_names = A_names, train_data, test_data, method_name_A, params_A)
    residuals_A_temp <- do.call(cbind.data.frame, residuals_A_temp)
    #colnames(residuals_A_temp) <- A_names
    
    # Store the residuals for this fold
    residuals_y_all <- c(residuals_y_all, residuals_y_temp)
    residuals_A_all <- rbind(residuals_A_all, residuals_A_temp)
    
    combined_residuals <- data.frame(residuals_y_all, residuals_A_all)
    colnames(combined_residuals)[1] <- "Y"
  }
  
  
  ##################scaling for neural network####################
  if (method_name_Y == "nnet"){
    continuous_A_names <- A_names[!sapply(dat[A_names], is.factor)]
    continuous_vars <- c("Y", continuous_A_names)
    
    ###z-score normalization
    # scaled_scale <- attr(dat_scaled_cols, "scaled:scale")
    # scaled_center <- attr(dat_scaled_cols, "scaled:center")  
    # 
    # combined_residuals$residuals_y_all <- combined_residuals$residuals_y_all * scaled_scale['Y'] + scaled_center['Y']
    # combined_residuals$inc <- combined_residuals$inc * scaled_scale['inc'] + scaled_center['inc']
    # combined_residuals$e401_inc_inter <- combined_residuals$e401_inc_inter * scaled_scale['e401_inc_inter'] + scaled_center['e401_inc_inter']
    
    ###min-max normalization
    for (var in continuous_vars) {
      combined_residuals[[var]] <- combined_residuals[[var]] * 
        (maxs[[var]] - mins[[var]]) + mins[[var]]
    }
  }
  
  
  ########################################
  ####get coefs and variance estimate
  ###########################################
  Ncols <- ncol(combined_residuals)
  Nsubjects <- nrow(combined_residuals)
  
  #################################
  #####get coefs estimate
  lm_fit <- lm(Y ~ ., data = combined_residuals)
  coefs <- lm_fit$coefficients[2:Ncols]
  
  ###################################
  #####get variance estimate
  Y_res <- as.vector(combined_residuals[,1])
  V <- t(as.matrix(combined_residuals[,2:Ncols])) # dimension: K * N
  
  #psi_a <- V%*%t(V)
  #psi_b <- V %*% Y_res
  
  ###all psi_a s
  psi_a_all <- array(0, dim = c(Ncols - 1, Ncols -1, Nsubjects))
  for (i in 1:Nsubjects){
    psi_a_all[,,i] <-  V[,i] %*% t(V[,i])
  }
  
  ###all psi_b s
  psi_b_all <- matrix(0, nrow = Ncols - 1, ncol = Nsubjects)  #dimension: K * N to store all subjects
  for (i in 1:Nsubjects){
    psi_b_all[,i] <- V[,i] * Y_res[i]
  }
  
  #solve(apply(psi_a_all, c(1,2), mean)) %*% rowMeans(psi_b_all) coefficient estimate
  
  ##all overall psi s
  psi_square_all <- array(0, dim = c(Ncols - 1, Ncols -1, Nsubjects))
  for (i in 1:Nsubjects){
    temp <- -psi_a_all[,,i] %*% coefs
    psi_i <- -psi_a_all[,,i] %*% coefs + psi_b_all[,i]
    psi_square_all[,,i] <- psi_i %*% t(psi_i)
  }
  
  ##
  #J * Y_res^2
  ##
  #print(coefs)
  #apply(temp, c(1,2), mean) 
  #print(psi_i)
  
  #mean(Y_res^2)
  
  #str(psi_b_all)
  #rowMeans(psi_i)
  #rowMeans(psi_b_all)
  
  
  J <- apply(psi_a_all, c(1,2), mean)
  #var_mat <- solve(J) %*% apply(psi_square_all, c(1,2), mean) %*% solve(J)
  var_mat <- solve(J) %*% apply(psi_square_all, c(1,2), mean) %*% t(solve(J))
  var_mat <- var_mat / Nsubjects     #normalized variance
  se_mat <- sqrt(var_mat)
  
  return(list(coefs = coefs, var_mat = var_mat, se_mat = se_mat))
}


#make different number of splits and then get a median estimates
get_final_median_estimates <- function(input_dat, folds_num, split_num, method_name_A, method_name_Y, A_names, params = list(), params_A = list()){
  #print(params_A)
  
  cl <- makeCluster(20)
  registerDoSNOW(cl)
  
  
  parallel::clusterExport(cl, varlist = c("DML_residual_reg", "gen_y_residuals", "gen_treatment_residuals"))
  
  # Set up a progress function to update the bar
  pb <- txtProgressBar(max = split_num, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  all_results <- foreach(i = 1:split_num, .combine = rbind, .packages = c("ranger", "glmnet", "gbm", "nnet"), .options.snow = opts) %dopar% {
    temp_res <- DML_residual_reg(input_dat, folds_num = folds_num, method_name_A, method_name_Y, A_names = A_names, params = params, params_A = params_A)
    list(coefs = temp_res$coefs, vars = diag(temp_res$var_mat), var_mat = temp_res$var_mat)
  }
  
  stopCluster(cl)
  close(pb)  # Also close the progress bar if needed
  gc() 
  
  
  all_coefs <- c()
  all_vars <- c()
  for (i in 1:nrow(all_results)) {
    res <- all_results[i, ]
    all_coefs <- rbind(all_coefs, res$coefs)  # Append coefficients
  }
  
  final_coefs <- apply(all_coefs, 2, median)
  
  var_mat_dim <- nrow(all_results[1, ]$var_mat)
  var_mats_array <- array(
    unlist(lapply(1:nrow(all_results), function(i) all_results[i, ]$var_mat)),
    dim = c(var_mat_dim, var_mat_dim, nrow(all_results))
  )
  
  final_varmat <- apply(var_mats_array, c(1, 2), median)
  final_se <- sqrt(diag(final_varmat))
  
  return(list(coefs = final_coefs, se = final_se, varmat = final_varmat))
}



sim_summaries <- function(sim_coefs, Truth){
  mat_subtract_truth <- sweep(as.matrix(sim_coefs[,2:4]), 2, Truth, '-')
  bias <- colMeans(mat_subtract_truth) 
  rMSE <- sqrt(colMeans(mat_subtract_truth**2))
  
  return (data.frame(bias = bias, rMSE = rMSE))
}






########################################
###########Evaluation function
###########################################

extract_coefs <- function(res_list, method_prefix) {
  do.call(rbind, lapply(res_list, function(res) res[, paste0(method_prefix, "_coefs")]))
}

extract_ses <- function(res_list, method_prefix) {
  do.call(rbind, lapply(res_list, function(res) res[, paste0(method_prefix, "_se")]))
}

calc_bias_rmse <- function(estimates, true_vals) {
  bias <- colMeans(estimates) - true_vals
  rmse <- sqrt(colMeans((t(t(estimates) - true_vals))^2))
  data.frame(Bias = bias, rMSE = rmse)
}

calc_confint <- function(estimates, ses, level = 0.95, true_vals = NULL) {
  z <- qnorm(1 - (1 - level) / 2)  # z-score for the CI
  lower <- estimates - z * ses
  upper <- estimates + z * ses
  
  result <- list(lower = lower, upper = upper)
  
  if (!is.null(true_vals)) {
    # Check if true value is between lower and upper bound
    coverage_mat <- (lower <= matrix(true_vals, nrow = nrow(lower), ncol = length(true_vals), byrow = TRUE)) & 
      (upper >= matrix(true_vals, nrow = nrow(lower), ncol = length(true_vals), byrow = TRUE))
    
    # Coverage probability for each parameter (column-wise mean)
    coverage_prob <- colMeans(coverage_mat)
    result$coverage <- coverage_prob
    result$ses <- ses
  }
  
  return(result)
}



extract_coefs_single <- function(res_list, method) {
  map_dfr(res_list, ~{
    tibble(
      variable = rownames(.x[[method]]$coefs),  # Preserve row names
      estimate = .x[[method]]$coefs[, "Estimate"],
      se = .x[[method]]$coefs[, "Std. Error"]
    )
  })
}




########################################
###########Other helper functions
###########################################


get_single_res <- function(A_names, dat, folds_num, split_num, method = "nnet", params = list(), params_A = list()){
  all_res<- c()
  
  # #RF results
  if (method == "rf"){
    method_name_Y = "rf"
    method_name_A = "rf"
    final_results <- get_final_median_estimates(dat, folds_num, split_num, method_name_A, method_name_Y, A_names, params = params, params_A = params_A)
    all_res <- cbind(all_res, final_results$coefs, final_results$se)
  }
  
  #GBM results
  if (method == "gbm"){
    method_name_Y = "gbm"
    method_name_A = "gbm"
    final_results <- get_final_median_estimates(dat, folds_num, split_num, method_name_A, method_name_Y, A_names, params = params, params_A = params_A)
    all_res <- cbind(all_res, final_results$coefs, final_results$se)
  }
  
  # #neural network
  if (method == "nnet"){
    method_name_Y = "nnet"
    method_name_A = "nnet"
    final_results <- get_final_median_estimates(dat, folds_num, split_num, method_name_A, method_name_Y, A_names, params = params, params_A = params_A)
    all_res <- cbind(all_res, final_results$coefs, final_results$se)
  }
  # 
  
  colnames(all_res) <- c("single_coefs", "single_se")
  return (all_res)
}




get_all_method_res <- function(A_names, dat, folds_num, split_num, params = list(), params_A = list()){
  all_res<- c()
  
  #regression model results
  method_name_Y = "regression"
  method_name_A = "regression"
  final_results <- get_final_median_estimates(dat, folds_num, split_num, method_name_A, method_name_Y, A_names)
  all_res <- cbind(all_res, final_results$coefs, final_results$se)
  
  
  #RF results
  method_name_Y = "rf"
  method_name_A = "rf"
  final_results <- get_final_median_estimates(dat, folds_num, split_num, method_name_A, method_name_Y, A_names)
  all_res <- cbind(all_res, final_results$coefs, final_results$se)
  
  #GBM results
  method_name_Y = "gbm"
  method_name_A = "gbm"
  final_results <- get_final_median_estimates(dat, folds_num, split_num, method_name_A, method_name_Y, A_names, params = params, params_A = params_A)
  all_res <- cbind(all_res, final_results$coefs, final_results$se)
  
  
  # #neural network
  method_name_Y = "nnet"
  method_name_A = "nnet"
  final_results <- get_final_median_estimates(dat, folds_num, split_num, method_name_A, method_name_Y, A_names)
  all_res <- cbind(all_res, final_results$coefs, final_results$se)
  
  colnames(all_res) <- c("lm_coefs", "lm_se", "rf_coefs", "rf_se", "gbm_coefs", "gbm_se", "nnet_coefs", "nnet_se")
  return (all_res)
}




#######################################
#######generalized grid search function
######################################



tune_method_general <- function(method_name, tune_grid, data_nums, all_dat,
                                A_names, folds_num, split_num) {
  
  # Save tuning grid for later evaluation
  
  # Process each parameter combination
  for (grid_idx in 1:nrow(tune_grid)) {
    params_row <- tune_grid[grid_idx, ]
    
    if(method_name == "gbm"){
      params_list <- list(
        binary = list(),
        regression = as.list(params_row))
    }else{
      params_list <- as.list(params_row)
    }
    
    # Create generic parameter ID
    param_id <- sprintf("grid_%03d", grid_idx)
    
    # Process all datasets
    for (i in 1:data_nums) {
      current_dat <- all_dat[[i]]
      current_dat$A_1 <- as.factor(current_dat$A_1)
      
      # Call modeling function
      single_res <- get_single_res(
        A_names = A_names,
        dat = current_dat,
        folds_num = folds_num,
        split_num = split_num,
        method = method_name,
        params = list(),
        #params = as.list(params_row),
        #params_A = list() # Pass all parameters as list
        params_A = params_list
      )
      
      # Save with generic ID
      saveRDS(
        single_res,
        file = file.path(
          dir_path, "res", "PLR_sims", "phase1",
          paste0("res_", method_name, "_", param_id, 
                 "_dataset_", sprintf("%04d", i), ".rds")
        )
      )
    }
  }
}
