rm(list = ls())

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
library(gtsummary)
library(mice)

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
source("PLR_functions.R")


all_dat <- read.csv(paste0(dir_path, "/data/Final_Nigeria_R3_Data.csv"))


######################################
######### data pre-processing
######################################
model_dat <- all_dat %>% dplyr::select(Age, Gender, Ethnicity, Risk.Alleles, Do_you_smoke_cigarettes, Current_ART_medication, Tenofovir, Current_Medication_besides_ART, Duration_on_ART, 
                                       Recent_CD4._count..cells.mm3., Recent_Viral_Load_Count..copies.mm3., DM, HTN, Others,  
                                       Mean_BMI, Mean_SBP, Mean_DBP, eGFR)  %>% 
  filter(Risk.Alleles %in%  c("0","1","2") ) %>%
  filter(Current_ART_medication != "ABC-3TC-AZT") %>%
  rename(smoking = Do_you_smoke_cigarettes, CD4_count = Recent_CD4._count..cells.mm3., Viral_load = Recent_Viral_Load_Count..copies.mm3., anti_HT = Current_Medication_besides_ART) 

model_dat$anti_HT[model_dat$anti_HT == 2] <- 0
model_dat$Tenofovir[model_dat$Tenofovir == 2] <- 0

model_dat <- model_dat %>%
  mutate(ART_3class = case_when(
    str_detect(str_to_lower(Current_ART_medication), "nvp|efv") ~ "NNRTI",
    str_detect(str_to_lower(Current_ART_medication), "atv/r|lpv/r|drv/r") ~ "Boosted-PI",
    str_detect(str_to_lower(Current_ART_medication), "dtg") ~ "DTG",
    TRUE ~ "Other"
  )) %>%
  dplyr::select(-Current_ART_medication) 


##make those variable into factors
factor_vars <- c("Gender", "Ethnicity", "Risk.Alleles", "smoking", "Viral_load", "Tenofovir", "anti_HT", "DM", "HTN", "ART_3class", "Others")

model_dat <- model_dat %>% 
  mutate(CD4_count = as.numeric(CD4_count), Viral_load = as.numeric(Viral_load)) %>%
  mutate(across(all_of(factor_vars), as.factor))

model_dat$ART_3class <- relevel(model_dat$ART_3class, ref = "NNRTI")


imputed_data <- mice(model_dat[, !names(model_dat) %in% "eGFR"], method = "pmm", m = 5)
imputed_data <- complete(imputed_data)

model_dat <- data.frame(imputed_data, Y = model_dat$eGFR)


####summary statistics
summary_table <- model_dat %>%
  tbl_summary(
    by = ART_3class,  # Group by ART level
    missing = "ifany",  # Show missing values if any
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",   # Median (Q1, Q3) for continuous variables
      all_categorical() ~ "{n} ({p}%)"  # n (%) for categorical variables
    )
  ) %>%
  add_overall() %>%  # Include an overall column
  bold_labels() %>%
  as_gt() %>%
  gt::gtsave(filename = paste0(dir_path, "/summary.docx"))



model_dat$ART_TDF_inter <- interaction(model_dat$ART_3class, model_dat$Tenofovir)


##create dummy variables
one_hot_formula <- as.formula(paste("~ ART_TDF_inter - 1"))
A_one_hot_coding <- model.matrix(one_hot_formula, data = model_dat)


##############################################
####Run DML with gbm as learners, and results
##############################################
method_name_Y = "gbm"
method_name_A = "gbm"
split_num = 100
folds_num = 5
A_names <- c("ART_3class", "Tenofovir", "anti_HT", "ART_TDF_inter")
final_results <- get_final_median_estimates(model_dat, folds_num, split_num, method_name_A, method_name_Y, A_names)


round(final_results$coefs, 1)
round(final_results$se, 1) 

lower_CIs = final_results$coefs - final_results$se  * 1.96
upper_CIs = final_results$coefs + final_results$se  * 1.96


names(final_results$coefs) <- c(
  "ART_bPI", "ART_DTG", "TDF", "HTN", "ARTxTDF_DTG", "ARTxTDF_bPI"
)
colnames(final_results$varmat) <- rownames(final_results$varmat) <- names(final_results$coefs)

beta_hat <- as.numeric(final_results$coefs)
names(beta_hat) <- names(final_results$coefs)
vcov_mat <- as.matrix(final_results$varmat)





# create final effect interpretation table
get_result <- function(art, tdf, htn) {
  contrast <- rep(0, length(beta_hat))
  names(contrast) <- names(beta_hat)
  
  if (art == "bPI") contrast["ART_bPI"] <- 1
  if (art == "DTG") contrast["ART_DTG"] <- 1
  if (tdf == 1) contrast["TDF"] <- 1
  if (htn == 1) contrast["HTN"] <- 1
  if (art == "bPI" && tdf == 1) contrast["ARTxTDF_bPI"] <- 1
  if (art == "DTG" && tdf == 1) contrast["ARTxTDF_DTG"] <- 1
  
  contrast <- as.numeric(contrast)
  est <- sum(contrast * beta_hat)
  se <- sqrt(t(contrast) %*% vcov_mat %*% contrast)
  ci <- est + c(-1.96, 1.96) * as.vector(se)
  
  return(data.frame(
    ART = art, TDF = tdf, HTN = htn,
    Estimate = round(est, 2),
    SE = round(se, 2),
    CI_L = round(ci[1], 2),
    CI_U = round(ci[2], 2)
  ))
}

# Explicitly evaluate all combinations
res1  <- get_result("NNRTI", 0, 0)
res2  <- get_result("NNRTI", 0, 1)
res3  <- get_result("NNRTI", 1, 0)
res4  <- get_result("NNRTI", 1, 1)

res5  <- get_result("bPI", 0, 0)
res6  <- get_result("bPI", 0, 1)
res7  <- get_result("bPI", 1, 0)
res8  <- get_result("bPI", 1, 1)

res9  <- get_result("DTG", 0, 0)
res10 <- get_result("DTG", 0, 1)
res11 <- get_result("DTG", 1, 0)
res12 <- get_result("DTG", 1, 1)


results <- rbind(
  res1, res2, res3, res4,
  res5, res6, res7, res8,
  res9, res10, res11, res12
)

print(results)


# Create a combined variable for Tenofovir and anti_HT
model_dat <- model_dat %>%
  mutate(TDF_HT = paste0("TDF", Tenofovir, "_HT", anti_HT))

# Generate a cross-tabulation of ART_3class vs. TDF_HT
model_dat %>%
  select(ART_3class, TDF_HT) %>%
  tbl_cross(
    row = TDF_HT ,
    col = ART_3class,
    percent = "cell"
  ) %>% 
  as_gt() %>% gt::gtsave(paste0(dir_path, "/trt_summary.docx"))


##True effect in linear regression
lm_res <- lm(Y~. + ART_3class * Tenofovir - ART_TDF_inter, data = model_dat_eGFR)
summary(lm_res)


model_dat_eGFR$ART_TDF_inter

