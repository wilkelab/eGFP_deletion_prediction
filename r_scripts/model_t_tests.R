library(data.table)
library(dplyr)


log_model_data <- fread("../data/log_regression_auc_cross_validated.csv",sep = ",")
log_model_data <- as.data.frame(log_model_data) %>% select(-V1)

svm_rsa_wcn_mean_score_model <- log_model_data %>% filter(Model_Type == "RSA_WCN_Mean_Score")
svm_rsa_wcn_ss_mean_score_model <- log_model_data %>% filter(Model_Type == "RSA_WCN_SS_Mean_Score")

t.test(log_rsa_wcn_mean_score_model$AUC_Cross_Validated, log_rsa_wcn_ss_mean_score_model$AUC_Cross_Validated)

svm_model_data <- fread("../data/svm_auc_cross_validated.csv",sep = ",")
svm_model_data <- as.data.frame(svm_model_data) %>% select(-V1)

svm_rsa_wcn_mean_score_model <- svm_model_data %>% filter(Model_Type == "RSA_WCN_Mean_Score")
svm_rsa_wcn_ss_mean_score_model <- svm_model_data %>% filter(Model_Type == "RSA_WCN_SS_Mean_Score")

t.test(svm_rsa_wcn_mean_score_model$AUC_Cross_Validated, svm_rsa_wcn_ss_mean_score_model$AUC_Cross_Validated)
