library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(caret)
library(e1071)


calc_ROC <- function(probabilities, known_truth, model.name=NULL)
{
  outcome <- as.numeric(factor(known_truth))-1
  pos <- sum(outcome) # total known positives
  neg <- sum(1-outcome) # total known negatives
  pos_probs <- outcome*probabilities # probabilities for known positives
  neg_probs <- (1-outcome)*probabilities # probabilities for known negatives
  true_pos <- sapply(probabilities,
                     function(x) sum(pos_probs>=x)/pos) # true pos. rate
  false_pos <- sapply(probabilities,
                      function(x) sum(neg_probs>=x)/neg)
  if (is.null(model.name))
    result <- data.frame(true_pos, false_pos)
  else
    result <- data.frame(true_pos, false_pos, model.name)
  result %>% arrange(false_pos, true_pos)
}


sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, mean_score)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "Mean_Score")
Mean_Score = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, RSA)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "RSA")
RSA = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, RSA, mean_score)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "RSA_Mean_Score")
RSA_Mean_Score = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, RSA, SS)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "RSA_SS")
RSA_SS = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, RSA, SS, mean_score)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "RSA_SS_Mean_Score")
RSA_SS_Mean_Score = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, RSA, wcnSC)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "RSA_WCN")
RSA_WCN = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, RSA, wcnSC, mean_score)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "RSA_WCN_Mean_Score")
RSA_WCN_Mean_Score = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, RSA, wcnSC, SS)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "RSA_WCN_SS")
RSA_WCN_SS = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))


sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, RSA, wcnSC, SS, mean_score)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "RSA_WCN_SS_Mean_Score")
RSA_WCN_SS_Mean_Score = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, SS)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "SS")
SS = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, SS, mean_score)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "SS_Mean_Score")
SS_Mean_Score = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, wcnSC)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "WCN")
WCN = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, wcnSC, mean_score)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "WCN_Mean_Score")
WCN_Mean_Score = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, wcnSC, SS)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "WCN_SS")
WCN_SS = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

sum_data = read.csv("../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(functional_status, wcnSC, SS, mean_score)
model_data_no_function = model_data %>% select(-functional_status)
model <- svm(functional_status~., data = model_data,  type = "C-classification", kernel = "radial", probability = TRUE)
pred_model <- predict(model,model_data_no_function, probability = TRUE)
model_prob = data.frame(attr(pred_model, "probabilities"))
ROC_model <- calc_ROC(probabilities = model_prob$tolerated, known_truth = model_data$functional_status, model.name = "WCN_SS_Mean_Score")
WCN_SS_Mean_Score = ROC_model %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

combined_AUC = rbind(Mean_Score,RSA,RSA_Mean_Score,RSA_SS,RSA_SS_Mean_Score,RSA_WCN,RSA_WCN_Mean_Score,RSA_WCN_SS,RSA_WCN_SS_Mean_Score,SS,SS_Mean_Score,WCN,WCN_Mean_Score,WCN_SS,WCN_SS_Mean_Score)

#names = rbind("Mean_Score","RSA","RSA_Mean_Score","RSA_SS","RSA_SS_Mean_Score","RSA_WCN","RSA_WCN_Mean_Score","RSA_WCN_SS","RSA_WCN_SS_Mean_Score","SS","SS_Mean_Score","WCN","WCN_Mean_Score", "WCN_SS","WCN_SS_Mean_Score")

#all_AUC_data = data.frame(model_names, combined_AUC)
#names(all_AUC_data) <- c("model_names", "Model_AUC")
all_AUC_data <- combined_AUC %>% arrange(desc(AUC))
#names(all_AUC_data) <- c("model_names", "Model_AUC")
write.csv(all_AUC_data, "../../data/svm_auc_of_model_data.csv")
