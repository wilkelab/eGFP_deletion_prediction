library(plyr)
library(dplyr)
library(caret)
library(ggplot2)
library(e1071)

#set.seed(111)

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

model_type = "WCN_SS_Mean_Score"
sum_data = read.csv("../../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(wcnSC, SS, mean_score, functional_status)
num_samples = 100

for (j in 1:num_samples){
  flds <- createFolds(model_data$functional_status, k = 10, list = TRUE, returnTrain = FALSE)
  flat_flds = unlist(unlist(flds, recursive = FALSE))
  num_datasets = length(flds)
  
  Test_Regression_df <- data.frame(Probability = double(), Known_Truth = character(), stringsAsFactors = FALSE)
  Train_AUC_df <- data.frame(num_dataset = double(), AUC_Train = double(), stringsAsFactors = FALSE)
  
  #Do num_replicate and Model_Type later
  for (i in 1:num_datasets){
    test_indices = as.vector(flds[[i]])
    train_data = model_data[-test_indices,] # training data
    test_data = model_data[test_indices,] # test data
    test_data_no_function = test_data %>% select(-functional_status)
    train_data_no_function = train_data %>% select(-functional_status)
    
    #Want probabilities that the deletions are tolerated
    model <- svm(functional_status~., data = train_data,  type = "C-classification", kernel = "radial", probability = TRUE)
    pred_train <- predict(model,train_data_no_function, probability = TRUE)
    pred_test <- predict(model,test_data_no_function, probability = TRUE)
    train_prob = data.frame(attr(pred_train, "probabilities"))
    test_prob = data.frame(attr(pred_test, "probabilities"))
    
    test_roc = data.frame(Probability = test_prob$tolerated, Known_Truth = test_data$functional_status, stringsAsFactors = FALSE)  
    ROC_train <- calc_ROC(probabilities = train_prob$tolerated, known_truth = train_data$functional_status, model.name = "RSA_WCN_SS_Mean_Score")
    AUC_train = ROC_train  %>% group_by(model.name)  %>% mutate(delta=false_pos-lag(false_pos)) %>% 
      summarize(AUC=sum(delta*true_pos, na.rm=T))
    train_auc <- data.frame(num_dataset = i, AUC_Train = AUC_train, stringsAsFactors = FALSE)
    
    Test_Regression_df = rbind(Test_Regression_df, test_roc) 
    Train_AUC_df = rbind(Train_AUC_df, train_auc)
  }
  
  Test_Regression_df = Test_Regression_df %>% mutate(Model_Type = model_type)
  Train_AUC_df = Train_AUC_df %>% mutate(Model_Type = model_type)
  
  ROC_test <- calc_ROC(probabilities = Test_Regression_df$Probability, known_truth = Test_Regression_df$Known_Truth, model.name = model_type)
  AUC_test = ROC_test  %>% group_by(model.name)  %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  #print(AUC_test)
  
  test_filename = paste("../model_data/test_roc_data/", model_type, "_", toString(j), ".csv", sep = "")
  write.csv(Test_Regression_df, test_filename)
  
  train_filename = paste("../model_data/training_auc_data/", model_type, "_", toString(j), ".csv", sep = "")
  write.csv(Train_AUC_df, train_filename)
}