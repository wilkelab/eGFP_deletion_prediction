library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

#set.seed(1) #Comment this out later for real thing

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


model_type = "WCN_Mean_Score"
sum_data = read.csv("../../../data/relax_full_egfp_summary_data.csv")
model_data = sum_data %>% select(RSA, wcnSC, mean_score, SS, functional_status)
num_samples = 100

for (j in 1:num_samples){
  flds <- createFolds(model_data$functional_status, k = 10, list = TRUE, returnTrain = FALSE)
  flat_flds = unlist(unlist(flds, recursive = FALSE))
  num_datasets = length(flds)
  
  Test_Regression_df <- data.frame(Probability = double(), Predictor = double(), Known_Truth = character(), stringsAsFactors = FALSE)
  Train_AUC_df <- data.frame(num_dataset = double(), AUC_Train = double(), stringsAsFactors = FALSE)
  
  for (i in 1:num_datasets){
    test_indices = as.vector(flds[[i]])
    train_data = model_data[-test_indices,] # training data
    test_data = model_data[test_indices,] # test data
    
    glm_train <-glm(functional_status ~ wcnSC + mean_score, data = train_data, family = binomial) #Logistic regression  
    test_prob = predict(glm_train, test_data, type = "response") #Use the trained model to predict functional status of test mutant
    test_predictor = predict(glm_train, test_data) #Use the trained model to predict functional status of test mutant
    train_prob = glm_train$fitted.values
    
    test_roc = data.frame(Probability = test_prob, Predictor = test_predictor, Known_Truth = test_data$functional_status, stringsAsFactors = FALSE)  
    ROC_train <- calc_ROC(probabilities = train_prob, known_truth = train_data$functional_status, model.name = model_type)
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
  
  train_filename = paste("../model_data/training_auc_data/", model_type,"_", toString(j), ".csv", sep = "")
  write.csv(Train_AUC_df, train_filename)
}
