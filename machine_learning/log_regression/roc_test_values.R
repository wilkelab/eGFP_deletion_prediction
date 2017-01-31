library(dplyr)

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



Test_AUC_Values = data.frame(Model_Type = character(), AUC_Cross_Validated = double(), stringsAsFactors = FALSE)

#models = c("RSA", "WCN")
models = c("RSA_WCN_Mean_Score","RSA_SS_Mean_Score", "RSA_SS",  "RSA_WCN_SS", "Mean_Score", "RSA_Mean_Score", "RSA_WCN", "RSA", "RSA_WCN_SS_Mean_Score","SS_Mean_Score", "SS",  "WCN_Mean_Score", "WCN", "WCN_SS", "WCN_SS_Mean_Score")
num_files = length(models)

num_samples = 100
for (i in 1:num_files){
  model_type = models[i]
  for (j in 1:num_samples){
    test_filename = paste("model_data/test_roc_data/", model_type, "_", toString(j),".csv", sep = "")
    Test_Regression_df = read.csv(test_filename)
    ROC_test <- calc_ROC(probabilities = Test_Regression_df$Probability, known_truth = Test_Regression_df$Known_Truth, model.name = model_type)
    AUC_test = ROC_test %>% group_by(model.name) %>% mutate(delta=false_pos-lag(false_pos)) %>% 
      summarize(AUC=sum(delta*true_pos, na.rm=T))
    #print (AUC_test)
    Test_AUC_Values = rbind(Test_AUC_Values, data.frame(Model_Type = AUC_test$model.name, AUC_Cross_Validated = AUC_test$AUC))
    #roc_filename = paste("../../log_regression/roc_curve_data/", model_type, ".csv", sep = "")
    #write.csv(ROC_test, roc_filename)
  }
}

#print(Test_AUC_Values)
write.csv(Test_AUC_Values, "../../data/log_regression_auc_cross_validated.csv")
summarized_auc = Test_AUC_Values %>% group_by(Model_Type) %>% summarize(Mean_AUC_Cross_Validated = mean(AUC_Cross_Validated), se_AUC_Cross_Validated = sd(AUC_Cross_Validated)/sqrt(num_samples)) %>% arrange(desc(Mean_AUC_Cross_Validated))
write.csv(summarized_auc, "../../data/log_regression_summarized_auc_cross_validated.csv")
