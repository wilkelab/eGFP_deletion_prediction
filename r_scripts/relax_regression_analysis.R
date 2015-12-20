library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(1) #Comment this out later for real thing

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

#PLOTS FOR TRAINING DATA
plot_prediction_curve <- function(my_data, glm_model, plot_name){  
  predict_data = data.frame(predictor = glm_model$linear.predictors, probability = glm_model$fitted.values, functional_status = my_data$functional_status)
  
  score_tolerated_data = predict_data %>% filter(functional_status == "tolerated")
  score_not_tolerated_data = predict_data %>% filter(functional_status == "non-tolerated")
  
  p = ggplot(data = NULL, aes(x = predictor, y = probability)) + 
    theme_classic() + 
    geom_point(data = predict_data, aes(color = functional_status)) +
    geom_point(data = score_tolerated_data, aes(x = predictor, y = 1, color = functional_status)) + 
    geom_point(data = score_not_tolerated_data, aes(x = predictor, y = 0, color = functional_status)) +
    guides(color=guide_legend(title="Functional Status")) + #Change legend title + 
    xlab("Predictor") +
    ylab("Probability the Mutant is Tolerated") + 
    theme(text = element_text(size = 18)) + 
    theme(axis.line = element_line(size = 1.05)) + 
    theme(axis.ticks = element_line(size = 1.05))
    ggsave(plot_name, plot = p, height = 6, width = 8, useDingbats = F)
    return(predict_data)
}

plot_density <-function(predicted_data, plot_name){
  predictor_density = ggplot(predict_data, aes(x = predictor, fill = functional_status)) + 
    geom_density(alpha = 0.4) + 
    theme_classic() +  
    guides(color=guide_legend(title="Functional Status")) + 
    theme(text = element_text(size = 18)) + 
    theme(axis.line = element_line(size = 1.05)) + 
    theme(axis.ticks = element_line(size = 1.05))
  ggsave(plot_name, plot = predictor_density, height = 6, width = 8, useDingbats = F)
} 


egfp_data = read.csv("../data/egfp_relax_model_scores.csv")
egfp_function_data = read.csv("../data/egfp_functional_data.csv")
egfp_structural_data = read.csv("../data/egfp_structural_data.csv")

sum_data = egfp_data %>% group_by(mutant) %>%summarize(var_score = var(total_score), mean_score = mean(total_score), max_score = max(total_score),  min_score = min(total_score))

sum_data = left_join(sum_data, egfp_structural_data)
sum_data = left_join(sum_data, egfp_function_data)
sum_data = sum_data %>% filter(RSA !="NA", SS != "C-terminus", SS != "Cro", SS != "Loop-Cro")
sum_data$SS = factor(sum_data$SS)

write.csv(sum_data, "../data/relax_full_egfp_summary_data.csv")

ggplot(data = sum_data, aes(x = functional_status, y = mean_score)) + geom_boxplot() + theme_classic()
ggplot(data = sum_data, aes(x = functional_status, y = var_score)) + geom_boxplot() + theme_classic()
ggplot(data = sum_data, aes(x = functional_status, y = max_score)) + geom_boxplot() + theme_classic()
ggplot(data = sum_data, aes(x = functional_status, y = min_score)) + geom_boxplot() + theme_classic()

#Check that the distributions of means and variances are the same/different?
all_tolerated = sum_data %>% filter(functional_status == "tolerated")
all_nontolerated = sum_data %>% filter(functional_status == "non-tolerated")
#t.test(all_tolerated$mean_score, all_nontolerated$mean_score) #Mean of means are different
#t.test(all_tolerated$var_score, all_nontolerated$var_score) #Mean of variances are the same
#t.test(all_tolerated$max_score, all_nontolerated$max_score) #Mean of the maxes are the same
#t.test(all_tolerated$min_score, all_nontolerated$min_score) #Mean of the mins are differebt


ggplot(data = sum_data, aes(x = mean_score, y=var_score, color = functional_status)) + 
  xlab("Mean of Total Rosetta Score") + 
  ylab("Variance of Total Rosetta Score") + 
  theme_classic() + 
  theme(text = element_text(size = 18)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  geom_point()
  

ggplot(data = sum_data, aes(x = mean_score, y=max_score, color = functional_status)) +
  geom_point() + 
  theme_classic() + 
  theme(text = element_text(size = 18)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  xlab("Mean of Total Rosetta Score") + 
  ylab("Max Total Rosetta Score") 


ggplot(data = sum_data, aes(x = mean_score, y=min_score, color = functional_status)) +
  geom_point() + 
  theme_classic() + 
  theme(text = element_text(size = 18)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  xlab("Mean of Total Rosetta Score") + 
  ylab("Min Total Rosetta Score") 

#Turns out that W75 had a rough positive model. Have no idea why.
cor.test(sum_data$mean_score, sum_data$var_score) #r = 0.46 between variance and mean
no_W57_data = sum_data %>% filter(mutant != "W57_del")
cor.test(no_W57_data$mean_score, no_W57_data$var_score) #r = 0.417 between variance and mean

AUC_df <- data.frame(num_replicate = character(), AUC_Type = character(), AUC_Value = double(), stringsAsFactors = FALSE)

j = 1
for (i in 1:100){
  AUC_df[i,1] = as.character(i)
  train_fraction = 0.5 #Percent of data for training the log regression model
  num_datapoints = nrow(sum_data)
  train_size = floor(num_datapoints*train_fraction) #Get the number of data points for training and Round down to an int

  train_indices = sample(1:num_datapoints, size = train_size) #Seperate data into training and test data
  #train_indices
  train_data = sum_data[train_indices,] # training data
  test_data = sum_data[-train_indices,] # test data

  #write.csv(train_data, "../data/new_training_data.csv")
  glm_out_train <-glm(functional_status ~ mean_score, data = train_data, family = binomial) #Logistic regression
  #summary(glm_out_train)
  glm_out_var <-glm(functional_status ~ var_score, data = train_data, family = binomial) #Logistic regression
  #summary(glm_out_var) #Not significant

  glm_max_train <-glm(functional_status ~ max_score, data = train_data, family = binomial) #Logistic regression
  #summary(glm_max_train) #Not significant

  glm_min_train <-glm(functional_status ~ min_score, data = train_data, family = binomial) #Logistic regression
  #summary(glm_min_train) #significant

  glm_min_mean_train <-glm(functional_status ~ min_score + mean_score, data = train_data, family = binomial) #Logistic regression
  #summary(glm_min_mean_train) #significant

  glm_rsa_train <-glm(functional_status ~ RSA, data = train_data, family = binomial) #Logistic regression
  #summary(glm_rsa_train)

  glm_wcn_train <-glm(functional_status ~ wcnSC, data = train_data, family = binomial) #Logistic regression
  #summary(glm_wcn_train)

  glm_score_rsa_train <-glm(functional_status ~ RSA + mean_score, data = train_data, family = binomial) #Logistic regression
  #summary(glm_score_rsa_train)

  glm_ss_train <-glm(functional_status ~ SS, data = train_data, family = binomial) #Logistic regression
  print(summary(glm_ss_train))
  
  
  test_functional_predictions = predict(glm_out_train, test_data, type = "response") #Use the test model to predict functional status of test mutants
  test_rsa_predictions = predict(glm_rsa_train, test_data, type = "response") #Use the test model to predict functional status of test mutants
  test_wcn_predictions = predict(glm_wcn_train, test_data, type = "response") #Use the test model to predict functional status of test mutants
  test_ss_predictions = predict(glm_ss_train, test_data, type = "response") #Use the test model to predict functional status of test mutants
  
  #ROC Values
  ROC_score_train <- calc_ROC(probabilities = glm_out_train$fitted.values, known_truth = train_data$functional_status, model.name = "Score Training Data set")
  AUC_df[j,2] <- "AUC_score_train"
  AUC_df[j,3] <- ROC_score_train %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  j = j + 1
  
  ROC_score_test <- calc_ROC(probabilities = test_functional_predictions, known_truth = test_data$functional_status, model.name = "Score Test Data set")
  AUC_df[j,2] <- "AUC_score_test" 
  AUC_df[j,3] <- ROC_score_test %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  j = j + 1
  
  ROC_rsa_train <- calc_ROC(probabilities = glm_rsa_train$fitted.values, known_truth = train_data$functional_status, model.name = "RSA Training Data set")
  AUC_df[j,2] <- "AUC_rsa_train"
  AUC_df[j,3] <- ROC_rsa_train %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  j = j + 1
  
  ROC_rsa_test <- calc_ROC(probabilities = test_rsa_predictions, known_truth = test_data$functional_status, model.name = "RSA Test Data set")
  AUC_df[j,2] <- "AUC_rsa_test"
  AUC_df[j,3] <- ROC_rsa_test %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  j = j + 1
  
  ROC_wcn_train <- calc_ROC(probabilities = glm_wcn_train$fitted.values, known_truth = train_data$functional_status, model.name = "WCN Training Data set")
  AUC_df[j,2] <- "AUC_wcn_train"
  AUC_df[j,3] <- ROC_wcn_train %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  j = j + 1
  
  ROC_wcn_test <- calc_ROC(probabilities = test_wcn_predictions, known_truth = test_data$functional_status, model.name = "WCN Test Data set")
  AUC_df[j,2] <- "AUC_wcn_test"
  AUC_df[j,3]<- ROC_wcn_test %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  j = j + 1
  
  ROC_ss_train <- calc_ROC(probabilities = glm_ss_train$fitted.values, known_truth = train_data$functional_status, model.name = "WCN Training Data set")
  AUC_df[j,2] <- "AUC_ss_train"
  AUC_df[j,3] <- ROC_ss_train %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  j = j + 1
  
  ROC_ss_test <- calc_ROC(probabilities = test_ss_predictions, known_truth = test_data$functional_status, model.name = "WCN Test Data set")
  AUC_df[j,2] <- "AUC_ss_test"
  AUC_df[j,3]<- ROC_ss_test %>% mutate(delta=false_pos-lag(false_pos)) %>% 
    summarize(AUC=sum(delta*true_pos, na.rm=T))
  j = j + 1
  
  
}

ggplot(data = AUC_df, aes(x = AUC_Type, y = AUC_Value)) + 
  theme(text = element_text(size = 18)) +
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  xlab("AUC Type") +
  ylab("AUC Value") + 
  scale_x_discrete(breaks = c("AUC_rsa_test", "AUC_rsa_train", "AUC_score_test", "AUC_score_train", "AUC_wcn_test", "AUC_wcn_train",  "AUC_ss_test", "AUC_ss_train"), labels = c("RSA Test", "RSA Train", "Mean Score Test", "Mean Score Test", "WCN Test", "WCN Train", "SS Test", "SS Train")) + 
  geom_boxplot() + 
  theme_classic()
