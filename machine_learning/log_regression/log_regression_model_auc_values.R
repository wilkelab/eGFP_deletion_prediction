library(plyr)
library(dplyr)
library(ggplot2)
#library(tidyr)
library(cowplot)

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
glm_rsa_mean_score <-glm(functional_status ~ RSA + mean_score, data = sum_data, family = binomial) #Logistic regression  
glm_rsa_ss_mean_score <-glm(functional_status ~ RSA + SS + mean_score, data = sum_data, family = binomial) #Logistic regression  
glm_rsa_ss <-glm(functional_status ~ RSA + SS, data = sum_data, family = binomial) #Logistic regression  
glm_rsa_wcn_mean_score <-glm(functional_status ~ RSA + wcnSC + mean_score, data = sum_data, family = binomial) #Logistic regression  
glm_rsa_wcn_ss_mean_score <-glm(functional_status ~ RSA + wcnSC + SS + mean_score, data = sum_data, family = binomial) #Logistic regression  
glm_rsa_wcn_ss <-glm(functional_status ~ RSA + wcnSC + SS, data = sum_data, family = binomial) #Logistic regression  
glm_rsa_wcn <-glm(functional_status ~ RSA + wcnSC, data = sum_data, family = binomial) #Logistic regression  
glm_rsa <-glm(functional_status ~ RSA, data = sum_data, family = binomial) #Logistic regression  

glm_ss_mean_score <-glm(functional_status ~ SS + mean_score, data = sum_data, family = binomial) #Logistic regression  
glm_ss <-glm(functional_status ~ SS, data = sum_data, family = binomial) #Logistic regression  

glm_wcn_mean_score <-glm(functional_status ~ wcnSC + mean_score, data = sum_data, family = binomial) #Logistic regression  
glm_ss_mean_score <-glm(functional_status ~ SS + mean_score, data = sum_data, family = binomial) #Logistic regression  

glm_wcn_mean_score <-glm(functional_status ~ wcnSC + mean_score, data = sum_data, family = binomial) #Logistic regression  
glm_wcn_ss_mean_score <-glm(functional_status ~ wcnSC + SS + mean_score, data = sum_data, family = binomial) #Logistic regression  
glm_wcn_ss <-glm(functional_status ~ wcnSC + SS, data = sum_data, family = binomial) #Logistic regression  
glm_wcn <-glm(functional_status ~ wcnSC, data = sum_data, family = binomial) #Logistic regression  

glm_mean_score <-glm(functional_status ~ wcnSC + mean_score, data = sum_data, family = binomial) #Logistic regression  

ROC_rsa_mean_score <- calc_ROC(probabilities = glm_rsa_mean_score$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_rsa_ss_mean_score <- calc_ROC(probabilities = glm_rsa_ss_mean_score$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_rsa_ss <- calc_ROC(probabilities = glm_rsa_ss$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_rsa_wcn_mean_score <- calc_ROC(probabilities = glm_rsa_wcn_mean_score$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_rsa_wcn_ss_mean_score <- calc_ROC(probabilities = glm_rsa_wcn_ss_mean_score$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_rsa_wcn_ss <- calc_ROC(probabilities = glm_rsa_wcn_ss$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_rsa_wcn <- calc_ROC(probabilities = glm_rsa_wcn$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_rsa <- calc_ROC(probabilities = glm_rsa$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_ss_mean_score <- calc_ROC(probabilities = glm_ss_mean_score$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_ss <- calc_ROC(probabilities = glm_ss$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_wcn_mean_score <- calc_ROC(probabilities = glm_wcn_mean_score$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_wcn_ss_mean_score <- calc_ROC(probabilities = glm_ss_mean_score$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_wcn_ss <- calc_ROC(probabilities = glm_wcn_ss$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_wcn <- calc_ROC(probabilities = glm_wcn$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")
ROC_mean_score <- calc_ROC(probabilities = glm_mean_score$fitted.values, known_truth = sum_data$functional_status, model.name = "RSA Training Data set")

AUC_rsa_mean_score = ROC_rsa_mean_score %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

AUC_rsa_ss_mean_score = ROC_rsa_ss_mean_score %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

AUC_rsa_ss = ROC_rsa_ss %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_rsa_wcn_mean_score = ROC_rsa_wcn_mean_score %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_rsa_wcn_ss_mean_score = ROC_rsa_wcn_ss_mean_score %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_rsa_wcn_ss = ROC_rsa_wcn_ss %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_rsa_wcn = ROC_rsa_wcn %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_rsa = ROC_rsa %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_ss_mean_score = ROC_ss_mean_score %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_ss = ROC_ss %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_wcn_mean_score = ROC_wcn_mean_score %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_wcn_ss_mean_score = ROC_wcn_ss_mean_score %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_wcn_ss = ROC_wcn_ss %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_wcn = ROC_wcn %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))
AUC_mean_score = ROC_mean_score %>% mutate(delta=false_pos-lag(false_pos)) %>% 
  summarize(AUC=sum(delta*true_pos, na.rm=T))

combined_AUC = rbind(AUC_rsa_mean_score,AUC_rsa_ss_mean_score,AUC_rsa_ss,AUC_rsa_wcn_mean_score,AUC_rsa_wcn_ss_mean_score,AUC_rsa_wcn_ss,AUC_rsa_wcn,AUC_rsa,AUC_ss_mean_score,AUC_ss,AUC_wcn_mean_score,AUC_wcn_ss_mean_score,AUC_wcn_ss,AUC_wcn,AUC_mean_score )

names = rbind("AUC_rsa_mean_score","AUC_rsa_ss_mean_score","AUC_rsa_ss","AUC_rsa_wcn_mean_score","AUC_rsa_wcn_ss_mean_score","AUC_rsa_wcn_ss","AUC_rsa_wcn","AUC_rsa","AUC_ss_mean_score","AUC_ss","AUC_wcn_mean_score","AUC_wcn_ss_mean_score","AUC_wcn_ss","AUC_wcn","AUC_mean_score")

all_AUC_data = data.frame(names, combined_AUC)
names(all_AUC_data) <- c("model_names", "Model_AUC")
all_AUC_data <- all_AUC_data %>% arrange(desc(Model_AUC))
#all_AUC_data = all_AUC_data %>% rename(c(names = "model_names", AUC = "AUC_Model"))
write.csv(all_AUC_data, "../../data/log_regression_auc_of_model_data.csv")
