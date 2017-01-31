library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)

svm_data = read.csv("../data/svm_summarized_auc_cross_validated.csv")
svm_data = svm_data %>% select(-X, -se_AUC_Cross_Validated)
colnames(svm_data)[2] = "SVM_Mean_AUC_Cross_Validated"
reg_data = read.csv("../data/log_regression_summarized_auc_cross_validated.csv")
reg_data = reg_data %>% select(-X, -se_AUC_Cross_Validated)
colnames(reg_data)[2] = "LR_Mean_AUC_Cross_Validated"

all_data = merge(reg_data, svm_data)

p = ggplot(all_data, aes(x = LR_Mean_AUC_Cross_Validated, y = SVM_Mean_AUC_Cross_Validated )) + 
  #theme_classic() +
  xlab("Log Regression Mean Cross-validated AUC") + 
  ylab("SVM Mean Cross-validated AUC") + 
  geom_point(size = 2 ) + 
  theme(text = element_text(size = 16)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  coord_cartesian(xlim = c(6, 1.0)) +
  coord_cartesian(ylim = c(0.55, 1.0)) +
  scale_x_continuous(breaks = seq(0.6, 1.0, 0.1 )) +
  scale_y_continuous(breaks = seq(0.6, 1.0, 0.1 )) +
  expand_limits(x = 0.55) + 
  expand_limits(x = 1.0) + 
  geom_abline(intercept=0,slope=1,lwd=1.05,lty=2,col="gray")# +
  #geom_smooth(data = all_data, method = "lm", se = FALSE, color = "black", size = 1.05) 

ggsave("../figures/svm_reg_comparison.pdf", plot = p, height = 5, width = 6, useDingbats = F)