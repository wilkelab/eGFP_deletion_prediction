library(dplyr)
library(ggplot2)
library(tidyr)


egfp_data = read.csv("../data/egfp_relax_model_scores.csv")
egfp_function_data = read.csv("../data/egfp_functional_data.csv")
egfp_structural_data = read.csv("../data/egfp_structural_data.csv")

sum_data = egfp_data %>% group_by(mutant) %>%summarize(var_score = var(total_score), mean_score = mean(total_score), max_score = max(total_score),  min_score = min(total_score))

sum_data = left_join(sum_data, egfp_structural_data)
sum_data = left_join(sum_data, egfp_function_data)
sum_data = sum_data %>% filter(RSA !="NA", SS != "C-terminus", SS != "Cro", SS != "Loop-Cro")
sum_data$SS = factor(sum_data$SS)

write.csv(sum_data, "../data/relax_full_egfp_summary_data.csv")
