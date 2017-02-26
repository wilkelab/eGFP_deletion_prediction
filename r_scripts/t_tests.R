library(dplyr)
library(ggplot2)
library(cowplot)

sum_data = read.csv("../data/relax_full_egfp_summary_data.csv")
tolerated_data = sum_data %>% filter(functional_status == "tolerated")
non_tolerated_data = sum_data %>% filter(functional_status == "non-tolerated")

sheet_data = sum_data %>% filter(SS == "Sheet")
helix_data = sum_data %>% filter(SS == "Helix")
loop_data = sum_data %>% filter(SS == "Loop")


t.test(tolerated_data$wcnSC, non_tolerated_data$wcnSC)
t.test(tolerated_data$RSA, non_tolerated_data$RSA)
t.test(tolerated_data$mean_score, non_tolerated_data$mean_score)


num_mutants = nrow(sum_data)

sheet_counts = sheet_data %>% group_by(functional_status) %>% tally() 
sheet_counts = sheet_counts %>% mutate(percent = n/sum(sheet_counts$n))
helix_counts = helix_data %>% group_by(functional_status) %>% tally() 
helix_counts = helix_counts %>% mutate(percent = n/sum(helix_counts$n))
loop_counts = loop_data %>% group_by(functional_status) %>% tally()
loop_counts = loop_counts %>% mutate(percent = n/sum(loop_counts$n))

tolerated_helix_data = sum_data %>% filter(SS == "Helix", functional_status == "tolerated")
non_tolerated_helix_data = sum_data %>% filter(SS == "Helix", functional_status == "non-tolerated")

tolerated_sheet_data = sum_data %>% filter(SS == "Sheet", functional_status == "tolerated")
non_tolerated_sheet_data = sum_data %>% filter(SS == "Sheet", functional_status == "non-tolerated")

tolerated_loop_data = sum_data %>% filter(SS == "Loop", functional_status == "tolerated")
non_tolerated_loop_data = sum_data %>% filter(SS == "Loop", functional_status == "non-tolerated")


t.test(sheet_data$wcnSC, loop_data$wcnSC)
t.test(helix_data$wcnSC, loop_data$wcnSC)
t.test(helix_data$wcnSC, sheet_data$wcnSC)

t.test(sheet_data$RSA, loop_data$RSA)
t.test(helix_data$RSA, loop_data$RSA)
t.test(helix_data$RSA, sheet_data$RSA)
