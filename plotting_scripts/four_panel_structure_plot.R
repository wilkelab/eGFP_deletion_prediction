library(dplyr)
library(ggplot2)
library(cowplot)

sum_data = read.csv("../data/relax_full_egfp_summary_data.csv")
tolerated_data = sum_data %>% filter(functional_status == "tolerated")
non_tolerated_data = sum_data %>% filter(functional_status == "non-tolerated")

structure_tally = sum_data %>% group_by(SS) %>% tally()
ss_counts = sum_data %>% group_by(SS, functional_status) %>% tally()

helix_count = ss_counts %>% filter(SS == "Helix") 
h_count = as.integer(structure_tally[1,2])
helix_data = helix_count %>% mutate(SS_freqs = n/h_count)

beta_count = ss_counts %>% filter(SS == "Sheet") 
b_count = as.integer(structure_tally[3,2])
beta_data = beta_count %>% mutate(SS_freqs = n/b_count)

loop_count = ss_counts %>% filter(SS == "Loop") 
l_count = as.integer(structure_tally[2,2])
loop_data = loop_count %>% mutate(SS_freqs = n/l_count)

all_data = rbind(helix_data, beta_data)
all_data = rbind(all_data, loop_data)
all_data$SS <- factor(all_data$SS, levels=c("Loop", "Helix", "Sheet"))

p1 = ggplot(sum_data, aes(x = functional_status, y = RSA, fill = functional_status)) +
  geom_boxplot() + 
  xlab("Functional Status") + 
  ylab("RSA of Deleted Residue") + 
  theme(text = element_text(size = 16)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  scale_x_discrete(breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated") ) +
  scale_fill_manual(name = c("Functional Status"), breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated"), values = c("dodgerblue4", "firebrick")) + 
  coord_cartesian(ylim = c(-0.1, 1.4)) + 
  scale_y_continuous(breaks = seq(0.0, 1.4, 0.2 )) 

p2 = ggplot(sum_data, aes(x = functional_status, y = wcnSC, fill = functional_status)) +
  geom_boxplot() + 
  xlab("Functional Status") + 
  ylab("WCN of Deleted Residue") + 
  theme(text = element_text(size = 16)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  scale_x_discrete(breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated") ) + 
  scale_fill_manual(name = c("Functional Status"), breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated"), values = c("dodgerblue4", "firebrick")) + 
  coord_cartesian(ylim = c(-0.1, 1.4)) + 
  scale_y_continuous(breaks = seq(0.0, 1.4, 0.2 )) 

p3 = ggplot(all_data, aes(x = SS, y = SS_freqs, fill = functional_status)) + 
  xlab("Secondary Structure") +
  ylab("Fraction of Residues in Structural Class") + 
  theme(text = element_text(size = 16)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(name = c("Functional Status"), breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated"), values = c("dodgerblue4", "firebrick"))

p4 = ggplot(sum_data, aes(x = functional_status, y = mean_score, fill = functional_status)) +
  geom_boxplot() + 
  xlab("Functional Status") + 
  ylab("Mean Score") + 
  theme(text = element_text(size = 16)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  scale_x_discrete(breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated") ) + 
  scale_fill_manual(name = c("Functional Status"), breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated"), values = c("dodgerblue4", "firebrick")) + 
  coord_cartesian(ylim = c(-365, -335)) + 
  scale_y_continuous(breaks = seq(-365, -335, 5 ))

p5 = ggplot(sum_data, aes(x = functional_status, y = mean_score, fill = functional_status)) +
  geom_boxplot() + 
  xlab("Functional Status") + 
  ylab("Mean Score") + 
  theme(text = element_text(size = 16, face = "bold")) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  scale_x_discrete(breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated") ) + 
  scale_fill_manual(name = c("Functional Status"), breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated"), values = c("dodgerblue4", "firebrick")) + 
  coord_cartesian(ylim = c(-365, -335)) + 
  scale_y_continuous(breaks = seq(-365, -335, 5 ))

combo_plot = plot_grid(p1, p2, nrow = 1, labels = c("A", "B", align = "hv"), label_size = 20)
ggsave("../figures/structure_boxplots.pdf", plot = combo_plot, height = 5, width = 12, useDingbats = F)
ggsave("../figures/mean_score_boxplot.pdf", plot = p5, height = 5, width = 6, useDingbats = F)
ggsave("../figures/ss_barplot.pdf", plot = p3, height = 5, width = 7, useDingbats = F)

four_plot = plot_grid(p1, p2, p3, p4, nrow = 2, labels = c("A", "B", "C", "D", align = "hv"), label_size = 20)
ggsave("../figures/four_panel_structure.pdf", plot = four_plot, height = 12, width = 14, useDingbats = F)

