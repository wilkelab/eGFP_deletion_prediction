library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(grid)

sum_data = read.csv("relax_full_egfp_summary_data.csv")
#test_data = sum_data %>% mutate(sheet = SS =="Sheet")
num_points = nrow(sum_data)
data = data.frame(Sheet = character())
for (i in 1:num_points){
  if (sum_data$SS[i] == "Sheet"){
      data = rbind(data, data.frame(Sheet = 1 ))
  }else{
    data = rbind(data, data.frame(Sheet = 0))
  }

}
sum_data = cbind(sum_data, data)

data = data.frame(Helix = character())
for (i in 1:num_points){
  if (sum_data$SS[i] == "Helix"){
    data = rbind(data, data.frame(Helix = 1))
  }else{
    data = rbind(data, data.frame(Helix = 0))
  }
  
}

sum_data = cbind(sum_data, data)

data = data.frame(Loop = character())
for (i in 1:num_points){
  if (sum_data$SS[i] == "Loop"){
    data = rbind(data, data.frame(Loop = 1 ))
  }else{
    data = rbind(data, data.frame(Loop = 0))
  }
  
}  
sum_data = cbind(sum_data, data)


pca_data = sum_data %>% select(-functional_status, -mutant, -SS, -X, -var_score, -max_score, -min_score, -SA, -wcnCA)

pca <- pca_data %>% #remove the mutant info, SS, and function
  scale %>% #scale to 0 mean and unit variance
  prcomp()  #do PCA

pca

all_pca_data = data.frame(pca$x, functional_status = sum_data$functional_status, secondary_structure = sum_data$SS)

p1 = ggplot(all_pca_data, aes(x=PC1, y=PC2)) + 
  geom_point(size = 2, aes(color = functional_status, shape = secondary_structure)) + 
  #theme(text = element_text(size = 18)) + 
  #theme(axis.line = element_line(size = 1.05)) + 
  #theme(axis.ticks = element_line(size = 1.05)) + 
  scale_color_manual(name = c("Functional Status"), breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated"), values = c("dodgerblue4", "firebrick")) + 
  scale_shape_discrete(name = "Secondary Structure") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) + 
  xlim(-4, 4) + ylim(-4, 4)
  #scale_x_continuous(breaks = seq(-5, 5, 2 )) +
  #coord_cartesian(ylim = c(-3, 3)) +
  #scale_y_continuous(breaks = seq(-5, 3, 1 ))


#ggsave("../figures/wcn_rsa_comparison_plot.pdf", plot = p2, height = 5, width = 7, useDingbats = F)


# capture the rotation matrix in a data frame
rotation_data <- data.frame(pca$rotation, variable=row.names(pca$rotation))
# define a pleasing arrow style
arrow_style <- arrow(length = unit(0.05, "inches"),
                     type = "closed")
# now plot, using geom_segment() for arrows and geom_text for labels
#labels =c("Mean Score", "RSA", "WCN", "Helix", "Sheet", "Loop") #Changed these labels...
labels = c("Mean Score", "RSA", "WCN", "Sheet", "Helix", "Loop")
p2 = ggplot(rotation_data) + 
  geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, arrow=arrow_style, size = .5) + 
  geom_text(aes(x=PC1, y=PC2, label=labels), hjust=c(1.05, -.1, 1.1, 1.1, -.1, -0.1), size=4, color='red') + 
  #theme(text = element_text(size = 24)) + 
  #theme(axis.text = element_text(size = 24)) + 
  #theme(axis.title = element_text(size = 24)) + 
  #theme(axis.line = element_line(size = 1.05)) + 
  #theme(axis.ticks = element_line(size = 1.05)) +
  xlim(-1.1,1.1) + 
  ylim(-1.1,1.1) +
  xlab("PC1 Rotation") + ylab("PC2 Rotation")
  #coord_fixed() # fix aspect ratio to 1:1

#labels_b =c("Stabilty Score", "RSA", "WCN", "Helix", "Sheet", "Loop")
labels_b =c("Stabilty Score", "RSA", "WCN", "Sheet", "Helix", "Loop")

p2_b = ggplot(rotation_data) + 
  geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, arrow=arrow_style, size = 1.05) + 
  geom_text(aes(x=PC1, y=PC2, label=labels_b), hjust=-0.15, size=5, color='red') + 
  theme(text = element_text(size = 16)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) +
  xlim(-1.,1.25) + 
  ylim(-1.,1.) +
  coord_fixed() # fix aspect ratio to 1:1

percent = 100*pca$sdev^2/sum(pca$sdev^2) #Find percent variance explained
perc_data <- data.frame(percent=percent, PC=1:length(percent))
perc_data_mod = perc_data[1:5,]

p3 = ggplot(perc_data_mod, aes(x=PC, y=percent)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label=round(percent, 2)), size=4.5, vjust=-.5) + 
  ylim(0, 80) + 
  ylab("Percent") + 
  theme(text = element_text(size = 16)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) 


p4 = ggplot(all_pca_data, aes(x=PC1, y=PC2)) + 
  geom_point(size = 2, aes(color = secondary_structure)) + 
  theme(text = element_text(size = 18)) + 
  theme(axis.line = element_line(size = 1.05)) + 
  theme(axis.ticks = element_line(size = 1.05)) + 
  #scale_color_manual(name = c("Functional Status"), breaks = c("non-tolerated", "tolerated"), labels = c("Non-Tolerated", "Tolerated"), values = c("dodgerblue4", "firebrick")) + 
  scale_x_continuous(breaks = seq(-5, 5, 2 )) +
  coord_cartesian(ylim = c(-3, 3)) +
  scale_y_continuous(breaks = seq(-5, 3, 1 )) 


combo_plot = plot_grid(p1, p2, nrow = 1, labels = c("A", "B"), align = "h", rel_widths = c(1.4, 1))
combo_plot_ss = plot_grid(p4, p2, nrow = 1, labels = c("A", "B", align = "hv"), label_size = 20)
combo_plot_b = plot_grid(p1, p2_b, nrow = 1, labels = c("A", "B", align = "hv"), label_size = 20)

#ggsave("mod_pca_plots.pdf", plot = combo_plot_b, height = 5, width = 14, useDingbats = F)
save_plot("pca_plots.pdf", plot = combo_plot, ncol = 2.4)
#ggsave("pca_plots_ss.pdf", plot = combo_plot_ss, height = 5, width = 14, useDingbats = F)

#ggsave(pca_variance_plot.pdf", plot = p3, height = 5, width = 7, useDingbats = F)
#ggsave(pca_arrow_plot.pdf", plot = p2, height = 8, width = 8, useDingbats = F)
