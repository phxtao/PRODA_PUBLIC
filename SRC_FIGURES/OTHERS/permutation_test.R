## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(interp)
library(gridExtra)
library(grid)
library(ModelMetrics)
library(R.matlab)
library(jcolors)
library(Rfast)
library(scales)

library(cowplot)

library(ggcorrplot)
library(corrplot)

# library(colorplaner)

dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))


#############################################################################
# function to increase vertical spacing between legend keys
#############################################################################
# @clauswilke
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3


#################################################################################
# Data Pathway
#################################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23' 
time_domain = 'whole_time'

data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'


#################################################################################
# Importance Cate to Cate
#################################################################################

permu_cate_mse = read.csv(paste(data_dir_output, 'neural_networking/permutation_test_mse_cate2cate_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_10.csv', sep = ''), header = FALSE)

var_category = c('Soil Physical', 'Soil Chemical', 'Climate', 'Vegetation', 'Geography')
para_category = c('Microbial CUE', 'Baseline Decomposition', 'Environmental Modifier', 'Vertical Transport', 'Input Allocation')

permu_cate_mse_std = permu_cate_mse[ , 6:10]
permu_cate_mse = permu_cate_mse[ , 1:5]
colnames(permu_cate_mse_std) = para_category
colnames(permu_cate_mse) = para_category

current_data = data.frame(array(NA, dim = c(length(var_category)*length(para_category), 4)))
colnames(current_data) = c('Importance', 'Parameter', 'Category', 'Std')

count_num = 0
for (ipara in 1:length(para_category)) {
  for (icate in 1:length(var_category)) {
    count_num = count_num + 1
    current_data[count_num, 'Importance'] = permu_cate_mse[icate, para_category[ipara]]/sum(permu_cate_mse[ , para_category[ipara]])
    current_data[count_num, 'Std'] = permu_cate_mse_std[icate, para_category[ipara]]/sum(permu_cate_mse[ , para_category[ipara]])
    current_data[count_num, 'Parameter'] = para_category[ipara]
    current_data[count_num, 'Category'] = var_category[icate]
  }
}

current_data$Parameter = as.character(current_data$Parameter)

ipara = 1
for (ipara in 1:length(para_category)) {
  
sub_current_data = current_data[current_data$Parameter == para_category[ipara], ]

# p =
  ggplot(data = sub_current_data, aes(x = reorder(Category, Importance), y = Importance)) + 
  geom_bar(width = 0.8, alpha = 1, linetype = 'solid', size = 0.7, color = 'black', fill = 'blue4', stat = 'identity', position = position_dodge(0.8)) +
  # geom_point(shape = 21, size = 7, color = 'white', stat = 'identity', position = position_dodge(0.8)) +
  geom_errorbar(aes(x = reorder(Category, Importance), ymin = Importance - Std, ymax = Importance + Std), width = 0.3, size = 0.7, stat = 'identity', position = position_dodge(0.8)) + 
  coord_flip() +
  theme_classic() +
  scale_x_discrete(name = '') +
  scale_y_continuous(limits=c(0, 0.5), oob = rescale_none) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
  theme(axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) +
  # add title
  labs(title = para_category[ipara], x = '', y = 'Relative Importance') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 20), axis.title = element_text(size = 20), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 


eval(parse(text = paste('p', ipara, ' = p', sep = '')))

}


jpeg(paste('./Ensemble/relative_importance_env_var_to_process.jpeg', sep = ''), width = 18, height = 10, units = 'in', res = 300)
plot_grid(p1, p2, p3, p4, p5,
          nrow = 2, ncol = 3)
dev.off()

