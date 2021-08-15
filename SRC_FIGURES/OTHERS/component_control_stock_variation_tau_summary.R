library(R.matlab)
library(jcolors)
library(maps)
library(ggplot2)
library(Metrics)
library(gridExtra)
library(grid)

library(cowplot)
dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))
discrete.colors <- c('#3A3A3A', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')

#############################################################################
# Data Path
#############################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'
cross_valid_num = 10

data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'

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
# Load Projected SOC (different depth)
#################################################################################
stat_matric = c('rmse', 'mae', 'coeff_efficiency', 'adjust_coeff_efficiency', 'index_agree', 'adjust_index_agree', 'rsr', 'pbias', 'concordance')
stat_matric_index = 3

component_name = c('G_NN', 'A_A', 'B_I', 'C_K', 'D_V', 'E_Xi', 'Fa_Homo', 'Fb_Default')

depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')

control_summary = data.frame(array(NA, dim = c(length(depth_name)*length(component_name), 6)))
colnames(control_summary) = c('var_soc', 'std_soc', 'var_stat', 'std_stat', 'depth', 'component')

icross_valid = 1
counter = 1

for (iquantile in 1:4) {
  for (icomponent in 1:length(component_name)) {
    control_matric_total_stock = array(NA, dim = c(cross_valid_num, 1))
    control_matric_stat = array(NA, dim = c(cross_valid_num, 1))
    
    for (icross_valid in 1:cross_valid_num) {
      # soc global stock
      control_total_stock = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_soc_total_stock_', model_name, '_', nn_exp_name, '_cross_valid_', as.character(icross_valid), '_control_test.mat', sep = ''))
      control_total_stock = control_total_stock$global.stock
      # stock 0 - 200 to 100 - 200cm
      control_total_stock[3, ] = control_total_stock[3, ] - control_total_stock[2, ] 
      # stock 0 - 100 to 30 - 100cm
      control_total_stock[2, ] = control_total_stock[2, ] - control_total_stock[1, ] 
      # record
      control_matric_total_stock[icross_valid] = control_total_stock[iquantile, icomponent]
      
      
      control_test = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_spatial_variation_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_', as.character(icross_valid), '.mat', sep = ''))
      control_test = control_test$compotent.matric # [iquantile, icomponent, imatric]
      control_matric_stat[icross_valid] = control_test[iquantile, icomponent, stat_matric_index]*100
      
    }
    
    control_summary[counter, 'var_soc'] = mean(control_matric_total_stock, na.rm = TRUE)
    control_summary[counter, 'std_soc'] = sd(control_matric_total_stock, na.rm = TRUE)
    
    control_summary[counter, 'var_stat'] = mean(control_matric_stat, na.rm = TRUE)
    control_summary[counter, 'std_stat'] = sd(control_matric_stat, na.rm = TRUE)
    
    control_summary[counter, 'depth'] = depth_name[iquantile]
    control_summary[counter, 'component'] = component_name[icomponent]
    
    if (icomponent == 1 & iquantile == 4) {
      estimates_proda  = control_matric_total_stock
    }
    
    if (icomponent == 2 & iquantile == 4) {
      estimates_A  = control_matric_total_stock
      print(mean(estimates_A - estimates_proda))
      print(sd(estimates_A - estimates_proda))
    }
    
    counter = counter + 1
    
  }
}

### plot figure component only
current_data_range = control_summary[which(control_summary$depth == 'D_full_depth' & control_summary$component == 'G_NN'), ]
current_data = control_summary[which(control_summary$depth == 'D_full_depth' & control_summary$component != 'Fa_Homo' & control_summary$component != 'Fb_Default' & control_summary$component != 'G_NN'), ]

abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:5]
mean(abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:5])

abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:5]
mean(abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:5])
# color_scheme = c('#1d3554', '#DFE07C', '#7F8E39', '#42858C', '#E48F1B', '#D33B44') #jcolors(palette = c("pal7"))
color_scheme = c('#E69F00', '#D55E00', '#F0E442', '#0072B2', '#56B4E9', NA)

jpeg('./Ensemble/component_control_stock_variation_global_summary_CESM2.jpeg', width = 10, height = 10, units = 'in', res = 300)
ggplot(data = current_data) + 
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = current_data_range$var_stat - current_data_range$std_stat, ymax =  current_data_range$var_stat + current_data_range$std_stat, fill = 'grey45', alpha = 0.5) +
  annotate('rect', ymin = -Inf, ymax = Inf, xmin = current_data_range$var_soc - current_data_range$std_soc, xmax =  current_data_range$var_soc + current_data_range$std_soc, fill = 'grey45', alpha = 0.5) +
  geom_hline(yintercept = current_data_range$var_stat, size = 2, color = 'grey45', alpha = 1) +
  geom_vline(xintercept = current_data_range$var_soc, size = 2, color = 'grey45', alpha = 1) +
  geom_point(aes(x = var_soc, y = var_stat, color = component, fill = component), shape = 15, size = 10) + 
  geom_errorbar(aes(x = var_soc, ymin = var_stat - std_stat, ymax = var_stat + std_stat, color = component), size = 2, width = 40, stat = 'identity') +
  geom_errorbarh(aes(y = var_stat, xmin = var_soc - std_soc, xmax = var_soc + std_soc, color = component), size = 2, height = 0.4, stat = 'identity') + 
  scale_color_manual(name = '', labels = c('Microbial CUE', 'Input Allocation', 'Baseline Decomposition', 'Vertical Transportation', 'Environmental Modifers', 'PRODA'), values = color_scheme) +
  # scale_x_continuous(limits = c(1750, 2800)) +
  theme_classic() +
  # change the legend properties
  theme(legend.position = 'none') +
  # add title
  labs(x = 'Global SOC Stock (Pg C)', y = expression(paste('Explained Variation (%)'), sep = '')) +
  # modify the position of title
  # modify the font size
  theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=35))

dev.off()
