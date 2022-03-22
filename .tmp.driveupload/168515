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
bootstrap_num = 200

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
# importance of mechanisms to soc storage and spatial variation
#################################################################################
stat_matric = c('rmse', 'mae', 'coeff_efficiency', 'adjust_coeff_efficiency', 'index_agree', 'adjust_index_agree', 'rsr', 'pbias', 'concordance')
stat_matric_index = 3

component_name = c('H_NN', 'A_A', 'F_I', 'D_K', 'E_V', 'B_Xi', 'C_NPP', 'Ga_Homo', 'Gb_Default')

depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')

control_summary = data.frame(array(NA, dim = c(length(depth_name)*length(component_name), 8)))
colnames(control_summary) = c('var_soc', 'upper_soc', 'lower_soc', 'var_stat', 'upper_stat', 'lower_stat', 'depth', 'component')

iquantile = 4
ibootstrap = 1
icomponent = 2
counter = 1
for (iquantile in 1:4) {
  for (icomponent in 1:length(component_name)) {
    control_matric_total_stock = array(NA, dim = c(bootstrap_num, 1))
    control_matric_stat = array(NA, dim = c(bootstrap_num, 1))
    
    for (ibootstrap in 1:bootstrap_num) {
      print(paste('processing quantile ', iquantile, ' component ', icomponent, ' bootstrap  ', ibootstrap,  sep = ''))
      # soc global stock
      control_total_stock = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_soc_total_stock_', model_name, '_', nn_exp_name, '_bootstrap_', as.character(ibootstrap), '_control_test.mat', sep = ''))
      control_total_stock = control_total_stock$var.data.middle
      # record
      control_matric_total_stock[ibootstrap] = control_total_stock[iquantile, icomponent]

      
      control_test = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_spatial_variation_', model_name, '_', time_domain, '_', nn_exp_name, '_bootstrap_', as.character(ibootstrap), '.mat', sep = ''))
      control_test = control_test$var.data.middle # [iquantile, icomponent, imatric]
      
      control_matric_stat[ibootstrap] = control_test[iquantile, icomponent, stat_matric_index]*100
      
    }
    # summary
    control_summary[counter, 'var_soc'] = median(control_matric_total_stock, na.rm = TRUE)
    control_summary[counter, 'upper_soc'] = quantile(control_matric_total_stock, probs = 0.84, na.rm = TRUE)
    control_summary[counter, 'lower_soc'] = quantile(control_matric_total_stock, probs = 0.16, na.rm = TRUE)

    control_summary[counter, 'var_stat'] = median(control_matric_stat[control_matric_stat > 0], na.rm = TRUE)
    control_summary[counter, 'upper_stat'] = quantile(control_matric_stat[control_matric_stat > 0], probs = 0.84, na.rm = TRUE)
    control_summary[counter, 'lower_stat'] = quantile(control_matric_stat[control_matric_stat > 0], probs = 0.16, na.rm = TRUE)
    
    # control_summary[counter, 'var_soc'] = median(control_matric_total_stock[control_matric_stat > 0], na.rm = TRUE)
    # control_summary[counter, 'upper_soc'] = median(control_matric_total_stock[control_matric_stat > 0], na.rm = TRUE) + sd(control_matric_total_stock[control_matric_stat > 0], na.rm = TRUE)
    # control_summary[counter, 'lower_soc'] = median(control_matric_total_stock[control_matric_stat > 0], na.rm = TRUE) - sd(control_matric_total_stock[control_matric_stat > 0], na.rm = TRUE)
    # 
    # control_summary[counter, 'var_stat'] = median(control_matric_stat[control_matric_stat > 0], na.rm = TRUE)
    # control_summary[counter, 'upper_stat'] = median(control_matric_stat[control_matric_stat > 0], na.rm = TRUE) + sd(control_matric_stat[control_matric_stat > 0], na.rm = TRUE)
    # control_summary[counter, 'lower_stat'] = median(control_matric_stat[control_matric_stat > 0], na.rm = TRUE) - sd(control_matric_stat[control_matric_stat > 0], na.rm = TRUE)
    
    control_summary[counter, 'depth'] = depth_name[iquantile]
    control_summary[counter, 'component'] = component_name[icomponent]
    
    if (icomponent == 1 & iquantile == 4) {
      estimates_proda  = control_matric_total_stock
    }
    
    if (icomponent == 2 & iquantile == 4) {
      estimates_A  = control_matric_total_stock
      print(mean(estimates_A - estimates_proda))
      print(quantile((estimates_A - estimates_proda), probs = c(0.16, 0.5, 0.84), na.rm = TRUE))
      print(sd(estimates_A - estimates_proda))
    }
    
    counter = counter + 1
    
  }
}

# # record from best guess model
# control_stock_best_guess = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_soc_total_stock_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(7), '_control_test.mat', sep = ''))
# control_stock_best_guess = control_stock_best_guess$var.data.middle
# 
# control_stat_best_guess = readMat(paste(data_dir_output, 'world_simulation_analyses/compontent_control_spatial_variation_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_0_', as.character(7), '.mat', sep = ''))
# control_stat_best_guess = control_stat_best_guess$var.data.middle
# 
# control_summary[which(control_summary$component == 'H_NN'), 'var_soc'] = control_stock_best_guess[ , 1]
# control_summary[which(control_summary$component == 'H_NN'), 'var_stat'] = control_stat_best_guess[ , 1, stat_matric_index]*100
write.csv(control_summary, paste(data_dir_output, 'world_simulation_analyses/compontent_control_summary_bootstrap.csv', sep = ','))

### plot figure component only
current_data_range = control_summary[which(control_summary$depth == 'D_full_depth' & control_summary$component == 'H_NN'), ]
current_data = control_summary[which(control_summary$depth == 'D_full_depth' & control_summary$component != 'Ga_Homo' & control_summary$component != 'Gb_Default' & control_summary$component != 'H_NN'), ]

abs(current_data_range$var_soc - current_data$var_soc)[1]/abs(current_data_range$var_soc - current_data$var_soc)[2:6]
mean(abs(current_data_range$var_soc - current_data$var_soc)[1]/abs(current_data_range$var_soc - current_data$var_soc)[2:6])

abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:6]
mean(abs(current_data_range$var_stat - current_data$var_stat)[1]/abs(current_data_range$var_stat - current_data$var_stat)[2:6])

color_scheme = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00')
line_label = c('Microbial CUE', 'Environmental modifers', 'Plant carbon inputs', 'Baseline decomposition', 'Vertical transport', 'Carbon input allocation')

p_importance =
  ggplot(data = current_data) + 
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = current_data_range$lower_stat, ymax =  current_data_range$upper_stat, fill = 'grey45', alpha = 0.5) +
  annotate('rect', ymin = -Inf, ymax = Inf, xmin = current_data_range$lower_soc, xmax =  current_data_range$upper_soc, fill = 'grey45', alpha = 0.5) +
  geom_hline(yintercept = current_data_range$var_stat, size = 2, color = 'grey45', alpha = 1) +
  geom_vline(xintercept = current_data_range$var_soc, size = 2, color = 'grey45', alpha = 1) +
  geom_point(aes(x = var_soc, y = var_stat, color = component, fill = component), shape = 15, size = 10) + 
  geom_errorbar(aes(x = var_soc, ymin = lower_stat, ymax = upper_stat, color = component), size = 2, width = 40, stat = 'identity') +
  geom_errorbarh(aes(y = var_stat, xmin = lower_soc, xmax = upper_soc, color = component), size = 2, height = 0.4, stat = 'identity') +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  # scale_x_continuous(limits = c(1750, 2800)) +
  theme_classic() +
  # change the legend properties
  theme(legend.position = 'None') +
  # add title
  labs(x = 'Global SOC stock (Pg C)', y = expression(paste('Explained spatial variation (%)'), sep = '')) +
  # modify the position of title
  # modify the font size
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


#################################################################################
# sensitivity curve
#################################################################################
process_list =  c('A', 'K', 'Xi', 'V', 'I', 'NPP')
process_label = c('A_A', 'D_K', 'B_Xi', 'E_V', 'F_I', 'C_NPP')

process_description = c('CUE', 'Baseline Decomposition', 'Environmental Impacts', 'Vertical Transport', 'Input Allocation', 'Carbon Input')
climate_list = c('A', 'B', 'C', 'D', 'E_all')

manage_proposal = seq(from = -0.2, to = 0.2, by = 0.04)*100
depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')

grid_lon_lat = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_', process_list[1], '_', model_name, '_', nn_exp_name, '_bootstrap_', as.character(1), '.mat', sep = ''))

grid_lon_lat = grid_lon_lat$var.data.middle
grid_lon_lat = grid_lon_lat[[1]][[1]][ , c(1:2)]


############################# global soc stock
resolution = 0.5
lat_seq = grid_lon_lat[ , 2]
# area of grid 
radius = 6371008.8
length_top = (2*pi*radius*cos(abs(lat_seq+resolution/2)/180*pi)/360)*resolution
length_down = (2*pi*radius*cos(abs(lat_seq-resolution/2)/180*pi)/360)*resolution
height = (pi*radius/180)*resolution
lat_grid_area = (length_top + length_down)*height/2

marginal_change_total = array(NA, dim = c(3, length(process_list), length(manage_proposal), bootstrap_num))

iprocess = 1
counter = 1
for (iprocess in 1:length(process_list)) {
  
  bulk_process = array(NA, dim = c(nrow(grid_lon_lat), length(manage_proposal), bootstrap_num))
  soc_stock_tau = array(NA, dim = c(nrow(grid_lon_lat), 2, length(manage_proposal), bootstrap_num))
  
  ibootstrap = 1
  for (ibootstrap in 1:bootstrap_num) {
    print(paste('processing process ', iprocess, ' bootstrap ', ibootstrap, sep = ''))
    global_simu = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_', as.character(process_list[iprocess]), '_', model_name, '_', nn_exp_name, '_bootstrap_', as.character(ibootstrap), '.mat', sep = ''))
    global_simu = global_simu$var.data.middle
    
    soc_stock_tau[ , 1, , ibootstrap] = Re(global_simu[[2]][[1]])# soc stock
    soc_stock_tau[ , 2, , ibootstrap] = Re(global_simu[[3]][[1]]) # soc tau
    bulk_process[ , , ibootstrap] = Re(global_simu[[3+iprocess]][[1]]) # corresponding managed process
    
    imanage = 1
    for (imanage in 1:length(manage_proposal)) {
      invalid_loc = which(soc_stock_tau[ , 1, imanage, ibootstrap] < 0
                          | soc_stock_tau[ , 1, imanage, ibootstrap] > 1000000
                          | soc_stock_tau[ , 2, imanage, ibootstrap] < 0
                          | soc_stock_tau[ , 2, imanage, ibootstrap] > 100000
                          | bulk_process[ , imanage, ibootstrap] < 0)
      soc_stock_tau[invalid_loc, , , ] = NA
      bulk_process[invalid_loc, , ] = NA
      ### calculation of soc stock and tau
      #changes of soc stock
      marginal_change_total[1, iprocess, imanage, ibootstrap] = sum(soc_stock_tau[ , 1, imanage, ibootstrap]*lat_grid_area, na.rm = TRUE)/(10**15) 
      # changes of tau
      marginal_change_total[2, iprocess, imanage, ibootstrap] = mean(soc_stock_tau[ , 2, imanage, ibootstrap], na.rm = TRUE)
      # changes of bulk process
      marginal_change_total[3, iprocess, imanage, ibootstrap] = mean((bulk_process[ , imanage, ibootstrap] - bulk_process[ , 6, ibootstrap])/bulk_process[ , 6, ibootstrap], na.rm = TRUE)
    }
    
  }
}


####################################

marginal_change_total_mean = apply(marginal_change_total, c(1, 2, 3), mean, na.rm = TRUE)
marginal_change_total_upper = apply(marginal_change_total, c(1, 2, 3), quantile, probs = 0.16, na.rm = TRUE)
marginal_change_total_lower = apply(marginal_change_total, c(1, 2, 3), quantile, probs = 0.84, na.rm = TRUE)

col_list = c('stock', 'tau', 'stock_upper', 'stock_lower', 'tau_upper', 'tau_lower', 'process', 'process_upper', 'process_lower', 'process_name')
marginal_change_summary = data.frame(array(NA, dim = c(length(process_list)*length(manage_proposal), length(col_list))))
colnames(marginal_change_summary) = col_list

counter = 1

for (iprocess in 1:length(process_list)) {
  
  for (imanage in 1:length(manage_proposal)) {
    marginal_change_summary[counter, 'process_name'] = process_label[iprocess]
    marginal_change_summary[counter, 'stock'] = marginal_change_total_mean[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_upper'] = marginal_change_total_upper[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_lower'] = marginal_change_total_lower[1, iprocess, imanage]
    
    marginal_change_summary[counter, 'tau'] = marginal_change_total_mean[2, iprocess, imanage]
    marginal_change_summary[counter, 'tau_upper'] = marginal_change_total_upper[2, iprocess, imanage]
    marginal_change_summary[counter, 'tau_lower'] = marginal_change_total_lower[2, iprocess, imanage]
    
    marginal_change_summary[counter, 'process'] = marginal_change_total_mean[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_upper'] = marginal_change_total_upper[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_lower'] = marginal_change_total_lower[3, iprocess, imanage]
    
    counter = counter + 1
  }
}


write.csv(marginal_change_summary, paste(data_dir_output, 'world_simulation_analyses/marginal_change_summary_bootstrap.csv', sep = ','))

# marginal_change_summary = read.table(paste(data_dir_output, 'world_simulation_analyses/marginal_change_summary_bootstrap.csv', sep = ','))

color_scheme = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00')
line_label = c('Microbial CUE', 'Environmental modifers', 'Plant carbon inputs', 'Baseline decomposition', 'Vertical transport', 'Carbon input allocation')

p_curve =
  ggplot(data = marginal_change_summary) +
  geom_ribbon(aes(x = process*100, ymin = stock_lower, ymax = stock_upper, fill = process_name), alpha = 0.15) +
  geom_line(aes(x = process*100, y = stock, color = process_name), alpha = 1, size = 2) +
  geom_hline(yintercept = marginal_change_summary$stock[6], color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  # geom_hline(yintercept = marginal_change_summary$stock[6]*(1+0.004)**30, color = 'red', size = 2, linetype = 'dashed', alpha = 1) +
  geom_hline(yintercept = marginal_change_summary$stock[6]*(1+0.1), color = 'red', size = 2, linetype = 'dashed', alpha = 1) +
  geom_vline(xintercept = 0, color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  scale_x_continuous(trans = 'identity', limits = c(-10, 10), n.breaks = 7) +
  scale_y_continuous(trans = 'identity', limits = c(1500, 3200), n.breaks = 7) +
  # annotate('text', x = -10, y = marginal_change_summary$stock[6]*(1+0.004)**30, label = '"4per1000" by 2050 target', hjust = 0, vjust = -0.5, size = 10, color = 'red') +
  annotate('text', x = -10, y = marginal_change_summary$stock[6]*(1+0.1), label = '10% increase of total stock', hjust = 0, vjust = -0.5, size = 12, color = 'red') +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') + 
  theme(legend.justification = 'right', legend.position = 'right', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.8, 'inch')) +
  # add title
  labs(y = 'Global SOC stock (Pg C)', x = paste('Proportional changes of mechanisms (%)')) +
  # modify the position of title
  # modify the font sizea
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

jpeg(paste('./Ensemble/revision1_main_fig4.jpeg', sep = ''), width = 26.5, height = 10, units = 'in', res = 300)
plot_grid(p_importance, p_curve,
          labels = c('a', 'b'),
          rel_widths = c(1, 1.65),
          label_size = 60,
          label_x = -0.01, label_y = 1,
          label_fontfamily = 'Arial',
          label_fontface = 'bold',
          nrow = 1)
dev.off()



###########################################
col_list = c('contribution_soc', 'contribution_soc_std', 'contribution_tau', 'contribution_tau_std', 'process_name')
contribution_summary = data.frame(array(NA, dim = c(length(process_list), length(col_list))))
colnames(contribution_summary) = col_list

contribution_middle = array(NA, dim = c(2, length(process_list), bootstrap_num))
double_soc_tau = array(NA, dim = c(2, length(process_list), bootstrap_num))
fourper1000_soc_tau = array(NA, dim = c(2, length(process_list), bootstrap_num))

ibootstrap = 1
for (ibootstrap in 1:bootstrap_num) {
  iprocess = 1
  for (iprocess in 1:length(process_list)) {
    
    regression_data_soc = marginal_change_total[1, iprocess, , ibootstrap]
    regression_data_tau = marginal_change_total[2, iprocess, , ibootstrap]
    regression_data_process = marginal_change_total[3, iprocess, , ibootstrap]
    
    lm_soc = lm(log(regression_data_soc/regression_data_soc[6]) ~ log(regression_data_process + 1))
    contribution_middle[1, iprocess, ibootstrap] = regression_data_soc[6]*(1.01**lm_soc$coefficients[2] - 1)
    double_soc_tau[1, iprocess, ibootstrap] = exp((log(2) - lm_soc$coefficients[1])/lm_soc$coefficients[2]) - 1
    # fourper1000_soc_tau[1, iprocess, ibootstrap] = exp((log((1.004)**30) - lm_soc$coefficients[1])/lm_soc$coefficients[2]) - 1
    fourper1000_soc_tau[1, iprocess, ibootstrap] = exp((log((1 + 0.1)) - lm_soc$coefficients[1])/lm_soc$coefficients[2]) - 1
    
    lm_tau = lm(log(regression_data_tau/regression_data_tau[6]) ~ log(regression_data_process + 1))
    contribution_middle[2, iprocess, ibootstrap] = regression_data_tau[6]*(1.01**lm_tau$coefficients[2] - 1)
    double_soc_tau[2, iprocess, ibootstrap] = exp(lm_tau$coefficients[2]*log(double_soc_tau[1, iprocess, ibootstrap] + 1) + lm_tau$coefficients[1])
    fourper1000_soc_tau[2, iprocess, ibootstrap] = exp(lm_tau$coefficients[2]*log(fourper1000_soc_tau[1, iprocess, ibootstrap] + 1) + lm_tau$coefficients[1])
  }
  
  contribution_middle[1, , ibootstrap] = contribution_middle[1, , ibootstrap] # /sum(contribution_middle[1, , ibootstrap])*regression_data_soc[6]
  contribution_middle[2, , ibootstrap] = contribution_middle[2, , ibootstrap] # /sum(contribution_middle[2, , ibootstrap])*regression_data_tau[6]
}

apply(double_soc_tau, c(1, 2), mean, na.rm = TRUE)
apply(double_soc_tau, c(1, 2), sd, na.rm = TRUE)

apply(fourper1000_soc_tau, c(1, 2), mean, na.rm = TRUE)
apply(fourper1000_soc_tau, c(1, 2), sd, na.rm = TRUE)



