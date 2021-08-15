## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(viridis)
dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R/Ensemble')

####################################################################################
# Load data
####################################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
exp_name = 'exp_pc_cesm2_23' # 'exp_server_2019_median_5_3'
time_domain = 'whole_time'

data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'

current_data = c()
icross_valid = 1
for (icross_valid in 1:10) { 
  
  data_middle = readMat(paste(data_dir_output, 'world_simulation_analyses/indi_project_', model_name, '_', time_domain, '_', exp_name, '_cross_valid_', icross_valid, '.mat', sep = ''))
  data_middle = data_middle$indi.project
  current_data = rbind(current_data, data_middle)
}

current_data = data.frame(current_data)/1000 # unit kg C /m3
colnames(current_data) = c('obs', 'default', 'da', 'proda')


ss_tot = sum((current_data$obs - mean(current_data$obs, na.rm = TRUE))**2, na.rm = TRUE)
ss_res = sum((current_data$proda - current_data$obs)**2, na.rm = TRUE)
coeff_efficiency_proda = 1 - ss_res/ss_tot;
coeff_efficiency_proda

ss_tot = sum((current_data$obs - mean(current_data$obs, na.rm = TRUE))**2, na.rm = TRUE)
ss_res = sum((current_data$da - current_data$obs)**2, na.rm = TRUE)
coeff_efficiency_da = 1 - ss_res/ss_tot;
coeff_efficiency_da

ss_tot = sum((current_data$obs - mean(current_data$obs, na.rm = TRUE))**2, na.rm = TRUE)
ss_res = sum((current_data$default - current_data$obs)**2, na.rm = TRUE)
coeff_efficiency_default = 1 - ss_res/ss_tot;
coeff_efficiency_default


valid_loc = which(current_data$proda > 0.05 & current_data$da > 0.05 & current_data$default > 0.05 & current_data$obs > 0.05 &
                    is.na(current_data$proda) == 0 & is.na(current_data$da) == 0 & is.na(current_data$default) == 0 & is.na(current_data$obs) == 0)
current_data = current_data[valid_loc, ]

####################################################################################
# SOC Obs v.s. Mod
####################################################################################
## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

p_proda =
  ggplot(data = current_data) + 
  stat_bin_hex(aes(x = obs, y = proda), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(17), trans = 'identity', limits = c(1, 300), oob = scales::squish) + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + 
  annotate(geom = 'text', x = 0.1, y = 300, label = paste('Explained Variation: 56%'), hjust = 0, size = 10) + 
  
  # axis limit
  coord_cartesian(ylim = c(0.05, 500), xlim = c(0.05, 500)) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'black') + 
  # change the background to black and white
  theme_classic() +
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background = element_rect(fill = NA)) +
  # add title
  labs(title = 'PRODA', x = expression(paste('Observations (kg C m'^'-3', ')', sep = '')), y = expression(paste('Prediction (kg C m'^'-3', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
  # modify the font size
  theme(axis.title = element_text(size = 30)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


p_da = ggplot(data = current_data) + 
  stat_bin_hex(aes(x = obs, y = da), bins = 100) +
  scale_fill_gradientn(colors = viridis(17), trans = 'identity', limits = c(1, 300), oob = scales::squish) + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + 
  annotate(geom = 'text', x = 0.1, y = 300, label = paste('Explained Variation: 33%'), hjust = 0, size = 10) + 
  
  # axis limit
  coord_cartesian(ylim = c(0.05, 500), xlim = c(0.05, 500)) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'black') + 
  # change the background to black and white
  theme_classic() +
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 2, barheight = 7)) + 
  theme(legend.position = 'None') +
  # add title
  labs(title = 'Data Assimilation', x = expression(paste('Observations (kg C m'^'-3', ')', sep = '')), y = expression(paste('Prediction (kg C m'^'-3', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
  # modify the font size
  theme(axis.title = element_text(size = 30)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

p_default = ggplot(data = current_data) + 
  stat_bin_hex(aes(x = obs, y = default), bins = 100) +
  scale_fill_gradientn(colors = viridis(17), trans = 'identity', limits = c(1, 300), oob = scales::squish) + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + 
  annotate(geom = 'text', x = 0.1, y = 300, label = paste('Explained Variation: 11%'), hjust = 0, size = 10) + 
  
  # axis limit
  coord_cartesian(ylim = c(0.05, 500), xlim = c(0.05, 500)) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'black') + 
  # change the background to black and white
  theme_classic() +
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 2, barheight = 7)) + 
  theme(legend.position = 'None') +
  # add title
  labs(title = 'Default CLM5', x = expression(paste('Observations (kg C m'^'-3', ')', sep = '')), y = expression(paste('Prediction (kg C m'^'-3', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
  # modify the font size
  theme(axis.title = element_text(size = 30)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


jpeg('./obs_vs_mod.jpeg', width = 21, height = 7, units = 'in', res = 300)
plot_grid(p_default, p_da, p_proda, nrow = 1)
dev.off()

####################################################################################
# bulk process
####################################################################################
bulk_process_list = c('A', 'K', 'V', 'Xi', 'I')

process_scale_option = c('identity', 'log10', 'log10', 'identity', 'identity')

process_axis_label_x = c('Site Retrievals', 
                         expression(paste('Site Retrievals (yr'^'-1', ')', sep = '')), 
                         expression(paste('Site Retrievals (yr'^'-1', ')', sep = '')),
                         'Site Retrievals', 
                         'Site Retrievals')
process_axis_label_y = c('PRODA Predictions', 
                         expression(paste('PRODA Predictions (yr'^'-1', ')', sep = '')), 
                         expression(paste('PRODA Predictions (yr'^'-1', ')', sep = '')),
                         'PRODA Predictions', 
                         'PRODA Predictions')

process_name =  c('Microbial CUE', 
                  'Baseline Decomposition', 
                  'Vertical Transport',
                  'Environmental Modifers', 
                  'Input Allocation')

limit_lower = c(0.275, 0.005, 0.005,   0, 0.4)
limit_upper = c(0.6,  1,     0.2, 1, 1)

label_x = c(0.3,  0.01, 0.005, 0.05, 0.4)
label_y = c(0.55, 0.7,  0.15, 0.8, 0.9)

iprocess = 1
for (iprocess in 1:length(bulk_process_list)) {
  current_data = c()
  for (icross_valid in 1:10) { 
    data_middle = readMat(paste(data_dir_output, 'world_simulation_analyses/bulk_process_site_da_vs_proda_', model_name, '_', time_domain, '_', exp_name, '_cross_valid_', icross_valid, '.mat', sep = ''))
    data_middle = data_middle$bulk.process.site.da.vs.proda
    current_data = rbind(current_data, cbind(data_middle[[2]][ , iprocess], data_middle[[1]][ , iprocess]))
  }
  current_data = data.frame(current_data)
  colnames(current_data) = c('obs', 'proda')
  
  ss_tot = sum((current_data$obs - mean(current_data$obs, na.rm = TRUE))**2, na.rm = TRUE)
  ss_res = sum((current_data$proda - current_data$obs)**2, na.rm = TRUE)
  coeff_efficiency_proda = 1 - ss_res/ss_tot;
  coeff_efficiency_proda
  
  corr_test = cor.test(current_data$obs, current_data$proda)
  
  valid_loc = which(current_data$obs > limit_lower[iprocess] & current_data$obs < limit_upper[iprocess]
                    & current_data$proda > limit_lower[iprocess] & current_data$proda < limit_upper[iprocess])
  current_data = current_data[valid_loc, ]
  
  p = 
    ggplot(data = current_data) + 
    stat_bin_hex(aes(x = obs, y = proda), bins = 100) +
    scale_fill_gradientn(name = 'Count', colors = viridis(17), trans = 'identity', limits = c(1, 100), oob = scales::squish) + 
    scale_x_continuous(trans = process_scale_option[iprocess]) + 
    scale_y_continuous(trans = process_scale_option[iprocess]) + 
    annotate(geom = 'text', x = label_x[iprocess], y = label_y[iprocess], label = paste('r = ', round(corr_test$estimate, 2), '\nP < 0.001'), hjust = 0, size = 10) +
    
    # axis limit
    coord_cartesian(ylim = c(limit_lower[iprocess], limit_upper[iprocess]), c(limit_lower[iprocess], limit_upper[iprocess])) +
    geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'black') + 
    # change the background to black and white
    theme_classic() +
    # change the legend properties
    guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
    theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
    theme(legend.justification = c(1, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) +
    # add title
    labs(title = process_name[iprocess], x = process_axis_label_x[iprocess], y = process_axis_label_y[iprocess], sep = '') + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
    # modify the font size
    theme(axis.title = element_text(size = 30)) + 
    # modify the margin
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
    theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 
  
  
  eval(parse(text = paste('p', iprocess, ' = p', sep = '')))
}


p_proda = p_proda + 
  labs(title = 'SOC', x = expression(paste('Observations (kg C m'^'-3', ')', sep = '')), y = expression(paste('PRODA Predictions (kg C m'^'-3', ')', sep = ''))) + 
  theme(legend.justification = c(1, 0), legend.position = 'None', legend.background = element_rect(fill = NA))


jpeg('./bulk_process_site_da_vs_proda.jpeg', width = 21, height = 14, units = 'in', res = 300)
plot_grid(p1, p2, p3, p4, p5, p_proda, nrow = 2)
dev.off()



