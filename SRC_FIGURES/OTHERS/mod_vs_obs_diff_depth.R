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
coeff_efficiency_cross_valid = array(NA, dim = c(10, 1))
icross_valid = 10
for (icross_valid in 1:10) { 
  
  # data_middle = readMat(paste(data_dir_output, 'world_simulation_analyses/nn_soc_validation_auto_corr_', model_name, '_', time_domain, '_', exp_name, '_cross_valid_', icross_valid, '.mat', sep = ''))
  data_middle = readMat(paste(data_dir_output, 'world_simulation_analyses/nn_soc_validation_', model_name, '_', time_domain, '_', exp_name, '_cross_valid_', icross_valid, '.mat', sep = ''))
  data_middle = data_middle$nn.soc.validation
  
  data_middle = cbind(as.vector(data_middle[[2]]), as.vector(data_middle[[1]]), as.vector(data_middle[[3]]))
  data_middle = data_middle[which(is.na(data_middle[ , 1]) == 0), ]
  current_data = rbind(current_data, data_middle)
  
  ss_tot = sum((data_middle[ , 1] - mean(data_middle[ , 1], na.rm = TRUE))**2, na.rm = TRUE)
  ss_res = sum((data_middle[ , 2] - data_middle[ , 1])**2, na.rm = TRUE)
  coeff_efficiency_cross_valid[icross_valid] = 1 - ss_res/ss_tot;

}
mean(coeff_efficiency_cross_valid, na.rm = TRUE)
sd(coeff_efficiency_cross_valid, na.rm = TRUE)

current_data[ , c(1, 2)] = current_data[ , c(1, 2)]/1000 # unit kg C /m3
current_data = data.frame(current_data)
colnames(current_data) = c('obs', 'proda', 'depth')


ss_tot = sum((current_data$obs - mean(current_data$obs, na.rm = TRUE))**2, na.rm = TRUE)
ss_res = sum((current_data$proda - current_data$obs)**2, na.rm = TRUE)
coeff_efficiency_proda = 1 - ss_res/ss_tot;
coeff_efficiency_proda


ss_tot = sum((current_data$obs[current_data$depth <= 0.3] - mean(current_data$obs[current_data$depth <= 0.3], na.rm = TRUE))**2, na.rm = TRUE)
ss_res = sum((current_data$proda[current_data$depth <= 0.3] - current_data$obs[current_data$depth <= 0.3])**2, na.rm = TRUE)
coeff_efficiency_proda = 1 - ss_res/ss_tot;
coeff_efficiency_proda

ss_tot = sum((current_data$obs[current_data$depth > 0.3 & current_data$depth <= 1] - mean(current_data$obs[current_data$depth > 0.3 & current_data$depth <= 1], na.rm = TRUE))**2, na.rm = TRUE)
ss_res = sum((current_data$proda[current_data$depth > 0.3 & current_data$depth <= 1] - current_data$obs[current_data$depth > 0.3 & current_data$depth <= 1])**2, na.rm = TRUE)
coeff_efficiency_proda = 1 - ss_res/ss_tot;
coeff_efficiency_proda

ss_tot = sum((current_data$obs[current_data$depth > 1] - mean(current_data$obs[current_data$depth > 1], na.rm = TRUE))**2, na.rm = TRUE)
ss_res = sum((current_data$proda[current_data$depth > 1] - current_data$obs[current_data$depth > 1])**2, na.rm = TRUE)
coeff_efficiency_proda = 1 - ss_res/ss_tot;
coeff_efficiency_proda


valid_loc = which(current_data$proda > 0.05 &
                    is.na(current_data$proda) == 0 & is.na(current_data$obs) == 0)
current_data = current_data[valid_loc, ]

####################################################################################
# SOC Obs v.s. Mod
####################################################################################
## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

p_30cm =
  ggplot(data = current_data[current_data$depth <= 0.3, ]) + 
  stat_bin_hex(aes(x = obs, y = proda), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(17), trans = 'identity', limits = c(1, 300), oob = scales::squish) + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + 
  annotate(geom = 'text', x = 0.1, y = 300, label = paste('Explained Variation: 43%'), hjust = 0, size = 10) + 
  
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
  labs(title = '0 - 30cm', x = expression(paste('Observations (kg C m'^'-3', ')', sep = '')), y = expression(paste('Prediction (kg C m'^'-3', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
  # modify the font size
  theme(axis.title = element_text(size = 30)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 


p_100cm =
  ggplot(data = current_data[current_data$depth > 0.3 & current_data$depth <= 1, ]) + 
  stat_bin_hex(aes(x = obs, y = proda), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(17), trans = 'identity', limits = c(1, 300), oob = scales::squish) + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + 
  annotate(geom = 'text', x = 0.1, y = 300, label = paste('Explained Variation: 24%'), hjust = 0, size = 10) + 
  
  # axis limit
  coord_cartesian(ylim = c(0.05, 500), xlim = c(0.05, 500)) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'black') + 
  # change the background to black and white
  theme_classic() +
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  # add title
  labs(title = '30 - 100cm', x = expression(paste('Observations (kg C m'^'-3', ')', sep = '')), y = expression(paste('Prediction (kg C m'^'-3', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
  # modify the font size
  theme(axis.title = element_text(size = 30)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 

p_200cm =
  ggplot(data = current_data[current_data$depth > 1, ]) + 
  stat_bin_hex(aes(x = obs, y = proda), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(17), trans = 'identity', limits = c(1, 300), oob = scales::squish) + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') + 
  annotate(geom = 'text', x = 0.1, y = 300, label = paste('Explained Variation: 11%'), hjust = 0, size = 10) + 
  
  # axis limit
  coord_cartesian(ylim = c(0.05, 500), xlim = c(0.05, 500)) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'black') + 
  # change the background to black and white
  theme_classic() +
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  # add title
  labs(title = '> 100cm', x = expression(paste('Observations (kg C m'^'-3', ')', sep = '')), y = expression(paste('Prediction (kg C m'^'-3', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 35)) + 
  # modify the font size
  theme(axis.title = element_text(size = 30)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 


jpeg('./obs_vs_mod_proda_depth.jpeg', width = 21, height = 7, units = 'in', res = 300)
plot_grid(p_30cm, p_100cm, p_200cm, nrow = 1)
dev.off()

