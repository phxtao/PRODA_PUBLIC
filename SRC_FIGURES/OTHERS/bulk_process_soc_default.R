## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(jcolors)
library(viridis)
library(cowplot)

dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))

#############################################################################
# Data Path
#############################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

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
# Load Projected SOC site by site
#################################################################################
grid_var_names = c('Lon', 'Lat', 'Date', 
                   'Rmean', 'Rmax', 'Rmin', 
                   'ESA_Land_Cover', 
                   'ET',
                   'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin',
                   'Veg_Cover', 
                   'Annual Mean Temperature', 'Mean Diurnal Range', 'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 
                   'Abs_Depth_to_Bedrock',
                   'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',
                   'CEC_0cm', 'CEC_30cm', 'CEC_100cm',
                   'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',
                   'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', 
                   'Depth_Bedrock_R', 
                   'Garde_Acid', 
                   'Occurrence_R_Horizon', 
                   'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', 
                   'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm',
                   'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', 
                   'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', 
                   'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', 
                   'USDA_Suborder', 
                   'WRB_Subgroup', 
                   'Drought',
                   'Elevation',
                   'Max_Depth', 
                   'Koppen_Climate_2018', 
                   'cesm2_npp', 'cesm2_npp_std',
                   'cesm2_gpp', 'cesm2_gpp_std',
                   'cesm2_vegc',
                   'nbedrock')

#################################################################################
# Load Projected SOC default
#################################################################################
bulk_process_list = c('A', 'I', 'K', 'V', 'Xi', 'NPP')

global_lat_lon = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_default_NPP_', model_name, '.mat', sep = ''))

global_lat_lon = global_lat_lon$Var.Decom.Grid[[1]]
global_lat_lon = global_lat_lon[ , 1:2]


bulk_process = array(NA, dim = c(nrow(global_lat_lon), length(bulk_process_list)))
soc_stock_tau = array(NA, dim = c(nrow(global_lat_lon), 2))

global_simu = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_default_NPP_', model_name, '.mat', sep = ''))

bulk_process_A = Re(global_simu$Var.Decom.Grid[[11]])[ , 6]
bulk_process_K = Re(global_simu$Var.Decom.Grid[[12]])[ , 6]
bulk_process_Xi = Re(global_simu$Var.Decom.Grid[[13]])[ , 6]
bulk_process_V = Re(global_simu$Var.Decom.Grid[[14]])[ , 6]
bulk_process_I = Re(global_simu$Var.Decom.Grid[[15]])[ , 6]
bulk_process_NPP = Re(global_simu$Var.Decom.Grid[[16]])[ , 6]

global_simu_process = cbind(bulk_process_A, bulk_process_I, bulk_process_K, bulk_process_V, bulk_process_Xi, bulk_process_NPP)
global_simu_soc = cbind(Re(global_simu$Var.Decom.Grid[[5]])[ , 6], Re(global_simu$Var.Decom.Grid[[9]])[ , 6])

soc_content_0_30cm = Re(global_simu$Var.Decom.Grid[[2]])[ , 6]
soc_tau_0_30cm = Re(global_simu$Var.Decom.Grid[[6]])[ , 6]

invalid_loc = which(global_simu_soc[ , 1] < 0
                    | global_simu_soc[ , 1] > 1000000
                    | global_simu_soc[ , 2] < 0
                    | global_simu_soc[ , 2] > 100000
                    # no process should less than 0 
                    | global_simu_process[ , 1] < 0
                    | global_simu_process[ , 2] < 0
                    | global_simu_process[ , 3] < 0
                    | global_simu_process[ , 4] < 0
                    | global_simu_process[ , 5] < 0
                    # CUE never larger than 1
                    | global_simu_process[ , 1] > 1
)

global_simu_soc[invalid_loc, ] = NA
global_simu_process[invalid_loc, ] = NA

soc_content_0_30cm[invalid_loc] = NA
soc_tau_0_30cm[invalid_loc] = NA

# global_simu_process[ , 3] = global_simu_process[ , 3]*global_simu_soc[ , 1]
# global_simu_process[ , 4] = global_simu_process[ , 4]*global_simu_soc[ , 1]

bulk_process[ , ] = global_simu_process
soc_stock_tau[ , ] = global_simu_soc


valid_grid_loc = read.csv(paste(data_dir_output, 'neural_networking/valid_grid_loc_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_', as.character(1), '.csv', sep = ''), header = FALSE)
valid_grid_loc = valid_grid_loc$V1
grid_env_info = readMat(paste(data_dir_input, 'data4nn/world_grid_envinfo_present.mat', sep = ''))
grid_env_info = grid_env_info$EnvInfo
colnames(grid_env_info) = grid_var_names


soc_content_0_30cm = soc_content_0_30cm/0.3/grid_env_info[valid_grid_loc, 'Bulk_Density_0cm'] # unit gC/kg

#################################################################################
# Plot Figures Residence time and cue (0 - 30cm)
#################################################################################
current_data = data.frame(cbind(global_simu_process[ , 1], soc_tau_0_30cm, grid_env_info[valid_grid_loc, 'Annual Mean Temperature']))
colnames(current_data) = c('cue', 'tau', 'mat')

lm_model = lm(tau ~ cue, data = current_data)
lm_summary = summary(lm_model)
lm_summary

corr_test = cor.test(current_data$cue, current_data$tau)
corr_test$estimate


text_data = data.frame('x_axis' = c(0.525), 
                       'y_axis' = c(1), 
                       'equation' = paste('r = ', round(corr_test$estimate, 2), '\np < 0.001', sep = ''))

jpeg(paste('./Ensemble/default_cue_tau_0_30cm.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

ggplot(data = current_data) + 
  stat_bin_hex(aes(x = cue, y = tau), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(limits = c(0.25, 0.65), trans = 'identity') +
  scale_y_continuous(limits = c(0.3, 3000), trans = 'log10') + 
  geom_smooth(data = current_data, aes(x = cue, y = tau), method = 'lm', color = 'black', fill = 'blue3', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'CUE', y =expression(paste('Residence Time (0 - 30cm, year)', sep = ''))) + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

dev.off()


#################################################################################
# Plot Figures SOC content and cue (0 - 30cm)
#################################################################################
current_data = data.frame(cbind(global_simu_process[ , 1][ , 1], soc_content_0_30cm, grid_env_info[valid_grid_loc, 'Annual Mean Temperature']))
colnames(current_data) = c('cue', 'soc', 'mat')

lm_model = lm(soc ~ cue, data = current_data)
lm_summary = summary(lm_model)
lm_summary

corr_test = cor.test(current_data$cue, current_data$soc)
corr_test$estimate


text_data = data.frame('x_axis' = c(0.525), 
                       'y_axis' = c(1), 
                       'equation' = paste('r = ', round(corr_test$estimate, 2), '\np < 0.001', sep = ''))


jpeg(paste('./Ensemble/default_cue_soc_0_30cm.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

ggplot(data = current_data) + 
  stat_bin_hex(aes(x = cue, y = soc), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(limits = c(0.25, 0.65), trans = 'identity') +
  scale_y_continuous(limits = c(0.3, 600), trans = 'log10') + 
  geom_smooth(data = current_data, aes(x = cue, y = soc), method = 'lm', color = 'black', fill = 'blue3', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'CUE', y =expression(paste('SOC (0 - 30cm g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

dev.off()

#################################################################################
# Plot Figures CUE vs. MAT
#################################################################################

lm_model = lm(cue ~ mat, data = current_data)
lm_summary = summary(lm_model)

corr_test = cor.test(current_data$cue, current_data$mat)


text_data = data.frame('x_axis' = c(15), 
                       'y_axis' = c(0.6), 
                       'equation' = paste('r = ', round(corr_test$estimate, 2), '\np < 0.001', sep = ''))

jpeg(paste('./Ensemble/default_mat_cue.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

ggplot(data = current_data) + 
  stat_bin_hex(aes(x = mat, y = cue), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(trans = 'identity') + 
  scale_y_continuous(limits = c(0.25, 0.65), trans = 'identity') +
  geom_smooth(data = current_data, aes(x = mat, y = cue), method = 'lm', color = 'black', fill = 'blue3', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'MAT (deg C)', y = 'CUE') + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

dev.off()
