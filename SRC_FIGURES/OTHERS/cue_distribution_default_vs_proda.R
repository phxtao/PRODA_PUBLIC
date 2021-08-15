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

cue_default = Re(global_simu$Var.Decom.Grid[[11]])[ , 6]

#################################################################################
# Load Projected SOC PRODA
#################################################################################
cross_valid_num = 10

bulk_process_list = c('A', 'I', 'K', 'V', 'Xi', 'NPP')

global_lat_lon = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_NPP_', model_name, '_', nn_exp_name, '_cross_valid_', as.character(1), '.mat', sep = ''))

global_lat_lon = global_lat_lon$Var.Decom.Grid[[1]]
global_lat_lon = global_lat_lon[ , 1:2]


bulk_process = array(NA, dim = c(nrow(global_lat_lon), length(bulk_process_list), cross_valid_num))
soc_stock_tau = array(NA, dim = c(nrow(global_lat_lon), 2, cross_valid_num))

icross_valid = 1
for (icross_valid in 1:cross_valid_num) {
  global_simu = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_NPP_', model_name, '_', nn_exp_name, '_cross_valid_', as.character(icross_valid), '.mat', sep = ''))
  
  bulk_process_A = Re(global_simu$Var.Decom.Grid[[11]])[ , 11]
  bulk_process_K = Re(global_simu$Var.Decom.Grid[[12]])[ , 11]
  bulk_process_Xi = Re(global_simu$Var.Decom.Grid[[13]])[ , 11]
  bulk_process_V = Re(global_simu$Var.Decom.Grid[[14]])[ , 11]
  bulk_process_I = Re(global_simu$Var.Decom.Grid[[15]])[ , 11]
  bulk_process_NPP = Re(global_simu$Var.Decom.Grid[[16]])[ , 11]
  
  global_simu_process = cbind(bulk_process_A, bulk_process_I, bulk_process_K, bulk_process_V, bulk_process_Xi, bulk_process_NPP)
  global_simu_soc = cbind(Re(global_simu$Var.Decom.Grid[[5]])[ , 11], Re(global_simu$Var.Decom.Grid[[9]])[ , 11])
  
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
  
  # global_simu_process[ , 3] = global_simu_process[ , 3]*global_simu_soc[ , 1]
  # global_simu_process[ , 4] = global_simu_process[ , 4]*global_simu_soc[ , 1]
  
  bulk_process[ , , icross_valid] = global_simu_process
  soc_stock_tau[ , , icross_valid] = global_simu_soc
  
}

bulk_process_mean = apply(bulk_process, c(1, 2), mean, na.rm = TRUE)

cue_proda = bulk_process_mean[ , 1]



#################################################################################
# distribution figures
#################################################################################
color_theme = c('#005AB5', '#DC3220')


current_data = rbind(cbind(cue_default, 'Default CLM5'), cbind(cue_proda, 'PRODA'))

current_data = data.frame(current_data)
colnames(current_data) = c('cue', 'source')
current_data$cue = as.numeric(current_data$cue)

jpeg(paste('./Ensemble/cue_distribution_default_vs_proda.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

ggplot() + 
  geom_histogram(data = current_data, aes(x = cue, y = ..density.., group = source, color = source, fill = source), alpha = 0.5, position = 'identity', bins = 30) +
  scale_color_manual(name = '', labels = c('Default CLM5', 'PRODA'), values = color_theme) +
  scale_fill_manual(name = '', labels = c('Default CLM5', 'PRODA'), values = color_theme) +
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', x = 'CUE', y = 'Density') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 
dev.off()









