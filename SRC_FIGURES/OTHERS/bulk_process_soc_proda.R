## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(jcolors)
library(gridExtra)
library(viridis)

library(GGally)

library(rgdal)
library(raster)
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
# Env. Var Names
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

soc_stock_tau_mean = apply(soc_stock_tau, c(1, 2), mean, na.rm = TRUE)
soc_stock_tau_std = apply(soc_stock_tau, c(1, 2), sd, na.rm = TRUE)

bulk_process_mean = apply(bulk_process, c(1, 2), mean, na.rm = TRUE)
bulk_process_std = apply(bulk_process, c(1, 2), sd, na.rm = TRUE)

## basic metrix
resolution = 0.5
lat_seq = global_lat_lon[ , 2]
# area of grid 
radius = 6371008.8
length_top = (2*pi*radius*cos(abs(lat_seq+resolution/2)/180*pi)/360)*resolution
length_down = (2*pi*radius*cos(abs(lat_seq-resolution/2)/180*pi)/360)*resolution
height = (pi*radius/180)*resolution
lat_grid_area = (length_top + length_down)*height/2
total_carbon_input = sum(bulk_process_mean[ , 6]*lat_grid_area, na.rm = TRUE)/10**15 # unit Pg/yr
mean_carbon_input = sum(bulk_process_mean[ , 6]*lat_grid_area/sum(lat_grid_area, na.rm = TRUE), na.rm = TRUE)

valid_grid_loc = read.csv(paste(data_dir_output, 'neural_networking/valid_grid_loc_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_', as.character(1), '.csv', sep = ''), header = FALSE)
valid_grid_loc = valid_grid_loc$V1
grid_env_info = readMat(paste(data_dir_input, 'data4nn/world_grid_envinfo_present.mat', sep = ''))
grid_env_info = grid_env_info$EnvInfo
colnames(grid_env_info) = grid_var_names

var_temp = grid_env_info[valid_grid_loc, 'Annual Mean Temperature']
var_prec = grid_env_info[valid_grid_loc, 'Annual Precipitation']
var_clay = apply(grid_env_info[valid_grid_loc, c('Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm')], 1, mean, na.rm = TRUE)
var_bd = grid_env_info[valid_grid_loc, 'Bulk_Density_0cm']*0.3/2 + grid_env_info[valid_grid_loc, 'Bulk_Density_0cm']*(1-0.3)/2 + grid_env_info[valid_grid_loc, 'Bulk_Density_100cm']*(2 - 1)/2

# load climate info
global_climate = grid_env_info[valid_grid_loc, 'Koppen_Climate_2018']
var_climate = array(NA, dim = c(length(global_climate), 1))

var_climate[which(global_climate >= 1 & global_climate <= 3)] = 'A'
var_climate[which(global_climate >= 4 & global_climate <= 7)] = 'B'
var_climate[which(global_climate >= 8 & global_climate <= 16)] = 'C'
var_climate[which(global_climate >= 17 & global_climate <= 30)] = 'D'

# load ESA land cover
global_landcover = grid_env_info[valid_grid_loc, 'ESA_Land_Cover']
var_landcover = array(NA, dim = c(length(global_landcover), 1))

var_landcover[which(global_landcover >= 1 & global_landcover <= 4)] = 'A' # agriculture
var_landcover[which(var_climate == 'D' & ((global_landcover >= 5 & global_landcover <= 10) | (global_landcover >= 16 & global_landcover <= 17)))] = 'B' # boreal forest
var_landcover[which(var_climate != 'D' & ((global_landcover >= 5 & global_landcover <= 10) | (global_landcover >= 16 & global_landcover <= 17)))] = 'C' # other forest
var_landcover[which(global_landcover >= 11 & global_landcover <= 15)] = 'D' # grassland & shrubland 
var_landcover[which(global_landcover == 18)] = 'E' # wetland
var_landcover[which(global_landcover >=19)] = 'F' # urban and other

round(apply(bulk_process_mean, 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)

round(apply(bulk_process_mean[var_landcover == 'A', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'A', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'B', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'B', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'C', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'C', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'D', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'D', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'E', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'E', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'F', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_landcover == 'F', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)

round(apply(bulk_process_mean[var_climate == 'A', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_climate == 'A', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_climate == 'B', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_climate == 'B', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_climate == 'C', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_climate == 'C', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_climate == 'D', ], 2, mean, na.rm = TRUE), digits = 2)
round(apply(bulk_process_mean[var_climate == 'D', ], 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE), digits = 2)

# slope info
slope_info = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/World_Slope_FAO_HWSD/slope_grid_half_deg.mat')
slope_info = slope_info$slope.info
slope_info = slope_info[valid_grid_loc, ]


#################################################################################
# Plot Figures
#################################################################################
bulk_process_list = c('A', 'I', 'K', 'V', 'Xi', 'NPP')

process_scale_option = c('identity', 'identity', 'log10', 'identity', 'identity', 'identity')
process_axis_label = c('Microbial CUE', 
                       'Input Allocation',
                       expression(paste('Baseline Decomposition (yr'^'-1', ')', sep = '')), 
                       expression(paste('Vertical Transport Rate (yr'^'-1', ')', sep = '')),
                       'Environmental Modifer', 
                       expression(paste('Input Carbon (g C yr'^'-1', ')', sep = '')))
process_name =  c('Microbial CUE', 
                  'Input Allocation', 
                  'Baseline Decomposition', 
                  'Vertical Transport Rate',
                  'Environmental Modifer', 
                  'Input Carbon')
process_unit = c('unitless',
                 'unitless', 
                 expression(paste('yr'^'-1', sep = '')),
                 expression(paste('yr'^'-1', sep = '')),
                 'unitless', 
                 expression(paste('g C yr'^'-1', sep = '')))


legend_limit_lower = c(0.3, 0.45, 0.008, 0.007, 0.0, 3) # apply(bulk_process_mean, 2, quantile, prob = 0.005, na.rm = TRUE)
legend_limit_upper = c(0.5, 0.85, 0.3,   0.14,  0.6, 1400) # apply(bulk_process_mean, 2, quantile, prob = 0.995, na.rm = TRUE)


legend_limit_lower_uncertain = c(0.,   0,   0.001,  0,     0.0,  0) # apply(bulk_process_mean, 2, quantile, prob = 0.005, na.rm = TRUE)
legend_limit_upper_uncertain = c(0.02, 0.2, 0.09,   0.02, 0.05, 100) # apply(bulk_process_mean, 2, quantile, prob = 0.995, na.rm = TRUE)

# legend_limit_lower = c(0.25, 0.45, 300,   300,   0.10, 3) # apply(bulk_process_mean, 2, quantile, prob = 0.005, na.rm = TRUE)
# legend_limit_upper = c(0.42, 0.85, 3000, 3000, 0.53, 1400) # apply(bulk_process_mean, 2, quantile, prob = 0.995, na.rm = TRUE)


temp_limit = c(quantile(var_temp, 0.005, na.rm = TRUE), quantile(var_temp, 0.995, na.rm = TRUE))
soc_stock_limit = c(quantile(soc_stock_tau_mean[ , 1], 0.005, na.rm = TRUE), quantile(soc_stock_tau_mean[ , 1], 0.995, na.rm = TRUE))/1000
soc_tau_limit = c(quantile(soc_stock_tau_mean[ , 2], 0.005, na.rm = TRUE), quantile(soc_stock_tau_mean[ , 2], 0.995, na.rm = TRUE))


##################################soc stock and Residence Time map
world_coastline = rgdal::readOGR(dsn='/Users/phoenix/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp',layer = 'ne_110m_land')
world_coastline <- spTransform(world_coastline, CRS('+proj=robin'))
world_coastline <- fortify(world_coastline)
Map.Using = world_coastline

ocean_left = cbind(rep(-180, 100), seq(from = 80, to = -56, by = -(80 + 56)/(100 -1)))
ocean_right = cbind(rep(180, 100), seq(from = -56, to = 80, by = (80 + 56)/(100 -1)))
ocean_top = cbind(seq(from = 180, to = -180, by = -(360)/(100 -1)), rep(80, 100))
ocean_bottom = cbind(seq(from = -180, to = 180, by = (360)/(100 -1)), rep(-56, 100))

world_ocean = rbind(ocean_left, ocean_bottom, ocean_right, ocean_top)
world_ocean = as.matrix(world_ocean)

world_ocean <- project(xy = world_ocean, proj = '+proj=robin') 

world_ocean = data.frame(world_ocean)
colnames(world_ocean) = c('long', 'lat')

color_scheme = c('#D55E00', '#E69F00', '#009E73', '#0072B2')

# soc map
CurrentData = data.frame(cbind(global_lat_lon, soc_stock_tau_mean, soc_stock_tau_std))
colnames(CurrentData) = c('Lon', 'Lat', 'soc', 'tau', 'std_soc', 'std_tau')
CurrentData$climate = var_climate

lon_lat_transfer = project(xy = as.matrix(CurrentData[ , c('Lon', 'Lat')]), proj = '+proj=robin') 
CurrentData[ , c('Lon', 'Lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_soc =
  ggplot() +
  geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = soc/1000), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(name = expression(paste('kg C m'^'-2', sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(10, 100), trans = 'log10', oob = scales::squish) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
  # change the background to black and white
  coord_equal() +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0.1, 0.36), legend.position = c(0.1, 0.36), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
  # add title
  labs(title = 'SOC Stock', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35))




p_std_soc =
  ggplot() +
  geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = std_soc/1000), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(name = expression(paste('kg C m'^'-2', sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(1, 50), trans = 'log10', oob = scales::squish) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
  # change the background to black and white
  coord_equal() +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0.1, 0.36), legend.position = c(0.1, 0.36), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
  # add title
  labs(title = 'Uncertainty SOC Stock', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35))


box_soc =
  ggplot() +
  geom_boxplot(size = 1, data = CurrentData, aes(x = climate, y = soc/1000, group = climate, fill = climate), outlier.colour = NA, outlier.shape = NA, outlier.size = NA) +
  scale_y_continuous(limits = c(2, 200), trans = 'log10', position = 'right') +
  scale_x_discrete(label = c('Tropical', 'Arid', 'Temperate', 'Boreal')) +
  scale_fill_manual(name = '', labels = '', values = color_scheme) +
  coord_flip() +
  # change the background to black and white
  theme_classic() +
  # change the legend properties
  theme(legend.position = 'none') +
  # add title
  labs(title = '', x = '', y = '') +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
  # modify the font size
  theme(axis.title = element_text(size = 30)) + 
  # modify the margin
  theme(plot.margin = unit(c(-0.4, 0.1, 0.5, -0.7), 'inch')) +
  theme(axis.text=element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +  
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch'))



p_tau =
  ggplot() +
  geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = tau), height = 60000, width = 60000, na.rm = TRUE) + 
  scale_fill_gradientn(name = 'year', colours = rev(viridis(15)), na.value="transparent", limits = c(10, 1000), trans = 'log10', oob = scales::squish) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
  # change the background to black and white
  coord_equal() +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0.1, 0.36), legend.position = c(0.1, 0.36), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
  # add title
  labs(title = 'SOC Residence Time', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35))
  


p_std_tau =
  ggplot() +
  geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = std_tau), height = 60000, width = 60000, na.rm = TRUE) + 
  scale_fill_gradientn(name = 'year', colours = rev(viridis(15)), na.value="transparent", limits = c(1, 500), trans = 'log10', oob = scales::squish) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
  # change the background to black and white
  coord_equal() +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0.1, 0.36), legend.position = c(0.1, 0.36), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
  # add title
  labs(title = 'Uncertainty Residence Time', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35))


box_tau =
  ggplot() +
  geom_boxplot(size = 1, data = CurrentData, aes(x = climate, y = tau, group = climate, fill = climate), outlier.colour = NA, outlier.shape = NA, outlier.size = NA) +
  scale_y_continuous(limits = c(2, 2000), trans = 'log10', position = 'right') +
  scale_x_discrete(label = c('Tropical', 'Arid', 'Temperate', 'Boreal')) +
  scale_fill_manual(name = '', labels = '', values = color_scheme) +
  coord_flip() +
  # change the background to black and white
  theme_classic() +
  # change the legend properties
  theme(legend.position = 'none') +
  # add title
  labs(title = '', x = '', y = '') +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
  # modify the font size
  theme(axis.title = element_text(size = 30)) + 
  # modify the margin
  theme(plot.margin = unit(c(-0.4, 0.1, 0.5, -0.7), 'inch')) +
  theme(axis.text=element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch'))



################################## process map

ipara = 1
## Projected para 
for (ipara in 1:length(process_name)){
  CurrentData = data.frame(cbind(global_lat_lon, bulk_process_mean[ , ipara]))
  colnames(CurrentData) = c('Lon', 'Lat', 'Project')
  CurrentData$climate = var_climate
  
  lon_lat_transfer = project(xy = as.matrix(CurrentData[ , c('Lon', 'Lat')]), proj = '+proj=robin') 
  CurrentData[ , c('Lon', 'Lat')] = lon_lat_transfer
  
  p =
    ggplot() +
    geom_boxplot(size = 1, data = CurrentData, aes(x = climate, y = Project, group = climate, fill = climate), outlier.colour = NA, outlier.shape = NA, outlier.size = NA) +
    scale_y_continuous(limits = c(legend_limit_lower[ipara], legend_limit_upper[ipara]), trans = process_scale_option[ipara], position = 'right') +
    scale_x_discrete(label = c('Tropical', 'Arid', 'Temperate', 'Boreal')) +
    scale_fill_manual(name = '', labels = '', values = color_scheme) +
    coord_flip() +
    # change the background to black and white
    theme_classic() +
    # change the legend properties
    theme(legend.position = 'none') +
    # add title
    labs(title = '', x = '', y = '') +
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
    # modify the font size
    theme(axis.title = element_text(size = 30)) + 
    # modify the margin
    theme(plot.margin = unit(c(-0.4, 0.1, 0.5, -0.7), 'inch')) +
    theme(axis.text=element_text(size = 30)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch'))
  
  
  
  eval(parse(text = paste('box_', ipara, ' = p', sep = '')))
  
  p =
    ggplot() +
    geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = Project), height = 60000, width = 60000, na.rm = TRUE) +
    scale_fill_gradientn(name = process_unit[ipara], colours = rev(viridis(15)), na.value="transparent", limits = c(legend_limit_lower[ipara], legend_limit_upper[ipara]), trans = process_scale_option[ipara], oob = scales::squish) +
    geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
    geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
    # change the background to black and white
    coord_equal() +
    # theme_map() +
    ylim(lat_limits_robin[ , 2]) +
    # change the legend properties
    # theme(legend.position = 'none') +
    theme(legend.justification = c(0.1, 0.36), legend.position = c(0.1, 0.36), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
    # change the size of colorbar
    guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
    theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
    # add title
    labs(title = process_name[ipara], x = '', y = '') + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
    # modify the font size
    theme(axis.title = element_text(size = 20)) + 
    theme(panel.background = element_rect(fill = NA, colour = NA)) +
    # modify the margin
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
    theme(axis.text=element_text(size = 35))
  
  eval(parse(text = paste('p', ipara, ' = p', sep = '')))
  
  
  CurrentData = data.frame(cbind(global_lat_lon, bulk_process_std[ , ipara]))
  colnames(CurrentData) = c('Lon', 'Lat', 'Project')
  
  lon_lat_transfer = project(xy = as.matrix(CurrentData[ , c('Lon', 'Lat')]), proj = '+proj=robin') 
  CurrentData[ , c('Lon', 'Lat')] = lon_lat_transfer
  
  
  p =
    ggplot() +
    geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = Project), height = 60000, width = 60000, na.rm = TRUE) +
    scale_fill_gradientn(name = process_unit[ipara], colours = rev(viridis(15)), na.value="transparent", limits = c(legend_limit_lower_uncertain[ipara], legend_limit_upper_uncertain[ipara]), trans = process_scale_option[ipara], oob = scales::squish) +
    geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
    geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
    # change the background to black and white
    coord_equal() +
    # theme_map() +
    ylim(lat_limits_robin[ , 2]) +
    # change the legend properties
    # theme(legend.position = 'none') +
    theme(legend.justification = c(0.1, 0.36), legend.position = c(0.1, 0.36), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
    # change the size of colorbar
    guides(fill = guide_colorbar(direction = 'vertical', barwidth = 2.5, barheight = 14, title.position = 'top', title.hjust = 0, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
    theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
    # add title
    labs(title = paste('Uncertainty ', process_name[ipara], sep = ''), x = '', y = '') + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
    # modify the font size
    theme(axis.title = element_text(size = 20)) + 
    theme(panel.background = element_rect(fill = NA, colour = NA)) +
    # modify the margin
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
    theme(axis.text=element_text(size = 35))
  
  eval(parse(text = paste('p_std', ipara, ' = p', sep = '')))
  
}

box_soc = box_soc + 
  scale_fill_manual(name = '', labels = c('Tropical', 'Arid', 'Temperate', 'Boreal'), values = color_scheme) +
  theme(legend.justification = c(1, 0.1), legend.position = 'None', legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.7, 'inch'))

box_4 = box_4 + 
  scale_fill_manual(name = '', labels = c('Tropical', 'Arid', 'Temperate', 'Boreal'), values = color_scheme) +
  theme(legend.justification = c(1, 0.1), legend.position = c(-1.5, -0.22), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  theme(legend.text = element_text(size = 50), legend.title = element_text(size = 50))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch'))

jpeg(paste('./Ensemble/Map_bulk_process_cross_valid_mean_proda_full_CESM2.jpeg', sep = ''), width = 35, height = 25, units = 'in', res = 300)
plot_grid(p_soc, box_soc, p_tau, box_tau, 
          p1, box_1, p6, box_6, 
          p2, box_2, p3, box_3, 
          p5, box_5, p4, box_4,
          NULL, NULL, NULL, NULL,
          nrow = 5, ncol = 4,
          rel_widths = c(3, 1, 3, 1),
          rel_heights = c(1, 1, 1, 1, 0.1),
          labels = c('a', ' ', 'b', ' ',
                     'c', ' ', 'd', ' ',
                     'e', ' ', 'f', ' ',
                     'g', ' ', 'h', ' ', 
                     ' ', ' ', ' ', ' '),
          label_size = 60,
          label_x = 0.1, label_y = 1,
          label_fontfamily = 'Arial',
          label_fontface = 'plain'
          )
dev.off()


box_soc = box_soc + 
  scale_fill_manual(name = '', labels = c('Tropical', 'Arid', 'Temperate', 'Boreal'), values = color_scheme) +
  theme(legend.justification = c(1, 0.1), legend.position = c(-0.3, -0.12), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.7, 'inch'))

box_4 = box_4 + 
  scale_fill_manual(name = '', labels = c('Tropical', 'Arid', 'Temperate', 'Boreal'), values = color_scheme) +
  theme(legend.justification = 'None' , legend.position = 'None', legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  theme(legend.text = element_text(size = 50), legend.title = element_text(size = 50))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch'))

jpeg(paste('./Ensemble/Map_bulk_process_cross_valid_mean_no_soc_proda_full_CESM2.jpeg', sep = ''), width = 35, height = 20, units = 'in', res = 300)
plot_grid(p1, box_1, p6, box_6, 
          p2, box_2, p3, box_3, 
          p5, box_5, p4, box_4, 
          nrow = 3, ncol = 4,
          rel_widths = c(3, 1, 3, 1),
          rel_heights = c(1, 1, 1, 1)
)
dev.off()


jpeg(paste('./Ensemble/Map_bulk_process_cross_valid_std_proda_CESM2.jpeg', sep = ''), width = 27, height = 25, units = 'in', res = 300)
plot_grid(p_std_soc, p_std_tau, 
          p_std1, p_std6,
          p_std2, p_std3,
          p_std5, p_std4, nrow = 4)
dev.off()


jpeg(paste('./Ensemble/Map_bulk_process_cross_valid_mean_proda_CESM2.jpeg', sep = ''), width = 45, height = 15, units = 'in', res = 300)
plot_grid(p1, p6, p2, p3, p5, p4, nrow = 2)
dev.off()


################################## land cover box plot
for (ipara in 1:length(process_name)){
  CurrentData = data.frame(cbind(global_lat_lon, bulk_process_mean[ , ipara]))
  colnames(CurrentData) = c('Lon', 'Lat', 'Project')
  CurrentData$landcover = var_landcover
  
  color_scheme = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00')
  p =
    ggplot(data = CurrentData, aes(x = landcover, y = Project, group = landcover, fill = landcover)) +
    geom_boxplot(size = 1, outlier.colour = 'blue3', outlier.shape = 16, outlier.size = 0.3) +
    stat_summary(fun = 'mean', geom = 'point', shape = 18, size = 7, color = 'white') +
    scale_y_continuous(limits = c(legend_limit_lower[ipara], legend_limit_upper[ipara]), trans = process_scale_option[ipara], position = 'left') +
    scale_x_discrete(label = c('Agriculture', 'Boreal Forest', 'Other Forest', 'Grassland & Shrubland', 'Wetland', 'Bare Land')) +
    scale_fill_manual(name = '', labels = '', values = color_scheme) +
    # change the background to black and white
    theme_classic() +
    # change the legend properties
    theme(legend.position = 'none') +
    # add title
    labs(title = process_axis_label[ipara], x = '', y = '') +
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
    # modify the font size
    theme(axis.title = element_text(size = 30)) + 
    # modify the margin
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
    theme(axis.text=element_text(size = 30)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch'))
  
  
  eval(parse(text = paste('box_landcover', ipara, ' = p', sep = '')))
  
}

jpeg(paste('./Ensemble/boxplot_landcover_bulk_process_proda_CESM2.jpeg', sep = ''), width = 25, height = 15, units = 'in', res = 300)
plot_grid(box_landcover1, box_landcover2, box_landcover3, box_landcover4, box_landcover5, box_landcover6, nrow = 2)
dev.off()



########################### land cover and climate zones distribution map
color_scheme = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00')

# soc map
CurrentData = data.frame(cbind(global_lat_lon))
colnames(CurrentData) = c('Lon', 'Lat')

CurrentData$var_landcover = var_landcover

lon_lat_transfer = project(xy = as.matrix(CurrentData[ , c('Lon', 'Lat')]), proj = '+proj=robin') 
CurrentData[ , c('Lon', 'Lat')] = lon_lat_transfer

map_landcover = 
  ggplot() +
  geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = var_landcover), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_manual(name = '', labels = c('Agriculture', 'Boreal Forest', 'Other Forest', 'Grassland & Shrubland', 'Wetland', 'Bare Land'), values = color_scheme) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
  # change the background to black and white
  coord_equal() +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0.1, 0.1), legend.position = c(0.1, 0.1), legend.background = element_rect(fill = NA), legend.direction = 'vertical') +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # change the size of colorbar
  theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
  # add title
  labs(title = 'Land Cover Types', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35))


color_scheme = c('#D55E00', '#E69F00', '#009E73', '#0072B2')

CurrentData$var_climate = var_climate
# climate zones map
map_climate_zones =
  ggplot() +
  geom_tile(data = CurrentData, aes(x = Lon, y = Lat, fill = var_climate), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_manual(name = '', labels = c('Tropical', 'Arid', 'Temperate', 'Boreal'), values = color_scheme) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 1) +
  # change the background to black and white
  coord_equal() +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0.1, 0.1), legend.position = c(0.1, 0.1), legend.background = element_rect(fill = NA), legend.direction = 'vertical') +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # change the size of colorbar
  theme(legend.text = element_text(size = 20, ), legend.title = element_text(size = 30)) +
  # add title
  labs(title = 'Climate Zones', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35))


jpeg(paste('./Ensemble/climate_zones_and_landcover_distribution.jpeg', sep = ''), width = 15, height = 14, units = 'in', res = 300)
plot_grid(map_climate_zones, map_landcover, nrow = 2)
dev.off()

jpeg(paste('./Ensemble/climate_zones_distribution.jpeg', sep = ''), width = 15, height = 7, units = 'in', res = 300)
map_climate_zones + 
  labs(title = '', x = '', y = '') + 
  theme(legend.justification = c(0.15, 0.1), legend.position = c(0.15, 0.1), legend.background = element_rect(fill = NA), legend.direction = 'vertical') 

dev.off()



