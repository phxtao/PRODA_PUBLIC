library(semPlot)
library(lavaan)

library(R.matlab)
library(ggplot2)
library(cowplot)
library(lme4)
library(lmerTest)
library(viridis)

library(pracma)

rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

##################################################
## load emperical data
##################################################
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

var_nn_list = c('Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',
                'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',
                'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm',
                'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm',
                'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm',
                'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', 
                'CEC_0cm', 'CEC_30cm', 'CEC_100cm',
                'Garde_Acid',
                'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm',
                'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm',
                'USDA_Suborder',
                'WRB_Subgroup',
                'Annual Mean Temperature', 'Mean Diurnal Range', 'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 'Koppen_Climate_2018',
                'ESA_Land_Cover', 
                'cesm2_npp', 'cesm2_npp_std', 
                'cesm2_vegc',
                'Lon', 'Lat',
                'Abs_Depth_to_Bedrock',
                'Occurrence_R_Horizon',
                'nbedrock',
                'Elevation')



soil_var_texture_list = c('Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',
                          'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',
                          'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm',
                          'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm',
                          'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm',
                          'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm')


soil_var_chemical_list = c('CEC_0cm', 'CEC_30cm', 'CEC_100cm',
                           'Garde_Acid',
                           'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm',
                           'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm',
                           'USDA_Suborder',
                           'WRB_Subgroup')

climate_var_list = c('Annual Mean Temperature', 'Mean Diurnal Range', 'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 'Koppen_Climate_2018')

vegetation_var_list =  c('ESA_Land_Cover', 
                         'cesm2_npp', 'cesm2_npp_std', 
                         'cesm2_vegc')

geography_var_list =  c('Lon', 'Lat',
                        'Abs_Depth_to_Bedrock',
                        'Occurrence_R_Horizon',
                        'nbedrock',
                        'Elevation')


env_info = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/CUE_Synthesis/cue_meta_envinfo.mat')
env_info = env_info$EnvInfo[ , 4:80]
colnames(env_info) = grid_var_names
cue_meta_info = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/CUE_Synthesis/cue_meta.mat')
cue_meta_info = cue_meta_info$meta.site.info
colnames(cue_meta_info) = c('site_id', 'lon', 'lat', 'mat', 'incubation_temp', 'map', 'depth', 'cue', 'mic', 'soc', 'pH', 'cue_method', 'source_id')

cue_meta_info[is.na(cue_meta_info[, 'mat']) == 1, 'mat'] = env_info[is.na(cue_meta_info[, 'mat']) == 1, 'Annual Mean Temperature']
cue_meta_info[is.na(cue_meta_info[, 'pH']) == 1, 'pH'] = env_info[is.na(cue_meta_info[, 'pH']) == 1, 'pH_Water_0cm']/10
cue_meta_info[cue_meta_info[ , 'cue_method'] < 3, 'cue_method'] = 1

valid_loc = which(cue_meta_info[ , 'source_id'] < 1000)

####################################################################
# Geographical Location
####################################################################
world_coastline = rgdal::readOGR(dsn='/Users/phoenix/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp',layer = 'ne_110m_land')
world_coastline <- fortify(world_coastline)
Map.Using = world_coastline

current_data_middle = cue_meta_info[valid_loc, ]
site_list = unique(current_data_middle[ , 1])
current_data = data.frame(array(NA, dim = c(length(site_list), 3)))
colnames(current_data) = c('lon', 'lat', 'num')

isite = 1
for (isite in 1:length(site_list)) {
  site_loc = which(current_data_middle[ , 1] == site_list[isite])
  current_data[isite, 'lon'] = current_data_middle[site_loc[1], 'lon']
  current_data[isite, 'lat'] = current_data_middle[site_loc[1], 'lat']
  current_data[isite, 'num'] = length(site_loc)
}


jpeg(paste('./Ensemble/meta_site_loc.jpeg', sep = ''), width = 10, height = 5, units = 'in', res = 300)

ggplot() +
  geom_point(data = current_data, aes(x = lon, y = lat, size = num), color = 'red3', alpha = 0.7, na.rm = TRUE) +
  scale_size_continuous(name = 'Record Number', breaks = c(1, 4, 8)) + 
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
  # change the background to black and white
  theme_bw() +
  ylim(-56, 80) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0.5, 0), legend.position = c(0.5, 0), legend.background = element_rect(fill = NA), legend.direction = 'horizontal') +
  # change the size of colorbar
  guides(fill = guide_colorbar(barwidth = 2.5, barheight = 14), reverse = FALSE) + 
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  # add title
  labs(title = '', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 35))

dev.off()

####################################################################
# Correlation plot CUE-SOC
####################################################################
current_data = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'incubation_temp', 'cue', 'soc', 'depth', 'pH', 'cue_method')])
colnames(current_data) = c('site', 'source', 'mat', 'incubation_temp', 'cue', 'soc', 'depth', 'pH', 'cue_method')
current_data$site = as.factor(current_data$site)
current_data$source = as.factor(current_data$source)
current_data$cue_method = as.factor(current_data$cue_method)

## soc ~ cue
mix_model = lmer(soc ~ cue + depth + (cue | source), data = current_data)
summary(mix_model)
anova(mix_model)

r_squared = summary(lm(log(current_data$soc) ~ predict(mix_model)))
r_squared = r_squared$r.squared
r_squared


coef_mix_model_c1314 = coef(mix_model)$cue_method[1, ]
coef_mix_model_o18 = coef(mix_model)$cue_method[2, ]
coef_mix_model_all = fixef(mix_model)

cue_seq = seq(from = min(current_data$cue), by = 0.01, to = max(current_data$cue))
pred_soc_c1314 = exp(coef_mix_model_c1314$`(Intercept)` + coef_mix_model_c1314$cue*cue_seq + coef_mix_model_c1314$depth*15)
pred_soc_o18 = exp(coef_mix_model_o18$`(Intercept)` + coef_mix_model_o18$cue*cue_seq + coef_mix_model_o18$depth*15)
pred_soc_all = exp(coef_mix_model_all['(Intercept)'] + coef_mix_model_all['cue']*cue_seq + coef_mix_model_all['depth']*15)

corr_test_all = cor.test(current_data$cue, current_data$soc)
corr_test_all

corr_test_o18 = cor.test(current_data$cue[current_data$cue_method == 3], current_data$soc[current_data$cue_method == 3])
corr_test_o18

corr_test_c1314 = cor.test(current_data$cue[current_data$cue_method == 1], current_data$soc[current_data$cue_method == 1])
corr_test_c1314


text_data_all = data.frame('x_axis' = c(0.05),
                             'y_axis' = c(50),
                             'equation' = paste('r = ', round(corr_test_all$estimate, 2), '\nP < 0.001', sep = ''))

text_data_o18 = data.frame('x_axis' = c(0.05),
                       'y_axis' = c(90),
                       'equation' = paste('r = ', round(corr_test_o18$estimate, 2), '\nP < 0.001', sep = ''))

text_data_c1314 = data.frame('x_axis' = c(0.05),
                           'y_axis' = c(50),
                           'equation' = paste('r = ', round(corr_test_c1314$estimate, 2), '\nP < 0.05', sep = ''))


color_scheme = c('#005AB5', '#DC3220')

jpeg(paste('./Ensemble/plot_cue_soc_log.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)
ggplot() +
  geom_point(data = current_data, aes(x = cue, y = soc, size = depth), color = 'black', alpha = 0.7, na.rm = TRUE) +
  # geom_point(data = current_data, aes(x = cue, y = soc, group = cue_method, color = cue_method, size = depth), alpha = 0.7, na.rm = TRUE) +
  # geom_smooth(data = current_data, aes(x = cue, y = soc, group = cue_method, color = cue_method, fill = cue_method), alpha = 0.2, method = 'lm', size = 2) +
  geom_smooth(data = current_data, aes(x = cue, y = soc), alpha = 0.2, method = 'lm', color = 'black', fill = 'grey45', size = 2) +
  scale_x_continuous(trans = 'identity') + 
  scale_y_continuous(limits = c(0.3, 600), trans = 'log10') + 
  scale_color_manual(name = 'CUE Methods', labels = c(expression(paste('Isotope C ('^'13', 'C/', ''^'14', 'C)', sep = '')), expression(paste(''^'18', 'O', sep = ''))),  values = color_scheme) +
  scale_fill_manual(name = 'CUE Methods',  labels = c(expression(paste('Isotope C ('^'13', 'C/', ''^'14', 'C)', sep = '')), expression(paste(''^'18', 'O', sep = ''))), values = color_scheme) +
  scale_size_continuous(name = 'Depth (cm)', range = c(3, 9), breaks = c(5, 15, 30)) + 
  # change the background to black and white
  geom_text(data = text_data_all, aes(x = 0.05, y = 400, label = equation), color = 'black', hjust = 0, size = 10, show.legend = FALSE) + 
  # geom_text(data = text_data_c1314, aes(x = 0.6, y = 400, label = equation), color = color_scheme[1], hjust = 0, size = 10, show.legend = FALSE) + 
  # geom_text(data = text_data_o18, aes(x = 0.325, y = 400, label = equation), color = color_scheme[2], hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', x = 'CUE', y = expression(paste('SOC (g C kg'^'-1', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(legend.direction = 'horizontal') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


dev.off()



####################################################################
# Correlation plot CUE-MAT
####################################################################
current_data = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'incubation_temp', 'cue', 'soc', 'depth', 'pH', 'cue_method')])
colnames(current_data) = c('site', 'source', 'mat', 'incubation_temp', 'cue', 'soc', 'depth', 'pH', 'cue_method')
current_data$site = as.factor(current_data$site)
current_data$source = as.factor(current_data$source)
current_data$cue_method = as.factor(current_data$cue_method)

## cue ~ mat
mix_model = lmer(cue ~ mat + (1 | cue_method), data = current_data)
summary(mix_model)



corr_test_all = cor.test(current_data$cue, current_data$mat)
corr_test_all
text_data_all = data.frame('x_axis' = c(0.05),
                           'y_axis' = c(90),
                           'equation' = paste('r = ', round(corr_test_all$estimate, 2), '\nP < 0.001', sep = ''))

corr_test_o18 = cor.test(current_data$cue[current_data$cue_method == 3], current_data$mat[current_data$cue_method == 3])
corr_test_o18
text_data_o18 = data.frame('x_axis' = c(0.05),
                           'y_axis' = c(90),
                           'equation' = paste('r = ', round(corr_test_o18$estimate, 2), '\nP < 0.01', sep = ''))

corr_test_o1314 = cor.test(current_data$cue[current_data$cue_method == 1], current_data$mat[current_data$cue_method == 1])
corr_test_o1314
text_data_c1314 = data.frame('x_axis' = c(0.05),
                             'y_axis' = c(50),
                             'equation' = paste('r = ', round(corr_test_o1314$estimate, 2), '\nP < 0.05', sep = ''))



color_scheme = c('#005AB5', '#DC3220')

jpeg(paste('./Ensemble/plot_mat_cue.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)


ggplot() +
  geom_point(data = current_data, aes(x = mat, y = cue, size = depth), color = 'black', alpha = 0.7, na.rm = TRUE) +
  # geom_point(data = current_data, aes(x = mat, y = cue, group = cue_method, color = as.factor(cue_method), size = depth), alpha = 0.7, na.rm = TRUE) +
  # geom_smooth(data = current_data, aes(x = mat, y = cue, group = cue_method, color = cue_method, fill = cue_method), alpha = 0.2, method = 'lm', size = 2) +
  geom_smooth(data = current_data, aes(x = mat, y = cue), method = 'lm', color = 'black', fill = 'grey45', alpha = 0.2, size = 2) +
  scale_x_continuous(trans = 'identity') + 
  scale_y_continuous(limits = c(0., 0.9), trans = 'identity') + 
  scale_color_manual(name = 'CUE Methods', labels = c(expression(paste('Isotope C ('^'13', 'C/', ''^'14', 'C)', sep = '')), expression(paste(''^'18', 'O', sep = ''))),  values = color_scheme) +
  scale_fill_manual(name = 'CUE Methods',  labels = c(expression(paste('Isotope C ('^'13', 'C/', ''^'14', 'C)', sep = '')), expression(paste(''^'18', 'O', sep = ''))), values = color_scheme) +
  scale_size_continuous(name = 'Depth (cm)', range = c(3, 9), breaks = c(5, 15, 30)) + 
  # change the background to black and white
  geom_text(data = text_data_all, aes(x = -3, y = 0.85, label = equation), color = 'black', hjust = 0, size = 10, show.legend = FALSE) + 
  # geom_text(data = text_data_c1314, aes(x = 15, y = 0.85, label = equation), color = color_scheme[1], hjust = 0, size = 10, show.legend = FALSE) + 
  # geom_text(data = text_data_o18, aes(x = 6, y = 0.85, label = equation), color = color_scheme[2], hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', x = 'deg C', y = 'CUE') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(legend.direction = 'horizontal') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

dev.off()
