## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(jcolors)
library(viridis)
library(cowplot)
library(lme4)
library(lmerTest)
library(scales)

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
# meta-analysis
#################################################################################
## load emperical data
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

# cue from proda 
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

global_lat_lon = readMat(paste( data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_NPP_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(7), '.mat', sep = ''))

global_lat_lon = global_lat_lon$var.data.middle[[1]][[1]]
global_lat_lon = global_lat_lon[ , 1:2]

cue_proda_grid = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_NPP_', model_name, '_', nn_exp_name, '_cross_valid_0_', as.character(7), '.mat', sep = ''))

cue_proda_grid = Re(cue_proda_grid$var.data.middle[[4]][[1]])[ , 11]

cue_proda_meta = array(NA, dim = c(nrow(cue_meta_info), 2))

isite  = 1
for (isite in 1:nrow(cue_meta_info)) {
  cue_proda_meta[isite, 1] = cue_meta_info[isite, 'cue']
  length_middle = (global_lat_lon[ , 1] - cue_meta_info[isite, 'lon'])**2 +  (global_lat_lon[ , 2] - cue_meta_info[isite, 'lat'])**2
  grid_loc = which(length_middle == min(length_middle, na.rm = TRUE))
  
  cue_proda_meta[isite, 2] = cue_proda_grid[grid_loc[1]]
}

valid_loc = which(cue_meta_info[ , 'source_id'] < 1000)

cor.test(cue_proda_meta[valid_loc, 1], cue_proda_meta[valid_loc, 2])
# plot(cue_proda_meta[valid_loc, 1], cue_proda_meta[valid_loc, 2])

#---------------------------------------------------
# CUE-SOC relationship
#---------------------------------------------------

#---------------mix model microbial vs non-microbial biomass---------------#
current_data_meta = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'cue', 'mic', 'soc', 'depth', 'cue_method')])
current_data_meta = current_data_meta[is.na(apply(current_data_meta, 1, sum, na.rm = FALSE)) == 0, ]

mix_model = lmer(mic/1000 ~ cue + depth + mat + (1 | source_id), data = current_data_meta)
mix_model_summary_meta = summary(mix_model)
mix_model_summary_meta
r_squared_mixed = summary(lm(current_data_meta$mic/1000 ~ predict(mix_model)))
r_squared_mixed


mix_model = lmer((soc - mic/1000) ~ cue + depth + mat + (1 | source_id), data = current_data_meta)
mix_model_summary_meta = summary(mix_model)
mix_model_summary_meta
r_squared_mixed = summary(lm((current_data_meta$soc - current_data_meta$mic/1000) ~ predict(mix_model)))
r_squared_mixed

#---------------mix model---------------#
current_data_meta = data.frame(cue_meta_info[valid_loc, c('site_id', 'source_id', 'mat', 'cue', 'soc', 'depth', 'cue_method')])
current_data_meta = current_data_meta[is.na(apply(current_data_meta, 1, sum, na.rm = FALSE)) == 0, ]

mix_model = lmer(log10(soc) ~ cue + depth + mat + (1 | source_id), data = current_data_meta)
mix_model_summary_meta = summary(mix_model)
mix_model_summary_meta

fit_function_meta = function(x) {10**(mix_model_summary_meta$coefficients[1, 1] + x*(mix_model_summary_meta$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data_meta$soc) ~ predict(mix_model)))
r_squared_mixed

text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(1000), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_meta$coefficients[1, 1], 2), ' + ', round(mix_model_summary_meta$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_meta$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_meta$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))

color_scheme = c('#005AB5', '#DC3220')

p_meta_cue_soc =
  ggplot() +
  geom_point(data = current_data_meta, aes(x = cue, y = soc, size = depth), color = 'black', alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_meta, size = 2, color = 'black') + 
  scale_x_continuous(trans = 'identity') + 
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  scale_size_continuous(name = 'Depth (cm)', range = c(3, 9), breaks = c(5, 15, 30)) + 
  # change the background to black and white
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() +
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon  (g C kg'^'-1', ')', sep = ''))) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  theme(legend.direction = 'horizontal') + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


#################################################################################
# mcrobial model
#################################################################################
version_num = 7
model_name = paste('cesm2_clm5_mic_vr_v', version_num, sep = '')
model_name_beta = paste(strsplit(model_name, split = '_')[[1]], collapse = '.')

data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'

var_names = c('ProfileNum', 'ProfileID', 'LayerNum', 'Lon', 'Lat', 'Date',
              'Rmean', 'Rmax', 'Rmin',
              'ESA_Land_Cover',
              'ET', 
              'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin',
              'Veg_Cover', 
              'BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', 
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
              'nbedrock',
              'R_Squared')

para_name = c('bio', 'cryo', 'q10', 'efolding', 'w_scaling', 
              'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', 
              'mm_const_assim', 'mm_const_decom', 
              'fcwdl2', 'fl1s1', 'fl2s1', 'fl3s4', 
              'mic_cue', 'pdeath2soc', 
              'beta', 
              'allo_slope_mic')


#----------------------------------------------------
# microbial profile info
#----------------------------------------------------
output_name_list = c('candidate_para_value', 'candidate_cost_value', 'candidate_steady_state', 'candidate_r2_hist', 'candidate_para_std_hist', 
                     'wosis_layer_obs', 'soc_mod_opt', 'soc_stock_opt', 'soc_layer_opt', 'soc_std_opt', 'para_mean', 'r2_opt')


npara = length(para_name)

profile_env_info = readMat(paste(data_dir_input, '/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat', sep = ''))
profile_env_info = profile_env_info$EnvInfo
colnames(profile_env_info) = var_names

# load climate info
global_climate = profile_env_info[ , 'Koppen_Climate_2018']
var_climate = array(NA, dim = c(length(global_climate), 1))
# koppen climatte class 2
var_climate = global_climate
var_climate[which(global_climate == 1)] = 101 # Af
var_climate[which(global_climate == 2)] = 102 # Am
var_climate[which(global_climate == 3)] = 103 # Aw
var_climate[which(global_climate >= 4 & global_climate <= 5)] = 104 # BwX
var_climate[which(global_climate >= 6 & global_climate <= 7)] = 105 # BsW
var_climate[which(global_climate >= 8 & global_climate <= 10)] = 106 # CsX
var_climate[which(global_climate >= 11 & global_climate <= 13)] = 107 # CwX
var_climate[which(global_climate >= 14 & global_climate <= 16)] = 108 # CfX
var_climate[which(global_climate >= 17 & global_climate <= 20)] = 109 # DsX
var_climate[which(global_climate >= 21 & global_climate <= 24)] = 110 # DwX
var_climate[which(global_climate >= 25 & global_climate <= 28)] = 111 # DfX
var_climate[which(global_climate >= 29 & global_climate <= 30)] = 112 # E


#----------------------------------------------------
# microbial load da summary
#----------------------------------------------------
col_list_relation_summary = c('profile_num',
                              'cue', 'mm_ratio_assim', 'mm_ratio_decom', 'tau4s2_enz', 'allo_slope_mic',
                              'soc_0_30cm', 'soc_30_100cm', 'soc_100_200cm', 'mic_0_30cm', 'mic_stock', 'enz_stock', 'doc_stock', 'poc_stock', 'soc_total', 'r2')

relation_summary = readMat(paste(data_dir_output, '/da_summary_', model_name, '/', model_name, '_da_summary_mic.mat', sep = ''))
relation_summary = relation_summary$da.summary.mic

colnames(relation_summary) = col_list_relation_summary

relation_summary[ , 'soc_0_30cm'] = relation_summary[ , 'soc_0_30cm']/0.3/profile_env_info[ , 'Bulk_Density_0cm'] # unit gc/kg
relation_summary[ , 'soc_30_100cm'] = relation_summary[ , 'soc_30_100cm']/0.7/profile_env_info[ , 'Bulk_Density_30cm'] # unit gc/kg
relation_summary[ , 'soc_100_200cm'] = relation_summary[ , 'soc_100_200cm']/1/profile_env_info[ , 'Bulk_Density_100cm'] # unit gc/kg
relation_summary[ , 'mic_0_30cm'] = relation_summary[ , 'mic_0_30cm']/0.3/profile_env_info[ , 'Bulk_Density_0cm'] # unit gc/kg


bulk_process_summary = readMat(paste(data_dir_output, '/da_summary_', model_name, '/', model_name, '_bulk_process_opt.mat', sep = ''))
bulk_process_summary = bulk_process_summary$bulk.process.opt
#-----------------------------------------------
# representative points
#-----------------------------------------------
sample_profile_matrix_num = readMat(paste(data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_representative_profiles.mat', sep = ''))
sample_profile_matrix_num = sample_profile_matrix_num$sample.profile.id

sample_profile_matrix_ecoregion_class1 = matrix(profile_env_info[as.vector(sample_profile_matrix_num), which(var_names == 'Koppen_Climate_2018')], nrow = 100, ncol = 50)
sample_profile_matrix_ecoregion_class2 = sample_profile_matrix_ecoregion_class1
sample_profile_matrix_ecoregion_class3 = sample_profile_matrix_ecoregion_class1

# koppen climatte class
koppen_climate_class = data.frame(array(NA, dim = c(30, 2)))
colnames(koppen_climate_class) = c('class1', 'class2')
# koppen climatte class 1
koppen_climate_class$class1 = seq(1, 30)
# koppen climatte class 2
koppen_climate_class$class2[which(koppen_climate_class$class1 == 1)] = 101 # Af
koppen_climate_class$class2[which(koppen_climate_class$class1 == 2)] = 102 # Am
koppen_climate_class$class2[which(koppen_climate_class$class1 == 3)] = 103 # Aw
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 4 & koppen_climate_class$class1 <= 5)] = 104 # BwX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 6 & koppen_climate_class$class1 <= 7)] = 105 # BsW
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 8 & koppen_climate_class$class1 <= 10)] = 106 # CsX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 11 & koppen_climate_class$class1 <= 13)] = 107 # CwX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 14 & koppen_climate_class$class1 <= 16)] = 108 # CfX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 17 & koppen_climate_class$class1 <= 20)] = 109 # DsX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 21 & koppen_climate_class$class1 <= 24)] = 110 # DwX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 25 & koppen_climate_class$class1 <= 28)] = 111 # DfX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 29 & koppen_climate_class$class1 <= 30)] = 112 # E

for (iclimate in 1:30) { 
  sample_profile_matrix_ecoregion_class2[which(sample_profile_matrix_ecoregion_class1 == koppen_climate_class$class1[iclimate])] = koppen_climate_class$class2[iclimate]
}


#-----------------------------------------------
# sample reorder
#-----------------------------------------------
valid_profile_summary = c()
for (isample in 1:ncol(sample_profile_matrix_num)) {
  valid_profile_summary = rbind(valid_profile_summary, cbind(sample_profile_matrix_num[ , isample], sample_profile_matrix_ecoregion_class2[ , isample],  isample))
}

colnames(valid_profile_summary) = c('profile_num', 'climate', 'sample')
sample_order = valid_profile_summary[ , 'sample']

is_valid_vector = relation_summary[valid_profile_summary[ , 1], 'cue']
is_valid_vector[is.na(is_valid_vector) == 0] = 1
is_valid_vector[is.na(is_valid_vector) == 1] = 0

isample = 1
for (isample in 1:nrow(valid_profile_summary)) {
  sample_climate = valid_profile_summary[isample, 'climate']
  middle_loc = which(valid_profile_summary[ , 'climate'] == valid_profile_summary[isample, 'climate'] & is_valid_vector == 1)
  middle_loc = middle_loc[middle_loc > isample]
  
  if (length(middle_loc) >=1) {
    middle_loc = middle_loc[1]
    valid_profile_summary[c(isample, middle_loc), ] = valid_profile_summary[c(middle_loc, isample), ]
    is_valid_vector[c(isample, middle_loc)] = is_valid_vector[c(middle_loc, isample)]
  }
}

valid_profile_summary[ , 'sample'] = sample_order

valid_profile_loc = valid_profile_summary[1:1000, ]
#----------------------------------------------------
# cue soc relation 0 - 30 cm
#----------------------------------------------------

#---------------mix model microbial vs non-microbial biomass---------------#
current_data_mic = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm', 'mic_0_30cm', 'cue')], valid_profile_loc))
# current_data_mic = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm')], bulk_process_summary[valid_profile_loc[ , 1], 1], valid_profile_loc))
colnames(current_data_mic) = c('soc', 'mic', 'cue', 'profile_num', 'climate', 'sample')

current_data_mic = current_data_mic[is.na(apply(current_data_mic, 1, sum, na.rm = FALSE)) == 0, ]

mix_model = lmer(mic ~ cue + (1 | climate), data = current_data_mic)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic

r_squared_mixed = summary(lm(current_data_mic$mic ~ predict(mix_model)))
r_squared_mixed 


mix_model = lmer((soc-mic) ~ cue + (1 | climate), data = current_data_mic)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic

r_squared_mixed = summary(lm((current_data_mic$soc-current_data_mic$mic) ~ predict(mix_model)))
r_squared_mixed 


#---------------mix model---------------#
current_data_mic = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm', 'cue')], valid_profile_loc))
# current_data_mic = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm')], bulk_process_summary[valid_profile_loc[ , 1], 1], valid_profile_loc))
colnames(current_data_mic) = c('soc', 'cue', 'profile_num', 'climate', 'sample')

current_data_mic = current_data_mic[is.na(apply(current_data_mic, 1, sum, na.rm = FALSE)) == 0, ]

mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data_mic)
mix_model_summary_mic = summary(mix_model)
mix_model_summary_mic
fit_function_mic = function(x) {10**(mix_model_summary_mic$coefficients[1, 1] + x*(mix_model_summary_mic$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data_mic$soc) ~ predict(mix_model)))
r_squared_mixed 

fix_effect_model = summary(mix_model_summary_mic)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data_mic$cue
r_squared_fixed = summary(lm(log10(current_data_mic$soc) ~ fix_effect_pred))
r_squared_fixed


text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(1), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_mic$coefficients[1, 1], 2), ' + ', round(mix_model_summary_mic$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_mic$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_mic$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))

color_scheme = c('#b30000', '#e34a33', '#fc8d59', '#fdbb84', '#fdd49e', '#fef0d9', '#f0f9e8', '#ccebc5', '#a8ddb5', '#7bccc4', '#43a2ca', '#0868ac')

p_mic_cue_soc =
  ggplot(data = current_data_mic) + 
  geom_point(data = current_data_mic, aes(x = cue, y = soc, color = climate), size = 6, alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_mic, size = 2, color = 'black') + 
  # geom_smooth(data = current_data_mic, aes(x = cue, y = soc_0_30cm), method = 'lm', color = 'black', fill = 'black', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  # stat_bin_hex(aes(x = cue, y = soc_0_30cm), bins = 100) +
  scale_color_stepsn(name = 'Climate', colors = color_scheme, trans = 'identity', limits = c(100, 112), breaks = c(100:112), labels = c('',  'Tropics', '', '', '', '', '', '', '', '', '', '', 'Polar')) +
  scale_x_continuous(limits = c(0, max(current_data_mic$cue, na.rm = TRUE)), trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(color = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.2, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 



#################################################################################
# PRODA results
#################################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

#-----------------------------------------------
# Load Projected SOC site by site
#-----------------------------------------------

para_names =  c('diffus', 'cryo', 
                'q10', 'efolding', 
                'taucwd', 'taul1', 'taul2', 'tau4s1', 'tau4s2', 'tau4s3', 
                'fs1l1', 'fs1l2', 'fs2l3', 'fs2s1', 'fs3s1', 'fs1s2', 'fs3s2', 'fs1s3', 'fl2cwd', 
                'w-scaling', 'beta')


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

ss_para_result = readMat(paste(data_dir_output, 'mcmc_summary_', model_name, '/', model_name, '_para_mean.mat', sep = ''))
ss_para_result = ss_para_result$para.mean
colnames(ss_para_result) = para_names

env_info_ss = readMat(paste(data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat', sep = ''))
env_info_ss = env_info_ss$EnvInfo[ , 4:80]
colnames(env_info_ss) = grid_var_names

# site by site soc stock
soc_stock_ss = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_site_A_', model_name, '.mat', sep = ''))
soc_stock_ss = soc_stock_ss$Var.Decom.Grid

soc_stock_tau_mean = cbind(Re(soc_stock_ss[[5]][ , 6]), Re(soc_stock_ss[[9]][ , 6]))
bulk_process_mean = cbind(Re(soc_stock_ss[[11]][ , 6]), Re(soc_stock_ss[[15]][ , 6]), Re(soc_stock_ss[[12]][ , 6]), Re(soc_stock_ss[[14]][ , 6]), Re(soc_stock_ss[[13]][ , 6]), Re(soc_stock_ss[[16]][ , 6]))

soc_content_0_30cm = Re(soc_stock_ss[[2]][ , 6])
soc_tau_0_30cm = Re(soc_stock_ss[[6]][ , 6])


soc_content_0_30cm = soc_content_0_30cm/0.3/env_info_ss[ , 'Bulk_Density_0cm'] # unit gC/kg, unit of bulk density:# unit kg/m3 

# load climate info
global_climate = env_info_ss[ , 'Koppen_Climate_2018']
var_climate = array(NA, dim = c(length(global_climate), 1))
# koppen climatte class 2
var_climate = global_climate
var_climate[which(global_climate == 1)] = 101 # Af
var_climate[which(global_climate == 2)] = 102 # Am
var_climate[which(global_climate == 3)] = 103 # Aw
var_climate[which(global_climate >= 4 & global_climate <= 5)] = 104 # BwX
var_climate[which(global_climate >= 6 & global_climate <= 7)] = 105 # BsW
var_climate[which(global_climate >= 8 & global_climate <= 10)] = 106 # CsX
var_climate[which(global_climate >= 11 & global_climate <= 13)] = 107 # CwX
var_climate[which(global_climate >= 14 & global_climate <= 16)] = 108 # CfX
var_climate[which(global_climate >= 17 & global_climate <= 20)] = 109 # DsX
var_climate[which(global_climate >= 21 & global_climate <= 24)] = 110 # DwX
var_climate[which(global_climate >= 25 & global_climate <= 28)] = 111 # DfX
var_climate[which(global_climate >= 29 & global_climate <= 30)] = 112 # E

#-----------------------------------------------
# representative points and corresponding CUE-SOC relationship
#-----------------------------------------------
# sample_profile_matrix_num = readMat(paste(data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_representative_profiles.mat', sep = ''))
# sample_profile_matrix_num = sample_profile_matrix_num$sample.profile.id
# sample_profile_matrix_num = readMat(paste(data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_representative_profiles.mat', sep = ''))
# sample_profile_matrix_num = sample_profile_matrix_num$sample.profile.id
# sample_profile_matrix_num = sample_profile_matrix_num[ , 1:50]
# 
# valid_loc = c()
# for (isample in 1:ncol(sample_profile_matrix_num)) {
#   valid_loc = rbind(valid_loc, cbind(sample_profile_matrix_num[ , isample], var_climate[sample_profile_matrix_num[ , isample]],  isample))
# }

valid_loc = current_data_mic[ , c('profile_num', 'climate', 'sample')]

current_data_proda_sample = data.frame(cbind(valid_loc,
                                             bulk_process_mean[valid_loc[ , 1], 1], 
                                             soc_content_0_30cm[valid_loc[ , 1]], 
                                             env_info_ss[valid_loc[ , 1], 'Annual Mean Temperature'] 
))

colnames(current_data_proda_sample) = c('profile_num', 'climate', 'sample', 'cue', 'soc', 'mat')

current_data_proda_sample = current_data_proda_sample[is.na(apply(current_data_proda_sample, 1, sum, na.rm = FALSE)) == 0, ]

#------------------------- mixed model
mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data_proda_sample)
mix_model_summary_proda_sample = summary(mix_model)
mix_model_summary_proda_sample
fit_function_proda_sample = function(x) {10**(mix_model_summary_proda_sample$coefficients[1, 1] + x*(mix_model_summary_proda_sample$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data_proda_sample$soc) ~ predict(mix_model)))
r_squared_mixed 

fix_effect_model = summary(mix_model_summary_proda_sample)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data_proda_sample$cue
r_squared_fixed = summary(lm(log10(current_data_proda_sample$soc) ~ fix_effect_pred))
r_squared_fixed

text_data = data.frame('x_axis' = c(0.25), 
                       'y_axis' = c(1), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_proda_sample$coefficients[1, 1], 2), ' + ', round(mix_model_summary_proda_sample$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_proda_sample$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_proda_sample$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


#----------------------- Plot Figures SOC content and cue (0 - 30cm)
p_proda_cue_soc_sample =
  ggplot(data = current_data_proda_sample) + 
  geom_point(data = current_data_proda_sample, aes(x = cue, y = soc, color = climate), size = 6, alpha = 0.7, na.rm = TRUE) +
  scale_color_stepsn(name = 'Climate', colors = color_scheme, trans = 'identity', limits = c(100, 112), breaks = c(100:112), labels = c('',  'Tropics', '', '', '', '', '', '', '', '', '', '', 'Polar')) +
  # stat_bin_hex(aes(x = cue, y = soc), bins = 100) +
  # scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(limits = c(0.25, 0.65), trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  geom_function(fun = fit_function_proda_sample, size = 2, color = 'black') + 
  # geom_smooth(data = current_data_proda_sample, aes(x = cue, y = soc), method = 'lm', color = 'black', fill = 'blue3', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(color = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.2, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


#-----------------------------------------------
# all valid points and corresponding CUE-SOC relationship
#-----------------------------------------------
valid_loc = readMat(paste(data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info_', model_name, '_', time_domain, '_maxmin_scaled.mat', sep = ''))
valid_loc = valid_loc$profile.env.info[ , 1]

current_data_proda_total = data.frame(cbind(valid_loc, var_climate[valid_loc],
                                            bulk_process_mean[valid_loc, 1], 
                                            soc_content_0_30cm[valid_loc], 
                                            env_info_ss[valid_loc, 'Annual Mean Temperature'] 
))

colnames(current_data_proda_total) = c('profile_num', 'climate', 'cue', 'soc', 'mat')


current_data_proda_total = current_data_proda_total[is.na(apply(current_data_proda_total, 1, sum, na.rm = FALSE)) == 0, ]

#--------------------------- mixed model
mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data_proda_total)
mix_model_summary_proda_total = summary(mix_model)
mix_model_summary_proda_total
fit_function_proda_total = function(x) {10**(mix_model_summary_proda_total$coefficients[1, 1] + x*(mix_model_summary_proda_total$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data_proda_total$soc) ~ predict(mix_model)))
r_squared_mixed 

fix_effect_model = summary(mix_model_summary_proda_total)
fix_effect_pred = fix_effect_model$coefficients[1, 1] + 
  fix_effect_model$coefficients[2, 1]*current_data_proda_total$cue
r_squared_fixed = summary(lm(log10(current_data_proda_total$soc) ~ fix_effect_pred))
r_squared_fixed


text_data = data.frame('x_axis' = c(0.25), 
                       'y_axis' = c(1), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_proda_total$coefficients[1, 1], 2), ' + ', round(mix_model_summary_proda_total$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_proda_total$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_proda_total$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


#---------------------------- Plot Figures SOC content and cue (0 - 30cm)
p_proda_cue_soc_total =
  ggplot(data = current_data_proda_total) + 
  stat_bin_hex(aes(x = cue, y = soc), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
  scale_x_continuous(limits = c(0.25, 0.65), trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  geom_function(fun = fit_function_proda_total, size = 2, color = 'black') + 
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#################################################################################
# final figure
#################################################################################

jpeg(paste('./Ensemble/revision1_main_fig2.jpeg', sep = ''), width = 17, height = 17, units = 'in', res = 300)

p_meta_cue_soc = p_meta_cue_soc + 
  labs(x = '  ')
p_mic_cue_soc = p_mic_cue_soc + 
  labs(x = '  ', y = '  ')
p_proda_cue_soc_total = p_proda_cue_soc_total + 
  labs(y = '  ')

plot_grid(p_meta_cue_soc, p_mic_cue_soc,
          p_proda_cue_soc_sample, p_proda_cue_soc_total,
          labels = c('a', 'b', 'c', 'd'),
          label_size = 60,
          label_x = 0.02, label_y = 1.03,
          label_fontfamily = 'Arial',
          label_fontface = 'bold',
          nrow = 2)
dev.off()

# jpeg(paste('./Ensemble/revision1_main_fig2_standard.jpeg', sep = ''), width = 20, height = 20, units = 'in', res = 300)
# plot_grid(p_meta_cue_soc_standard, p_mic_cue_soc_standard,
#           p_proda_cue_soc_sample_standard, p_proda_cue_soc_total_standard,
#           labels = c('a', 'b', 'c', 'd'),
#           label_size = 60,
#           label_x = 0.02, label_y = 1,
#           label_fontfamily = 'Arial',
#           label_fontface = 'bold',
#           nrow = 2)
# dev.off()

#################################################################################
# variance
#################################################################################

current_data_variance = data.frame(rbind(cbind(current_data_meta$cue, 1), 
                                         cbind(current_data_proda_total$cue, 2), 
                                         cbind(current_data_mic$cue, 3)))
current_data_mean = data.frame(rbind(cbind(mean(current_data_meta$cue, na.rm = TRUE), 1), 
                                     cbind(mean(current_data_proda_total$cue, na.rm = TRUE), 2),
                                     cbind(mean(current_data_mic$cue, na.rm = TRUE), 3)))


rbind(c(mean(current_data_meta$cue, na.rm = TRUE), sd(current_data_meta$cue, na.rm = TRUE)), 
      c(mean(current_data_proda_total$cue, na.rm = TRUE), sd(current_data_proda_total$cue)),
      c(mean(current_data_mic$cue, na.rm = TRUE), sd(current_data_mic$cue, na.rm = TRUE)))


colnames(current_data_variance) = c('cue', 'source')
colnames(current_data_mean) = c('cue', 'source')

color_scheme = c('#1E88E5', '#FFC107', '#004D40')

jpeg(paste('./Ensemble/cue_variance.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

ggplot() +
  geom_histogram(data = current_data_variance, aes(x = cue, y = ..density.., color = as.factor(source), fill = as.factor(source)), position = 'identity', alpha = 0.3, na.rm = TRUE, size = 1, bins = 30) + 
  geom_density(data = current_data_variance, aes(x = cue, color = as.factor(source)), fill = NA, alpha = 0.5, na.rm = TRUE, size = 3) +
  geom_vline(data = current_data_mean, aes(xintercept = cue, color = as.factor(source)), linetype = 'longdash', size = 3) +
  scale_x_continuous(trans = 'identity') +
  scale_fill_manual(name = '', values = color_scheme, labels = c('Field measurements', 'CLM5', 'Microbial model')) +
  scale_color_manual(name = '', values = color_scheme, labels = c('Field measurements', 'CLM5', 'Microbial model')) +
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = 'Density') + 
  # change the legend properties
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

dev.off()


#################################################################################
# Intra-CUE
#################################################################################
current_data_intra_cue = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], 'cue'], bulk_process_mean[valid_profile_loc[ , 1], 1], var_climate[valid_profile_loc[ , 1]]))
colnames(current_data_intra_cue) = c('mic', 'cen', 'climate')

current_data_intra_cue = current_data_intra_cue[is.na(apply(current_data_intra_cue, 1, sum, na.rm = FALSE)) == 0, ]

cor.test(current_data_intra_cue$mic, current_data_intra_cue$cen)

