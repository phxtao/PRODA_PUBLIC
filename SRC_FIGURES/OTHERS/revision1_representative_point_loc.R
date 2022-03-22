## Packages
library(R.matlab)
library(maps)
library(ggplot2)
library(R.matlab)
library(cowplot)
library(jcolors)
library(viridis)
library(cowplot)
library(scales)

library(lme4)
library(lmerTest)


library(corrplot)
dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')


######################################################
# Data path
######################################################
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
                              'soc_0_30cm', 'soc_30_100cm', 'soc_100_200cm', 'mic_stock', 'enz_stock', 'doc_stock', 'poc_stock', 'soc_total', 'r2')

relation_summary = readMat(paste(data_dir_output, '/da_summary_', model_name, '/', model_name, '_da_summary_mic.mat', sep = ''))
relation_summary = relation_summary$da.summary.mic

colnames(relation_summary) = col_list_relation_summary

relation_summary[ , 'soc_0_30cm'] = relation_summary[ , 'soc_0_30cm']/0.3/profile_env_info[ , 'Bulk_Density_0cm'] # unit gc/kg
relation_summary[ , 'soc_30_100cm'] = relation_summary[ , 'soc_30_100cm']/0.7/profile_env_info[ , 'Bulk_Density_30cm'] # unit gc/kg
relation_summary[ , 'soc_100_200cm'] = relation_summary[ , 'soc_100_200cm']/1/profile_env_info[ , 'Bulk_Density_100cm'] # unit gc/kg
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


valid_profile_loc = valid_profile_summary[1:1000, ]


#############################################################################
# figure
#############################################################################
world_coastline = rgdal::readOGR(dsn='/Users/phoenix/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp',layer = 'ne_110m_land')
world_coastline <- fortify(world_coastline)
Map.Using = world_coastline

current_data = data.frame(profile_env_info[valid_profile_loc[ , 1], ])

pdf('./Ensemble/wosis_representative_profile_loc.pdf', width = 10, height = 5)
ggplot() +
  geom_point(data = current_data, aes(x = Lon, y = Lat), color = 'black', size = 1, alpha = 0.7, na.rm = TRUE) +
  geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.3) +
  ylim(c(-56, 80)) +
  # change the background to black and white
  theme_bw() +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 35), legend.title = element_text(size = 35))  +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 20)) +
  # add title
  labs(title = '', x = '', y = '') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 40)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30))




dev.off()

