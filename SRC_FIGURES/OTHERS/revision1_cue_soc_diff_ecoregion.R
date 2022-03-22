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

#############################################################################
# microbial data path
#############################################################################
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


#############################################################################
# microbial profile info
#############################################################################
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
var_climate[which(global_climate == 2)] = 101 # Am
var_climate[which(global_climate == 3)] = 101 # Aw
var_climate[which(global_climate >= 4 & global_climate <= 5)] = 102 # BwX
var_climate[which(global_climate >= 6 & global_climate <= 7)] = 102 # BsW
var_climate[which(global_climate >= 8 & global_climate <= 10)] = 103 # CsX
var_climate[which(global_climate >= 11 & global_climate <= 13)] = 103 # CwX
var_climate[which(global_climate >= 14 & global_climate <= 16)] = 103 # CfX
var_climate[which(global_climate >= 17 & global_climate <= 20)] = 104 # DsX
var_climate[which(global_climate >= 21 & global_climate <= 24)] = 114 # DwX
var_climate[which(global_climate >= 25 & global_climate <= 28)] = 114 # DfX
var_climate[which(global_climate >= 29 & global_climate <= 30)] = 114 # E

# load ESA land cover
global_landcover = profile_env_info[ , 'ESA_Land_Cover']
var_landcover = array(NA, dim = c(length(global_landcover), 1))

var_landcover[which(global_landcover >= 1 & global_landcover <= 4)] = 101 # agriculture
var_landcover[which((var_climate <= 111 & var_climate >= 109) & ((global_landcover >= 5 & global_landcover <= 10) | (global_landcover >= 16 & global_landcover <= 17)))] = 102 # boreal forest
var_landcover[which((var_climate <= 108 | var_climate >= 112) & ((global_landcover >= 5 & global_landcover <= 10) | (global_landcover >= 16 & global_landcover <= 17)))] = 102 # other forest
var_landcover[which(global_landcover >= 11 & global_landcover <= 15)] = 103 # grassland & shrubland 
var_landcover[which(global_landcover == 18)] = 104 # wetland
var_landcover[which(global_landcover >=19)] = 101 # urban and other

# load soil texture
global_texture = profile_env_info[ , 'Texture_USDA_0cm']
var_texture = array(NA, dim = c(length(global_texture), 1))

var_texture[which(global_texture >= 1 & global_texture <= 3)] = 101 # clay
var_texture[which(global_texture >= 4 & global_texture <= 9)] = 102 # loam
var_texture[which(global_texture >= 10 & global_texture <= 12)] = 103 # sand

# load soil pH
global_ph = profile_env_info[ , 'pH_Water_0cm']
var_ph = array(NA, dim = c(length(global_ph), 1))

var_ph[which(global_ph >= 0 & global_ph <= 61)] = 101 # acid
var_ph[which(global_ph >= 61 & global_ph <= 78)] = 102 # neutral
var_ph[which(global_ph >= 78 & global_ph <= 100)] = 103 # alkaline
#############################################################################
# microbial load da summary
#############################################################################
col_list_relation_summary = c('profile_num',
                              'cue', 'mm_ratio_assim', 'mm_ratio_decom', 'tau4s2_enz', 'allo_slope_mic',
                              'soc_0_30cm', 'soc_30_100cm', 'soc_100_200cm', 'mic_0_30cm', 'mic_stock', 'enz_stock', 'doc_stock', 'poc_stock', 'soc_total', 'r2')

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
koppen_climate_class$class2[which(koppen_climate_class$class1 == 2)] = 101 # Am
koppen_climate_class$class2[which(koppen_climate_class$class1 == 3)] = 101 # Aw
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 4 & koppen_climate_class$class1 <= 5)] = 102 # BwX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 6 & koppen_climate_class$class1 <= 7)] = 102 # BsW
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 8 & koppen_climate_class$class1 <= 10)] = 103 # CsX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 11 & koppen_climate_class$class1 <= 13)] = 103 # CwX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 14 & koppen_climate_class$class1 <= 16)] = 103 # CfX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 17 & koppen_climate_class$class1 <= 20)] = 104 # DsX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 21 & koppen_climate_class$class1 <= 24)] = 104 # DwX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 25 & koppen_climate_class$class1 <= 28)] = 104 # DfX
koppen_climate_class$class2[which(koppen_climate_class$class1 >= 29 & koppen_climate_class$class1 <= 30)] = 104 # E

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
# regressed model summary
#############################################################################

# climate 
current_data = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm', 'cue')], var_climate[valid_profile_loc[ , 1]]))
colnames(current_data) = c('soc', 'cue', 'eco_region')
current_data = current_data[which(is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0), ]

fit_summary_climate_mic = array(NA, dim = c(4, 5))

iceo = 1
for (iceo in 1:4) {
  fit_model = summary(lm(log10(soc) ~ cue, data = current_data[current_data$eco_region == (100+iceo), ]))
  fit_summary_climate_mic[iceo, ] = c(fit_model$coefficients[2, c(1, 2, 4)], fit_model$df[2], iceo)
  
}

# land cover 
current_data = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm', 'cue')], var_landcover[valid_profile_loc[ , 1]]))
colnames(current_data) = c('soc', 'cue', 'eco_region')
current_data = current_data[which(is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0), ]

fit_summary_veg_mic = array(NA, dim = c(4, 5))

iceo = 1
for (iceo in 1:4) {
  fit_model = summary(lm(log10(soc) ~ cue, data = current_data[current_data$eco_region == (100+iceo), ]))
  fit_summary_veg_mic[iceo, ] = c(fit_model$coefficients[2, c(1, 2, 4)], fit_model$df[2], iceo)
}

# soil texture
current_data = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm', 'cue')], var_texture[valid_profile_loc[ , 1]]))
colnames(current_data) = c('soc', 'cue', 'eco_region')
current_data = current_data[which(is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0), ]

fit_summary_texture_mic = array(NA, dim = c(3, 5))

iceo = 1
for (iceo in 1:3) {
  fit_model = summary(lm(log10(soc) ~ cue, data = current_data[current_data$eco_region == (100+iceo), ]))
  fit_summary_texture_mic[iceo, ] = c(fit_model$coefficients[2, c(1, 2, 4)], fit_model$df[2], iceo)
}


# soil pH
current_data = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm', 'cue')], var_ph[valid_profile_loc[ , 1]]))
colnames(current_data) = c('soc', 'cue', 'eco_region')
current_data = current_data[which(is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0), ]

fit_summary_ph_mic = array(NA, dim = c(3, 5))

iceo = 4
for (iceo in 1:3) {
  fit_model = summary(lm(log10(soc) ~ cue, data = current_data[current_data$eco_region == (100+iceo), ]))
  fit_summary_ph_mic[iceo, ] = c(fit_model$coefficients[2, c(1, 2, 4)], fit_model$df[2], iceo)
}

#############################################################################
# linear data path
#############################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'


#################################################################################
# Load Projected SOC site by site
#################################################################################
# site by site soc stock
soc_stock_ss = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_site_A_', model_name, '.mat', sep = ''))
soc_stock_ss = soc_stock_ss$Var.Decom.Grid

soc_stock_tau_mean = cbind(Re(soc_stock_ss[[5]][ , 6]), Re(soc_stock_ss[[9]][ , 6]))
bulk_process_mean = cbind(Re(soc_stock_ss[[11]][ , 6]), Re(soc_stock_ss[[15]][ , 6]), Re(soc_stock_ss[[12]][ , 6]), Re(soc_stock_ss[[14]][ , 6]), Re(soc_stock_ss[[13]][ , 6]), Re(soc_stock_ss[[16]][ , 6]))

soc_content_0_30cm = Re(soc_stock_ss[[2]][ , 6])/0.3/profile_env_info[ , 'Bulk_Density_0cm'] # unit gC/kg
soc_tau_0_30cm = Re(soc_stock_ss[[6]][ , 6])

soc_content_30_100cm = (Re(soc_stock_ss[[3]][ , 6]) - Re(soc_stock_ss[[2]][ , 6]))/0.7/profile_env_info[ , 'Bulk_Density_30cm']  # unit gC/kg
soc_tau_30_100cm = Re(soc_stock_ss[[7]][ , 6]) - Re(soc_stock_ss[[6]][ , 6])

soc_content_100_200cm = (Re(soc_stock_ss[[4]][ , 6]) - Re(soc_stock_ss[[3]][ , 6]))/1/profile_env_info[ , 'Bulk_Density_100cm'] # unit gC/kg
soc_tau_100_200cm = Re(soc_stock_ss[[8]][ , 6]) - Re(soc_stock_ss[[7]][ , 6])


valid_loc = readMat(paste(data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info_', model_name, '_', time_domain, '_maxmin_scaled.mat', sep = ''))
valid_loc = valid_loc$profile.env.info[ , 1]

# valid_loc = valid_profile_loc[ , 1]



#############################################################################
# regressed model summary
#############################################################################

# climate 
current_data = data.frame(cbind(bulk_process_mean[valid_loc, 1], soc_content_0_30cm[valid_loc], var_climate[valid_loc]))
colnames(current_data) = c('cue', 'soc', 'eco_region')
current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]

fit_summary_climate_clm5 = array(NA, dim = c(4, 5))

iceo = 1
for (iceo in 1:4) {
  fit_model = summary(lm(log10(soc) ~ cue, data = current_data[current_data$eco_region == (100+iceo), ]))
  fit_summary_climate_clm5[iceo, ] = c(fit_model$coefficients[2, c(1, 2, 4)], fit_model$df[2], iceo)
  
}

# land cover 
current_data = data.frame(cbind(bulk_process_mean[valid_loc, 1], soc_content_0_30cm[valid_loc], var_landcover[valid_loc]))
colnames(current_data) = c('cue', 'soc', 'eco_region')
current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]

fit_summary_veg_clm5 = array(NA, dim = c(4, 5))

iceo = 1
for (iceo in 1:4) {
  fit_model = summary(lm(log10(soc) ~ cue, data = current_data[current_data$eco_region == (100+iceo), ]))
  fit_summary_veg_clm5[iceo, ] = c(fit_model$coefficients[2, c(1, 2, 4)], fit_model$df[2], iceo)
}

# soil texture
current_data = data.frame(cbind(bulk_process_mean[valid_loc, 1], soc_content_0_30cm[valid_loc], var_texture[valid_loc]))
colnames(current_data) = c('cue', 'soc', 'eco_region')
current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]

fit_summary_texture_clm5 = array(NA, dim = c(3, 5))

iceo = 1
for (iceo in 1:3) {
  fit_model = summary(lm(log10(soc) ~ cue, data = current_data[current_data$eco_region == (100+iceo), ]))
  fit_summary_texture_clm5[iceo, ] = c(fit_model$coefficients[2, c(1, 2, 4)], fit_model$df[2], iceo)
}

# soil pH
current_data = data.frame(cbind(bulk_process_mean[valid_loc, 1], soc_content_0_30cm[valid_loc], var_ph[valid_loc]))
colnames(current_data) = c('cue', 'soc', 'eco_region')
current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]

fit_summary_ph_clm5 = array(NA, dim = c(3, 5))

iceo = 1
for (iceo in 1:3) {
  fit_model = summary(lm(log10(soc) ~ cue, data = current_data[current_data$eco_region == (100+iceo), ]))
  fit_summary_ph_clm5[iceo, ] = c(fit_model$coefficients[2, c(1, 2, 4)], fit_model$df[2], iceo)
}

#################################################################################
# plot figure
#################################################################################
# color_scheme = c('#DC3220', '#005AB5')
color_scheme = c('#005AB5')

# soil texture
# current_data = data.frame(rbind(cbind(fit_summary_texture_clm5, 1), 
#                                 cbind(fit_summary_texture_mic, 2)))
current_data = data.frame(rbind(cbind(fit_summary_texture_clm5, 1)))
colnames(current_data) = c('slope', 'var', 'pvalue', 'df', 'eco_region', 'model')
current_data$p_indice = current_data$pvalue
current_data$p_indice[current_data$pvalue < 0.005] = 1
current_data$p_indice[current_data$pvalue >= 0.005] = 2

p_texture =
  ggplot(data = current_data, aes(x = reorder(as.factor(eco_region), slope), y = slope, group = as.factor(model), linetype = as.factor(p_indice), fill = as.factor(model))) + 
  geom_bar(color = 'black', width = 0.6, alpha = 1, size = 1, stat = 'identity', position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = slope - var, ymax = slope + var), width = 0.3, size = 1, stat = 'identity', position = position_dodge(0.9)) +
  geom_text(aes(y = slope + var, label = df), size = 7, hjust = -0.2, vjust = 0.5, position = position_dodge(0.9)) +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(limits = c(NA, 8)) +
  scale_linetype_manual(name = '', values = c('solid', 'dashed'), labels = c('P < 0.005', 'P > 0.005')) + 
  scale_x_discrete(name = '', labels = c('Clay', '    Loam    ', 'Sand')) +
  scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5', 'MIC')) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
  theme(axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) +
  # add title
  labs(title = 'Soil Texture', x = '', y = 'Slope of CUE-SOC regression') + 
  # legend
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20))  +
  theme(legend.justification = c(1, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 20), axis.title = element_text(size = 20), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 


# land cover
current_data = data.frame(rbind(cbind(fit_summary_veg_clm5, 1)))
# current_data = data.frame(rbind(cbind(fit_summary_veg_clm5, 1), 
                                # cbind(fit_summary_veg_mic, 2)))
colnames(current_data) = c('slope', 'var', 'pvalue', 'df', 'eco_region', 'model')
current_data$p_indice = current_data$pvalue
current_data$p_indice[current_data$pvalue < 0.005] = 1
current_data$p_indice[current_data$pvalue >= 0.005] = 2

p_veg =
  ggplot(data = current_data, aes(x = reorder(as.factor(eco_region), slope), y = slope, group = as.factor(model), linetype = as.factor(p_indice), fill = as.factor(model))) + 
  geom_bar(color = 'black', width = 0.6, alpha = 1, size = 1, stat = 'identity', position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = slope - var, ymax = slope + var), width = 0.3, size = 1, stat = 'identity', position = position_dodge(0.9)) +
  geom_text(aes(y = slope + var, label = df), size = 7, hjust = -0.2, vjust = 0.5, position = position_dodge(0.9)) +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(limits = c(NA, 10)) +
  scale_linetype_manual(name = '', values = c('solid', 'dashed'), labels = c('P < 0.005', 'P > 0.005')) + 
  scale_x_discrete(name = '', labels = c('Agri. & Urban', 'Grass & Shrub', 'Forest', 'Wetland')) +
  scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5', 'MIC')) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
  theme(axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) +
  # add title
  labs(title = 'Land Cover', x = '', y = 'Slope of CUE-SOC regression') + 
  # legend
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 20), axis.title = element_text(size = 20), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 


# climate
# current_data = data.frame(rbind(cbind(fit_summary_climate_clm5, 1), 
#                                 cbind(fit_summary_climate_mic, 2)))
current_data = data.frame(rbind(cbind(fit_summary_climate_clm5, 1)))
colnames(current_data) = c('slope', 'var', 'pvalue', 'df', 'eco_region', 'model')
current_data$p_indice = current_data$pvalue
current_data$p_indice[current_data$pvalue < 0.005] = 1
current_data$p_indice[current_data$pvalue >= 0.005] = 2

# p_climate =
#   ggplot(data = current_data, aes(x = reorder(as.factor(eco_region), slope), y = slope, group = as.factor(model), linetype = as.factor(p_indice), fill = as.factor(model))) + 
#   geom_bar(color = 'black', width = 0.6, alpha = 1, size = 1, stat = 'identity', position = position_dodge(0.9)) +
#   geom_errorbar(aes(ymin = slope - var, ymax = slope + var), width = 0.3, size = 1, stat = 'identity', position = position_dodge(0.9)) +
#   geom_text(aes(y = slope + var, label = df), size = 7, hjust = -0.2, vjust = 0.5, position = position_dodge(0.9)) +
#   theme_classic() +
#   coord_flip() +
#   scale_y_continuous(limits = c(NA, 8)) +
#   scale_linetype_manual(name = '', values = c('solid', 'dashed'), labels = c('P < 0.005', 'P > 0.005')) + 
#   scale_x_discrete(name = '', labels = c('Boreal', 'Arid', 'Temperate', 'Tropics')) +
#   scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5', 'MIC')) +
#   theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
#   theme(axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) +
#   # add title
#   labs(title = 'Climate', x = '', y = 'Regressed CUE-SOC slope') + 
#   # legend
#   theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20))  +
#   theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background = element_rect(fill = NA)) + 
#   theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
#   # modify the position of title
#   theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
#   # modify the margin
#   theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
#   theme(axis.text=element_text(size = 20), axis.title = element_text(size = 20), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 

p_climate =
  ggplot(data = current_data, aes(x = reorder(as.factor(eco_region), slope), y = slope, group = as.factor(model), linetype = as.factor(p_indice))) + 
    geom_bar(color = 'black', fill = color_scheme, width = 0.6, alpha = 1, size = 1, stat = 'identity', position = position_dodge(0.9)) +
    geom_errorbar(aes(ymin = slope - var, ymax = slope + var), width = 0.3, size = 1, stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(y = slope + var, label = df), size = 7, hjust = -0.2, vjust = 0.5, position = position_dodge(0.9)) +
    theme_classic() +
    coord_flip() +
    scale_y_continuous(limits = c(NA, 8)) +
    scale_linetype_manual(name = '', values = c('solid', 'dashed'), labels = c('P < 0.005', 'P > 0.005')) + 
    scale_x_discrete(name = '', labels = c('Boreal', 'Arid', 'Temperate', 'Tropics')) +
    scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5', 'MIC')) +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
    theme(axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) +
    # add title
    labs(title = 'Climate', x = '', y = 'Slope of CUE-SOC regression') + 
    # legend
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20))  +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background = element_rect(fill = NA)) + 
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
    # modify the margin
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
    theme(axis.text=element_text(size = 20), axis.title = element_text(size = 20), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 
  
# pH
# current_data = data.frame(rbind(cbind(fit_summary_ph_clm5, 1), 
#                                 cbind(fit_summary_ph_mic, 2)))
current_data = data.frame(rbind(cbind(fit_summary_ph_clm5, 1)))
colnames(current_data) = c('slope', 'var', 'pvalue', 'df', 'eco_region', 'model')
current_data$p_indice = current_data$pvalue
current_data$p_indice[current_data$pvalue < 0.005] = 1
current_data$p_indice[current_data$pvalue >= 0.005] = 2

p_ph =
  ggplot(data = current_data, aes(x = reorder(as.factor(eco_region), slope), y = slope, group = as.factor(model), linetype = as.factor(p_indice), fill = as.factor(model))) + 
  geom_bar(color = 'black', width = 0.6, alpha = 1, size = 1, stat = 'identity', position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = slope - var, ymax = slope + var), width = 0.3, size = 1, stat = 'identity', position = position_dodge(0.9)) +
  geom_text(aes(y = slope + var, label = df), size = 7, hjust = -0.2, vjust = 0.5, position = position_dodge(0.9)) +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(limits = c(NA, 8)) +
  scale_linetype_manual(name = '', values = c('solid', 'dashed'), labels = c('P < 0.005', 'P > 0.005')) + 
  scale_x_discrete(name = '', labels = c('Alkaline', '      Neutral      ', 'Acid')) +
  scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5', 'MIC')) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
  theme(axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) +
  # add title
  labs(title = 'Soil Chemical', x = '', y = 'Slope of CUE-SOC regression') + 
  # legend
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 20), axis.title = element_text(size = 20), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.1, 'inch')) 

p_texture = p_texture + 
  labs(y = '  ')
p_ph = p_ph + 
  labs(y = '  ')

jpeg(paste('./Ensemble/cue_soc_slope_eco_regions.jpeg', sep = ''), width = 15, height = 12, units = 'in', res = 300)
plot_grid( p_texture, p_ph, p_climate, p_veg,
           labels = c('a', 'b', 'c', 'd'),
           label_size = 35,
           label_x = 0.02, label_y = 1,
           label_fontfamily = 'Arial',
           label_fontface = 'bold',
           nrow = 2
)
dev.off()
