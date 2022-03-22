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

library(ppcor)
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

para_summary = readMat(paste(data_dir_output, '/da_summary_', model_name, '/', model_name, '_da_summary_mic_para.mat', sep = ''))
para_summary = para_summary$da.summary.mic.para
colnames(para_summary) = para_name
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
# microbial figure
#############################################################################
current_data = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_0_30cm', 'cue')], valid_profile_loc))
colnames(current_data) = c('soc', 'cue', 'profile_num', 'climate', 'sample')
#---------------mix model---------------#
mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data)
mix_model_summary_1 = summary(mix_model)
mix_model_summary_1
fit_function_1 = function(x) {10**(mix_model_summary_1$coefficients[1, 1] + x*(mix_model_summary_1$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc[is.na(current_data$soc) == 0]) ~ predict(mix_model)))
r_squared_mixed 
text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(1), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_1$coefficients[1, 1], 2), ' + ', round(mix_model_summary_1$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_1$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_1$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))

color_scheme = c('#b30000', '#e34a33', '#fc8d59', '#fdbb84', '#fdd49e', '#fef0d9', '#f0f9e8', '#ccebc5', '#a8ddb5', '#7bccc4', '#43a2ca', '#0868ac')


p_mic_cue_soc_0_30cm =
ggplot(data = current_data) + 
  geom_point(data = current_data, aes(x = cue, y = soc, color = climate), size = 6, alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_1, size = 2, color = 'black') + 
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
  # stat_bin_hex(aes(x = cue, y = soc_0_30cm), bins = 100) +
  scale_color_stepsn(name = 'Climate', colors = color_scheme, trans = 'identity', limits = c(100, 112), breaks = c(100:112), labels = c('', 'Tropics', '', '', '', '', '', '', '', '', '', '', 'Polar')) +
  scale_x_continuous(limits = c(0, max(current_data$cue, na.rm = TRUE)), trans = 'identity') +
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
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#----------------------------------------------------
# cue soc relation 30 - 100 cm
#----------------------------------------------------
current_data = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_30_100cm', 'cue')], valid_profile_loc))
colnames(current_data) = c('soc', 'cue', 'profile_num', 'climate', 'sample')
#---------------mix model---------------#
mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data)
mix_model_summary_2 = summary(mix_model)
mix_model_summary_2
fit_function_2 = function(x) {10**(mix_model_summary_2$coefficients[1, 1] + x*(mix_model_summary_2$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc[is.na(current_data$soc) == 0]) ~ predict(mix_model)))
r_squared_mixed 
text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(1000), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_2$coefficients[1, 1], 2), ' + ', round(mix_model_summary_2$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_2$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_2$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


p_mic_cue_soc_30_100cm =
  ggplot(data = current_data) + 
  geom_point(data = current_data, aes(x = cue, y = soc, color = climate), size = 6, alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_2, size = 2, color = 'black') + 
  # geom_smooth(data = current_data, aes(x = cue, y = soc_0_30cm), method = 'lm', color = 'black', fill = 'black', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
  # stat_bin_hex(aes(x = cue, y = soc_0_30cm), bins = 100) +
  scale_color_stepsn(name = 'Climate', colors = color_scheme, trans = 'identity', limits = c(100, 112), breaks = c(100:112), labels = c('', 'Tropics', '', '', '', '', '', '', '', '', '', '', 'Polar')) +
  scale_x_continuous(limits = c(0, max(current_data$cue, na.rm = TRUE)), trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(color = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'top', title.hjust = 0, label.hjust = 0, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#----------------------------------------------------
# cue soc relation 100 - max cm
#----------------------------------------------------
current_data = data.frame(cbind(relation_summary[valid_profile_loc[ , 1], c('soc_100_200cm', 'cue')], valid_profile_loc))
colnames(current_data) = c('soc', 'cue', 'profile_num', 'climate', 'sample')
#---------------mix model---------------#
mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data)
mix_model_summary_3 = summary(mix_model)
mix_model_summary_3
fit_function_3 = function(x) {10**(mix_model_summary_3$coefficients[1, 1] + x*(mix_model_summary_3$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc[is.na(current_data$soc) == 0]) ~ predict(mix_model)))
r_squared_mixed 
text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(1000), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_3$coefficients[1, 1], 2), ' + ', round(mix_model_summary_3$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_3$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_3$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


p_mic_cue_soc_100_200cm =
  ggplot(data = current_data) + 
  geom_point(data = current_data, aes(x = cue, y = soc, color = climate), size = 6, alpha = 0.7, na.rm = TRUE) +
  geom_function(fun = fit_function_3, size = 2, color = 'black') + 
  # geom_smooth(data = current_data, aes(x = cue, y = soc_0_30cm), method = 'lm', color = 'black', fill = 'black', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
  # stat_bin_hex(aes(x = cue, y = soc_0_30cm), bins = 100) +
  scale_color_stepsn(name = 'Climate', colors = color_scheme, trans = 'identity', limits = c(100, 112), breaks = c(100:112), labels = c('', 'Tropics', '', '', '', '', '', '', '', '', '', '', 'Polar')) +
  scale_x_continuous(limits = c(0, max(current_data$cue, na.rm = TRUE)), trans = 'identity') +
  scale_y_continuous(limits = c(0.01, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(color = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'top', title.hjust = 0, label.hjust = 0, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#----------------------------------------------------
# cue mat relation
#----------------------------------------------------
current_data = data.frame(relation_summary[valid_profile_loc[ , 1], c('cue')], profile_env_info[valid_profile_loc[ , 1], 'BIO1'])
colnames(current_data) = c('cue', 'mat')

corr_test = cor.test(current_data$mat, current_data$cue)
corr_test
text_data = data.frame('x_axis' = c(min(current_data$mat, na.rm = TRUE)), 
                       'y_axis' = c(0.6), 
                       'equation' = paste('r = ', round(corr_test$estimate, 2), '\nP = ', round(corr_test$p.value, 2), sep = ''))

summary(lm(current_data$cue ~ current_data$mat))

p_mic_cue_mat =
  ggplot(data = current_data) + 
  geom_point(data = current_data, aes(x = mat, y = cue), color = 'black', size = 6, alpha = 0.7, na.rm = TRUE) +
  geom_smooth(data = current_data, aes(x = mat, y = cue), method = 'lm', color = 'black', fill = 'black', size = 2) +
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 12, show.legend = FALSE) + 
  scale_x_continuous(limits = c(min(current_data$mat, na.rm = TRUE), max(current_data$mat, na.rm = TRUE)), trans = 'identity') +
  scale_y_continuous(limits = c(0, max(current_data$cue, na.rm = TRUE)), trans = 'identity') +
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Mean annual temperature (°C)', y = 'Carbon use efficiency') + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#----------------------------------------------------
# michaelis menten ratios histograms
#----------------------------------------------------
color_scheme = c('#DC3220', '#005AB5')

jpeg(paste('./Ensemble/michaelis_menten_ratio_', model_name, '.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)
current_data = data.frame(rbind(cbind(relation_summary[valid_profile_loc[ , 1], 'mm_ratio_assim'], 1), cbind(cbind(relation_summary[valid_profile_loc[ , 1], 'mm_ratio_decom'], 2))))

colnames(current_data) = c('ratio', 'class')
ggplot() +
  geom_histogram(data = current_data, aes(x = ratio, y = ..density.., color = as.factor(class), fill = as.factor(class)), position = 'identity', alpha = 0.3, na.rm = TRUE, size = 1, bins = 30) + 
  geom_density(data = current_data, aes(x = ratio, color = as.factor(class)), fill = NA, alpha = 0.5, na.rm = TRUE, size = 3) + 
  geom_vline(xintercept = c(0.01, 100), linetype = 'longdash', size = 2, color = 'grey') + 
  scale_x_continuous(trans = 'log10', breaks = c(0.01, 1, 100, 10000), labels = trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(name = '', values = color_scheme, labels = c('Assimilation', 'Decomposition')) +
  scale_color_manual(name = '', values = color_scheme, labels = c('Assimilation', 'Decomposition')) +
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Km/[Substrate]', y = 'Density') + 
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


#################################################################################
# linear plot figure
#################################################################################
#----------------------------------------------------
# Plot Figures soc and cue (0 - 30cm)
#----------------------------------------------------
current_data = data.frame(cbind(bulk_process_mean[valid_loc, 1], soc_content_0_30cm[valid_loc], var_climate[valid_loc]))
colnames(current_data) = c('cue', 'soc', 'climate')
current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]
current_data = current_data[current_data$soc != 0, ]

#------------------------- mixed model
mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data)
mix_model_summary_soc_linear_1 = summary(mix_model)
mix_model_summary_soc_linear_1
fit_function_soc_linear_1 = function(x) {10**(mix_model_summary_soc_linear_1$coefficients[1, 1] + x*(mix_model_summary_soc_linear_1$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc) ~ predict(mix_model)))
r_squared_mixed 
text_data = data.frame('x_axis' = c(0.25), 
                       'y_axis' = c(1), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_soc_linear_1$coefficients[1, 1], 2), ' + ', round(mix_model_summary_soc_linear_1$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_soc_linear_1$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_soc_linear_1$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))

p_linear_cue_soc_0_30cm =
  ggplot(data = current_data) + 
  stat_bin_hex(aes(x = cue, y = soc), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
  scale_x_continuous(limits = c(0.25, 0.65), trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  geom_function(fun = fit_function_soc_linear_1, size = 2, color = 'black') +  
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
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
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


#----------------------------------------------------
# Plot Figures soc and cue (30 - 100cm)
#----------------------------------------------------
current_data = data.frame(cbind(bulk_process_mean[valid_loc, 1], soc_content_30_100cm[valid_loc], var_climate[valid_loc]))
colnames(current_data) = c('cue', 'soc', 'climate')
current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]
current_data = current_data[current_data$soc != 0, ]

#------------------------- mixed model
mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data)
mix_model_summary_soc_linear_2 = summary(mix_model)
mix_model_summary_soc_linear_2
fit_function_soc_linear_2 = function(x) {10**(mix_model_summary_soc_linear_2$coefficients[1, 1] + x*(mix_model_summary_soc_linear_2$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc) ~ predict(mix_model)))
r_squared_mixed 
text_data = data.frame('x_axis' = c(0.25), 
                       'y_axis' = c(1000), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_soc_linear_2$coefficients[1, 1], 2), ' + ', round(mix_model_summary_soc_linear_2$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_soc_linear_2$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_soc_linear_2$coefficients[2, 5], 3), '\nExplained Variation. = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


p_linear_cue_soc_30_100cm =
  ggplot(data = current_data) + 
  stat_bin_hex(aes(x = cue, y = soc), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
  scale_x_continuous(limits = c(0.25, 0.65), trans = 'identity') +
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  geom_function(fun = fit_function_soc_linear_2, size = 2, color = 'black') + 
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#----------------------------------------------------
# Plot Figures soc and cue (100 - 200cm)
#----------------------------------------------------
current_data = data.frame(cbind(bulk_process_mean[valid_loc, 1], soc_content_100_200cm[valid_loc], var_climate[valid_loc]))
colnames(current_data) = c('cue', 'soc', 'climate')
current_data = current_data[is.na(apply(current_data, 1, sum, na.rm = FALSE)) == 0, ]
current_data = current_data[current_data$soc != 0, ]

#------------------------- mixed model
mix_model = lmer(log10(soc) ~ cue + (1 | climate), data = current_data)
mix_model_summary_soc_linear_3 = summary(mix_model)
mix_model_summary_soc_linear_3
fit_function_soc_linear_3 = function(x) {10**(mix_model_summary_soc_linear_3$coefficients[1, 1] + x*(mix_model_summary_soc_linear_3$coefficients[2, 1]))}

r_squared_mixed = summary(lm(log10(current_data$soc) ~ predict(mix_model)))
r_squared_mixed 
text_data = data.frame('x_axis' = c(0.25), 
                       'y_axis' = c(1000), 
                       'equation' = paste('log10(SOC) = ', round(mix_model_summary_soc_linear_3$coefficients[1, 1], 2), ' + ', round(mix_model_summary_soc_linear_3$coefficients[2, 1], 2), '*CUE', '\nP(Intercept) = ', round(mix_model_summary_soc_linear_3$coefficients[1, 5], 3), ', P(CUE) = ', round(mix_model_summary_soc_linear_3$coefficients[2, 5], 3), '\nExplained Variation = ', round(r_squared_mixed$r.squared*100, 2), '%', sep = ''))


p_linear_cue_soc_100_200cm =
  ggplot(data = current_data) + 
  stat_bin_hex(aes(x = cue, y = soc), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
  scale_x_continuous(limits = c(0.25, 0.65), trans = 'identity') +
  scale_y_continuous(limits = c(0.01, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  geom_function(fun = fit_function_soc_linear_3, size = 2, color = 'black') + 
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 8.5, show.legend = FALSE) + 
  theme_classic() + 
  # add title
  labs(title = '', x = 'Carbon use efficiency', y = expression(paste('Soil organic carbon (g C kg'^'-1', ')', sep = ''))) + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


#################################################################################
# Plot Figures CUE vs. MAT
#################################################################################
current_data = data.frame(cbind(bulk_process_mean[valid_loc, 1], profile_env_info[valid_loc, 'BIO1']))
colnames(current_data) = c('cue', 'mat')
lm_model = lm(cue ~ mat, data = current_data)
lm_summary = summary(lm_model)
lm_summary

corr_test = cor.test(current_data$cue, current_data$mat)


text_data = data.frame('x_axis' = c(15), 
                       'y_axis' = c(0.6), 
                       'equation' = paste('r = ', round(corr_test$estimate, 2), '\nP = ', round(corr_test$p.value, 3), sep = ''))

p_linear_cue_mat = 
  ggplot(data = current_data) + 
  stat_bin_hex(aes(x = mat, y = cue), bins = 100) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(trans = 'identity') + 
  scale_y_continuous(trans = 'identity') +
  geom_smooth(data = current_data, aes(x = mat, y = cue), method = 'lm', color = 'black', fill = 'blue3', size = 2) + 
  geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = 'Mean annual temperature (°C)', y = 'Carbon use efficiency') + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


jpeg(paste('./Ensemble/cue_mat_relation.jpeg', sep = ''), width = 20, height = 10, units = 'in', res = 300)
plot_grid(p_mic_cue_mat, p_linear_cue_mat,
          labels = c('a', 'b'),
          label_size = 60,
          label_x = 0.02, label_y = 1,
          label_fontfamily = 'Arial',
          label_fontface = 'bold',
          nrow = 1)
dev.off()



jpeg(paste('./Ensemble/cue_soc_relation_depth.jpeg', sep = ''), width = 24, height = 16, units = 'in', res = 300)
p_mic_cue_soc_0_30cm = p_mic_cue_soc_0_30cm + 
  labs(title = '0 - 30cm', x = '  ')
p_mic_cue_soc_30_100cm = p_mic_cue_soc_30_100cm + 
  labs(title = '30 - 100cm', x = '  ', y = '  ')
p_mic_cue_soc_100_200cm  = p_mic_cue_soc_100_200cm + 
  labs(title = '> 100cm', x = '  ', y = '  ')

p_linear_cue_soc_30_100cm = p_linear_cue_soc_30_100cm + 
  labs(y = '  ')
p_linear_cue_soc_100_200cm = p_linear_cue_soc_100_200cm + 
  labs(y = '  ')


plot_grid(  p_mic_cue_soc_0_30cm, p_mic_cue_soc_30_100cm, p_mic_cue_soc_100_200cm,
            p_linear_cue_soc_0_30cm, p_linear_cue_soc_30_100cm, p_linear_cue_soc_100_200cm,
            labels = c('a', 'b', 'c', 
                       'd', 'e', 'f'),
            label_size = 60,
            label_x = 0.07, label_y = 1.02,
            label_fontfamily = 'Arial',
            label_fontface = 'bold',
            nrow = 2)
dev.off()

