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

max_depth = env_info_ss[ , 'Max_Depth']

stat_r2 = readMat(paste(data_dir_output, 'mcmc_summary_', model_name, '/', model_name, '_stat_r2.mat', sep = ''))
stat_r2 = apply(stat_r2$stat.r2, 1, max, na.rm = TRUE)

# site by site results
soc_stock_ss = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_site_A_', model_name, '.mat', sep = ''))
soc_stock_ss = soc_stock_ss$Var.Decom.Grid

soc_stock_tau_mean = cbind(Re(soc_stock_ss[[5]][ , 6]), Re(soc_stock_ss[[9]][ , 6]))
bulk_process_mean = cbind(Re(soc_stock_ss[[11]][ , 6]), Re(soc_stock_ss[[15]][ , 6]), Re(soc_stock_ss[[12]][ , 6]), Re(soc_stock_ss[[14]][ , 6]), Re(soc_stock_ss[[13]][ , 6]), Re(soc_stock_ss[[16]][ , 6]))


# site level predictions by PRODA
icross_valid = 1
proda_predict = c()
for (icross_valid in 1:10) { 
  data_middle = readMat(paste(data_dir_output, 'world_simulation_analyses/bulk_process_site_da_vs_proda_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_', icross_valid, '.mat', sep = ''))
  data_middle = data_middle$bulk.process.site.da.vs.proda
  
  nn_site_loc = read.table(paste(data_dir_output, 'neural_networking/nn_site_loc_', model_name, '_', time_domain, '_', nn_exp_name, '_cross_valid_', icross_valid, '.csv', sep = ''), header = FALSE)
  nn_site_loc = nn_site_loc$V1
  proda_predict = rbind(proda_predict, cbind(data_middle[[1]][ , 1], nn_site_loc, icross_valid))
}

valid_loc = proda_predict[ , 2]
bulk_process_mean = bulk_process_mean[valid_loc, ]

# valid_loc = which(is.na(ss_para_result[ , 1]) == 0 & stat_r2 > 0.75 & max_depth > 50
#                   & soc_stock_tau_mean[ , 1] > 1000
#                   & soc_stock_tau_mean[ , 1] < 1000000
#                   & soc_stock_tau_mean[ , 2] > 0
#                   & soc_stock_tau_mean[ , 2] < 10000)
# 
# bulk_process_mean = bulk_process_mean[valid_loc, ]


#################################################################################
# Plot Figures SOC content and Coarse_Fragments_v_0cm
#################################################################################
current_data = data.frame(cbind(env_info_ss[valid_loc, 'Coarse_Fragments_v_0cm'], bulk_process_mean[ , 1], proda_predict[ , 1]))
colnames(current_data) = c('var', 'cue', 'proda')

lm_model = lm(cue ~ var, data = current_data)
lm_summary = summary(lm_model)
lm_summary

corr_test = cor.test(current_data$cue, current_data$var)
corr_test$estimate

p_fragments =
  ggplot() + 
  stat_bin_hex(data = current_data, aes(x = var, y = cue), bins = 100, alpha = 0.7) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
  geom_smooth(data = current_data, aes(x = var, y = cue, ), method = 'auto', linetype = 'solid', color = '#005AB5', fill = '#005AB5', size = 3) +
  geom_smooth(data = current_data, aes(x = var, y = proda), method = 'auto', linetype = 'dashed', color = '#DC3220', fill = '#DC3220', size = 3) +
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = expression(paste('Volumetric Coarse Fragments (%)')), y = 'CUE') + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 

#################################################################################
# Plot Figures SOC content and silt fraction
#################################################################################
current_data = data.frame(cbind(env_info_ss[valid_loc, 'Silt_Content_0cm'], bulk_process_mean[ , 1], proda_predict[ , 1]))
colnames(current_data) = c('var', 'cue', 'proda')

lm_model = lm(cue ~ var, data = current_data)
lm_summary = summary(lm_model)
lm_summary

corr_test = cor.test(current_data$cue, current_data$var)
corr_test$estimate

p_silt =
  ggplot(data = current_data) + 
  # geom_boxplot(aes(x = as.factor(var), y = cue)) 
  stat_bin_hex(aes(x = var, y = cue), bins = 100, alpha = 0.7) +
  # geom_point(aes(x = var, y = cue), size = 0.1, color = 'red4') +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(trans = 'identity') + 
  geom_smooth(data = current_data, aes(x = var, y = cue), method = 'auto', linetype = 'solid', color = '#005AB5', fill = '#005AB5', size = 3) +
  geom_smooth(data = current_data, aes(x = var, y = proda), method = 'auto', linetype = 'dashed', color = '#DC3220', fill = '#DC3220', size = 3) +
  # geom_text(data = text_data, aes(x = 1250, y = 0.6, label = paste('r = ', round(corr_test$estimate, 2), '\np < 0.001', sep = '')), hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = expression(paste('Silt Fraction (%)')), y = 'CUE') + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 



#################################################################################
# Plot Figures SOC content and sand fraction
#################################################################################
current_data = data.frame(cbind(env_info_ss[valid_loc, 'Sand_Content_0cm'], bulk_process_mean[ , 1], proda_predict[ , 1]))
colnames(current_data) = c('var', 'cue', 'proda')


lm_model = lm(cue ~ var, data = current_data)
lm_summary = summary(lm_model)
lm_summary

corr_test = cor.test(current_data$cue, current_data$var)
corr_test$estimate

p_sand =
  ggplot(data = current_data) + 
  # geom_boxplot(aes(x = as.factor(var), y = cue)) 
  stat_bin_hex(aes(x = var, y = cue), bins = 100, alpha = 0.7) +
  # geom_point(aes(x = var, y = cue), size = 0.1, color = 'red4') +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(trans = 'identity') + 
  geom_smooth(data = current_data, aes(x = var, y = cue), method = 'auto', linetype = 'solid', color = '#005AB5', fill = '#005AB5', size = 3) +
  geom_smooth(data = current_data, aes(x = var, y = proda), method = 'auto', linetype = 'dashed', color = '#DC3220', fill = '#DC3220', size = 3) +
  # geom_text(data = text_data, aes(x = 1250, y = 0.6, label = paste('r = ', round(corr_test$estimate, 2), '\np < 0.001', sep = '')), hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = expression(paste('Sand Fraction (%)')), y = 'CUE') + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 



#################################################################################
# Plot Figures SOC content and clay fraction
#################################################################################
current_data = data.frame(cbind(env_info_ss[valid_loc, 'Clay_Content_0cm'], bulk_process_mean[ , 1], proda_predict[ , 1]))
colnames(current_data) = c('var', 'cue', 'proda')


lm_model = lm(cue ~ var, data = current_data)
lm_summary = summary(lm_model)
lm_summary

corr_test = cor.test(current_data$cue, current_data$var)
corr_test$estimate

p_clay =
  ggplot(data = current_data) + 
  # geom_boxplot(aes(x = as.factor(var), y = cue)) 
  stat_bin_hex(aes(x = var, y = cue), bins = 100, alpha = 0.7) +
  # geom_point(aes(x = var, y = cue), size = 0.1, color = 'red4') +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(trans = 'identity') + 
  geom_smooth(data = current_data, aes(x = var, y = cue), method = 'auto', linetype = 'solid', color = '#005AB5', fill = '#005AB5', size = 3) +
  geom_smooth(data = current_data, aes(x = var, y = proda), method = 'auto', linetype = 'dashed', color = '#DC3220', fill = '#DC3220', size = 3) +
  # geom_text(data = text_data, aes(x = 1250, y = 0.6, label = paste('r = ', round(corr_test$estimate, 2), '\np < 0.001', sep = '')), hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = expression(paste('Clay Fraction (%)')), y = 'CUE') + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


#################################################################################
# Plot Figures SOC content and bulk density
#################################################################################
current_data = data.frame(cbind(env_info_ss[valid_loc, 'Bulk_Density_0cm'], bulk_process_mean[ , 1], proda_predict[ , 1]))
colnames(current_data) = c('var', 'cue', 'proda')

lm_model = lm(cue ~ var, data = current_data)
lm_summary = summary(lm_model)
lm_summary

corr_test = cor.test(current_data$cue, current_data$var)
corr_test$estimate

p_bd =
  ggplot(data = current_data) + 
  # geom_boxplot(aes(x = as.factor(var), y = cue)) 
  stat_bin_hex(aes(x = var, y = cue), bins = 100, alpha = 0.7) +
  # geom_point(aes(x = var, y = cue), size = 0.1, color = 'red4') +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
  scale_x_continuous(trans = 'identity') +
  scale_y_continuous(trans = 'identity') + 
  geom_smooth(data = current_data, aes(x = var, y = cue), method = 'auto', linetype = 'solid', color = '#005AB5', fill = '#005AB5', size = 3) +
  geom_smooth(data = current_data, aes(x = var, y = proda), method = 'auto', linetype = 'dashed', color = '#DC3220', fill = '#DC3220', size = 3) +
  # geom_text(data = text_data, aes(x = 1250, y = 0.6, label = paste('r = ', round(corr_test$estimate, 2), '\np < 0.001', sep = '')), hjust = 0, size = 10, show.legend = FALSE) + 
  # change the background to black and white
  theme_classic() + 
  # add title
  labs(title = '', x = expression(paste('Bulk Density (kg', ' m'^'-3', ')')), y = 'CUE') + 
  # change the legend properties
  guides(fill = guide_colorbar(barwidth = 3, barheight = 10)) + 
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))  +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) 


jpeg(paste('./Ensemble/cue_vs_soil_var_CESM2.jpeg', sep = ''), width = 24, height = 16, units = 'in', res = 300)
plot_grid(p_bd, p_fragments, p_clay, p_sand, p_silt, 
          nrow = 2, ncol = 3)
dev.off()

