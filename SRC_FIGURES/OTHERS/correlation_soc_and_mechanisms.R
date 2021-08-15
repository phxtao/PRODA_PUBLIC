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

# site by site soc stock
soc_stock_ss = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_site_A_', model_name, '.mat', sep = ''))
soc_stock_ss = soc_stock_ss$Var.Decom.Grid

soc_stock_tau_mean = cbind(Re(soc_stock_ss[[5]][ , 6]), Re(soc_stock_ss[[9]][ , 6]))
bulk_process_mean = cbind(Re(soc_stock_ss[[11]][ , 6]), Re(soc_stock_ss[[15]][ , 6]), Re(soc_stock_ss[[12]][ , 6]), Re(soc_stock_ss[[14]][ , 6]), Re(soc_stock_ss[[13]][ , 6]), Re(soc_stock_ss[[16]][ , 6]))

soc_content_0_30cm = Re(soc_stock_ss[[2]][ , 6])
soc_tau_0_30cm = Re(soc_stock_ss[[6]][ , 6])

soc_content_30_100cm = Re(soc_stock_ss[[3]][ , 6]) - Re(soc_stock_ss[[2]][ , 6])
soc_tau_30_100cm = Re(soc_stock_ss[[7]][ , 6]) - Re(soc_stock_ss[[6]][ , 6])

soc_content_100cm_max = Re(soc_stock_ss[[5]][ , 6]) - Re(soc_stock_ss[[2]][ , 6])
soc_tau_100cm_max = Re(soc_stock_ss[[9]][ , 6]) - Re(soc_stock_ss[[6]][ , 6])

# load climate info
global_climate = env_info_ss[ , 'Koppen_Climate_2018']
var_climate = array(NA, dim = c(length(global_climate), 1))

var_climate[which(global_climate >= 1 & global_climate <= 3)] = 'A'
var_climate[which(global_climate >= 4 & global_climate <= 7)] = 'B'
var_climate[which(global_climate >= 8 & global_climate <= 16)] = 'C'
var_climate[which(global_climate >= 17 & global_climate <= 30)] = 'D'


valid_loc = which(is.na(ss_para_result[ , 1]) == 0 & stat_r2 > 0.75 & max_depth > 50
                  & is.na(var_climate) == 0
                  & soc_stock_tau_mean[ , 1] > 1000
                  & soc_stock_tau_mean[ , 1] < 1000000
                  & soc_stock_tau_mean[ , 2] > 0
                  & soc_stock_tau_mean[ , 2] < 10000)

soc_stock_tau_mean = soc_stock_tau_mean[valid_loc, ]
bulk_process_mean = bulk_process_mean[valid_loc, ]
soc_content_0_30cm = soc_content_0_30cm[valid_loc] # unit gC/m2
soc_tau_0_30cm = soc_tau_0_30cm[valid_loc]

soc_content_30_100cm = soc_content_30_100cm[valid_loc] # unit gC/m2
soc_tau_30_100cm = soc_tau_30_100cm[valid_loc]

soc_content_100cm_max = soc_content_100cm_max[valid_loc] # unit gC/m2
soc_tau_100cm_max = soc_tau_100cm_max[valid_loc]

var_climate = var_climate[valid_loc]

var_temp = env_info_ss[valid_loc, 'Annual Mean Temperature']
var_clay = apply(env_info_ss[valid_loc, c('Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm')], 1, mean, na.rm = TRUE)
var_bd = env_info_ss[valid_loc, 'Bulk_Density_0cm'] # unit kg/m3


soc_content_0_30cm = soc_content_0_30cm/0.3/var_bd # unit gC/kg
soc_content_30_100cm = soc_content_30_100cm/0.7/env_info_ss[valid_loc, 'Bulk_Density_30cm'] # unit gC/kg
soc_content_100cm_max = soc_content_100cm_max/7.4/env_info_ss[valid_loc, 'Bulk_Density_100cm'] # unit gC/kg



#################################################################################
# Plot Figures SOC content and cue (0 - 30cm)
#################################################################################
process_scale_option = c('identity', 'identity', 'log10', 'log10', 'identity', 'identity')
process_axis_label = c('Microbial CUE', 
                       'Input Allocation',
                       expression(paste('Baseline Decomposition (yr'^'-1', ')', sep = '')), 
                       expression(paste('Vertical Transport (yr'^'-1', ')', sep = '')),
                       'Environmental Modifers', 
                       expression(paste('Input Carbon (g C yr'^'-1', ')', sep = '')))

legend_limit_lower = apply(bulk_process_mean, 2, quantile, prob = 0.0001, na.rm = TRUE)
legend_limit_upper = apply(bulk_process_mean, 2, quantile, prob = 0.9999, na.rm = TRUE)


ivar = 1
for (ivar in 1:6) {
  current_data = data.frame(cbind(bulk_process_mean[ , ivar], soc_content_0_30cm, var_temp))
  colnames(current_data) = c('var', 'soc', 'mat')
  
  lm_model = lm(soc ~ var, data = current_data)
  lm_summary = summary(lm_model)
  lm_summary
  
  corr_test = cor.test(current_data$var, current_data$soc)
  corr_test$estimate
  
  
  text_data = data.frame('x_axis' = legend_limit_lower[ivar], 
                         'y_axis' = c(0.1), 
                         'equation' = paste('r = ', round(corr_test$estimate, 2), '\np < 0.001', sep = ''))
  
  p = ggplot(data = current_data) + 
    stat_bin_hex(aes(x = var, y = soc), bins = 100) +
    scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) + 
    scale_x_continuous(limits = c(legend_limit_lower[ivar], legend_limit_upper[ivar]), trans = process_scale_option[ivar]) +
    scale_y_continuous(limits = c(0.05, 600), trans = 'log10') + 
    geom_smooth(data = current_data, aes(x = var, y = soc), method = 'lm', color = 'black', fill = 'blue3', size = 2) +
    geom_text(data = text_data, aes(x = x_axis, y = y_axis, label = equation), hjust = 0, size = 10, show.legend = FALSE) + 
    # change the background to black and white
    theme_classic() + 
    # add title
    labs(title = '', x = process_axis_label[ivar], y =expression(paste('SOC (0 - 30cm g C kg'^'-1', ')', sep = ''))) + 
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
  
  eval(parse(text = paste('p', ivar, ' = p', sep = '')))
  
  
}


p1 = p1 + 
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA))
  
jpeg(paste('./Ensemble/corr_soc_mechanism_ss_0_30cm.jpeg', sep = ''), width = 30, height = 20, units = 'in', res = 300)
plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2)

dev.off()










