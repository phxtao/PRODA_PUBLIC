library(R.matlab)
library(jcolors)
library(ggplot2)
library(cowplot)

library(viridis)

dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))

discrete.colors <- c('#3A3A3A', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')

# world map
world_coastline = rgdal::readOGR(dsn='/Users/phoenix/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp',layer = 'ne_110m_land')
world_coastline <- fortify(world_coastline)
Map.Using = world_coastline

#############################################################################
# Data Path
#############################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'
cross_valid_num = 10

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
# Load Projected SOC
#################################################################################
process_list =  c('A', 'K', 'Xi', 'V', 'I', 'NPP')
process_label = c('A_A', 'D_K', 'B_Xi', 'E_V', 'F_I', 'C_NPP')

process_description = c('CUE', 'Baseline Decomposition', 'Environmental Impacts', 'Vertical Transport', 'Input Allocation', 'Carbon Input')
climate_list = c('A', 'B', 'C', 'D', 'E_all')

manage_proposal = seq(from = -0.2, to = 0.2, by = 0.02)*100
depth_name = c('C_30cm', 'B_100cm', 'A_200cm', 'D_full_depth')

grid_lon_lat = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_', process_list[1], '_', model_name, '_', nn_exp_name, '_cross_valid_', as.character(1), '.mat', sep = ''))

grid_lon_lat = grid_lon_lat$Var.Decom.Grid
grid_lon_lat = grid_lon_lat[[1]][ , c(1:2)]


############################# global soc stock and turnover
resolution = 0.5
lat_seq = grid_lon_lat[ , 2]
# area of grid 
radius = 6371008.8
length_top = (2*pi*radius*cos(abs(lat_seq+resolution/2)/180*pi)/360)*resolution
length_down = (2*pi*radius*cos(abs(lat_seq-resolution/2)/180*pi)/360)*resolution
height = (pi*radius/180)*resolution
lat_grid_area = (length_top + length_down)*height/2

marginal_change_total = array(NA, dim = c(3, length(process_list), length(manage_proposal), cross_valid_num))

process_contribution_soc_map = array(NA, dim = c(nrow(grid_lon_lat), length(process_list), cross_valid_num))
process_contribution_tau_map = array(NA, dim = c(nrow(grid_lon_lat), length(process_list), cross_valid_num))
proda_soc_stock_tau = array(NA, dim = c(2, nrow(grid_lon_lat), cross_valid_num))


iprocess = 1
counter = 1
for (iprocess in 1:length(process_list)) {
  
  bulk_process = array(NA, dim = c(nrow(grid_lon_lat), length(manage_proposal), cross_valid_num))
  soc_stock_tau = array(NA, dim = c(nrow(grid_lon_lat), 2, length(manage_proposal), cross_valid_num))
  
  icross_valid = 1
  for (icross_valid in 1:cross_valid_num) {
    global_simu = readMat(paste(data_dir_output, 'world_simulation_analyses/marginal_sensitivity_proda_', as.character(process_list[iprocess]), '_', model_name, '_', nn_exp_name, '_cross_valid_', as.character(icross_valid), '.mat', sep = ''))
    global_simu = global_simu$Var.Decom.Grid
    
    soc_stock_tau[ , 1, , icross_valid] = Re(global_simu[[5]])# soc stock
    soc_stock_tau[ , 2, , icross_valid] = Re(global_simu[[9]]) # soc tau
    bulk_process[ , , icross_valid] = Re(global_simu[[10+iprocess]]) # corresponding managed process
    
    imanage = 1
    for (imanage in 1:length(manage_proposal)) {
      invalid_loc = which(soc_stock_tau[ , 1, imanage, icross_valid] < 0
                          | soc_stock_tau[ , 1, imanage, icross_valid] > 1000000
                          | soc_stock_tau[ , 2, imanage, icross_valid] < 0
                          | soc_stock_tau[ , 2, imanage, icross_valid] > 100000
                          | bulk_process[ , imanage, icross_valid] < 0)
      soc_stock_tau[invalid_loc, , , ] = NA
      bulk_process[invalid_loc, , ] = NA
      ### calculation of soc stock and tau
      #changes of soc stock
      marginal_change_total[1, iprocess, imanage, icross_valid] = sum(soc_stock_tau[ , 1, imanage, icross_valid]*lat_grid_area, na.rm = TRUE)/(10**15) 
      # changes of tau
      marginal_change_total[2, iprocess, imanage, icross_valid] = mean(soc_stock_tau[ , 2, imanage, icross_valid], na.rm = TRUE)
      # changes of bulk process
      marginal_change_total[3, iprocess, imanage, icross_valid] = mean((bulk_process[ , imanage, icross_valid] - bulk_process[ , 11, icross_valid])/bulk_process[ , 11, icross_valid], na.rm = TRUE)
    }
    
    ## grid level contribution
    if (iprocess == 1) {
      proda_soc_stock_tau[1, , icross_valid] = soc_stock_tau[ , 1, 11, icross_valid]
      proda_soc_stock_tau[2, , icross_valid] = soc_stock_tau[ , 2, 11, icross_valid]
    }
    
    coef_log_slope_soc = log(soc_stock_tau[ , 1, , icross_valid]/soc_stock_tau[ , 1, 11, icross_valid])/log(bulk_process[ , , icross_valid]/bulk_process[ , 11, icross_valid])
    coef_log_slope_soc = apply(coef_log_slope_soc, 1, mean, na.rm = TRUE)
    process_contribution_soc_map[ , iprocess, icross_valid] = soc_stock_tau[ , 1, 11, icross_valid]*(1.01**coef_log_slope_soc - 1)
    
    coef_log_slope_tau = log(soc_stock_tau[ , 2, , icross_valid]/soc_stock_tau[ , 2, 11, icross_valid])/log(bulk_process[ , , icross_valid]/bulk_process[ , 11, icross_valid])
    coef_log_slope_tau = apply(coef_log_slope_tau, 1, mean, na.rm = TRUE)
    process_contribution_tau_map[ , iprocess, icross_valid] = soc_stock_tau[ , 2, 11, icross_valid]*(1.01**coef_log_slope_tau - 1)
  }
}

####################################
process_contribution_soc_map_contribution = array(NA, dim = dim(process_contribution_soc_map))

for (icross_valid in 1:cross_valid_num) {
  process_contribution_soc_map_contribution[ , , icross_valid] = process_contribution_soc_map[ , , icross_valid]/apply(process_contribution_soc_map[ , , icross_valid], c(1), sum, na.rm = TRUE)*proda_soc_stock_tau[1, , icross_valid]
  
}

process_contribution_soc_map_contribution = apply(process_contribution_soc_map_contribution, c(1, 2), mean, na.rm = TRUE)/1000

for (iprocess in 1:length(process_list)) {
  current_data_map = data.frame(cbind(grid_lon_lat, process_contribution_soc_map_contribution[ , iprocess]))
  colnames(current_data_map) = c('lon', 'lat', 'contribution')
  
  p_map =
    ggplot() +
    geom_tile(data = current_data_map, aes(x = lon, y = lat, fill = contribution), na.rm = TRUE) +
    scale_fill_gradient2(name = 'Contribution to\nSOC Stock', low = '#005AB5', mid = 'white', high = '#DC3220', midpoint = 0, limits = c(-100, 100), na.value = NA, trans = 'identity', oob = scales::squish) +
    geom_polygon(data = Map.Using, aes(x = long, y = lat, group = group), fill = NA, color = 'black', size = 0.5) +
    # change the background to black and white
    theme_bw() +
    ylim(-56, 80) +
    # change the legend properties
    theme(legend.position = 'none') +
    # theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA)) +
    # change the size of colorbar
    guides(fill = guide_colorbar(barwidth = 2.5, barheight = 14)) + 
    theme(legend.text = element_text(size = 25), legend.title = element_text(size = 30)) +
    # add title
    labs(title = process_description[iprocess], x = '', y = '') + 
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
    # modify the font size
    theme(axis.title = element_text(size = 20)) + 
    # modify the margin
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
    theme(axis.text=element_text(size = 35))
  
  eval(parse(text = paste('p_map', iprocess, ' = p_map', sep = '')))
  
}

jpeg('./Ensemble/bulk_process_contribution_map_proda_CESM2.jpeg', width = 45, height = 15, units = 'in', res = 300)
p_map1 = p_map1 + 
  theme(legend.justification = c(0, 0), legend.position = c(0, 0), legend.background = element_rect(fill = NA)) 

plot_grid(p_map1, p_map3, p_map4, p_map2, p_map5, p_map6, nrow = 2)
dev.off()


####################################
col_list = c('contribution_soc', 'contribution_soc_std', 'contribution_tau', 'contribution_tau_std', 'process_name')
contribution_summary = data.frame(array(NA, dim = c(length(process_list), length(col_list))))
colnames(contribution_summary) = col_list

contribution_middle = array(NA, dim = c(2, length(process_list), cross_valid_num))
double_soc_tau = array(NA, dim = c(2, length(process_list), cross_valid_num))
fourper1000_soc_tau = array(NA, dim = c(2, length(process_list), cross_valid_num))

icross_valid = 1
for (icross_valid in 1:cross_valid_num) {
  iprocess = 1
  for (iprocess in 1:length(process_list)) {
    
    regression_data_soc = marginal_change_total[1, iprocess, , icross_valid]
    regression_data_tau = marginal_change_total[2, iprocess, , icross_valid]
    regression_data_process = marginal_change_total[3, iprocess, , icross_valid]
    
    lm_soc = lm(log(regression_data_soc/regression_data_soc[11]) ~ log(regression_data_process + 1))
    contribution_middle[1, iprocess, icross_valid] = regression_data_soc[11]*(1.01**lm_soc$coefficients[2] - 1)
    double_soc_tau[1, iprocess, icross_valid] = exp((log(2) - lm_soc$coefficients[1])/lm_soc$coefficients[2]) - 1
    # fourper1000_soc_tau[1, iprocess, icross_valid] = exp((log((1.004)**30) - lm_soc$coefficients[1])/lm_soc$coefficients[2]) - 1
    fourper1000_soc_tau[1, iprocess, icross_valid] = exp((log((1 + 0.1)) - lm_soc$coefficients[1])/lm_soc$coefficients[2]) - 1
    
    lm_tau = lm(log(regression_data_tau/regression_data_tau[11]) ~ log(regression_data_process + 1))
    contribution_middle[2, iprocess, icross_valid] = regression_data_tau[11]*(1.01**lm_tau$coefficients[2] - 1)
    double_soc_tau[2, iprocess, icross_valid] = exp(lm_tau$coefficients[2]*log(double_soc_tau[1, iprocess, icross_valid] + 1) + lm_tau$coefficients[1])
    fourper1000_soc_tau[2, iprocess, icross_valid] = exp(lm_tau$coefficients[2]*log(fourper1000_soc_tau[1, iprocess, icross_valid] + 1) + lm_tau$coefficients[1])
  }
  
  contribution_middle[1, , icross_valid] = contribution_middle[1, , icross_valid] # /sum(contribution_middle[1, , icross_valid])*regression_data_soc[11]
  contribution_middle[2, , icross_valid] = contribution_middle[2, , icross_valid] # /sum(contribution_middle[2, , icross_valid])*regression_data_tau[11]
}

apply(double_soc_tau, c(1, 2), mean, na.rm = TRUE)
apply(double_soc_tau, c(1, 2), sd, na.rm = TRUE)

apply(fourper1000_soc_tau, c(1, 2), mean, na.rm = TRUE)
apply(fourper1000_soc_tau, c(1, 2), sd, na.rm = TRUE)



contribution_summary$process_name = process_label
contribution_summary$contribution_soc = apply(contribution_middle, c(1, 2), mean, na.rm = TRUE)[1, ]
contribution_summary$contribution_tau = apply(contribution_middle, c(1, 2), mean, na.rm = TRUE)[2, ]
contribution_summary$contribution_soc_std = apply(contribution_middle, c(1, 2), sd, na.rm = TRUE)[1, ]
contribution_summary$contribution_tau_std = apply(contribution_middle, c(1, 2), sd, na.rm = TRUE)[2, ]


# color_scheme = c('#1d3554', '#E48F1B', '#570D32', '#7F8E39', '#42858C', '#DFE07C') #jcolors(palette = c("pal7"))
color_scheme = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00')
line_label = c('Microbial CUE', 'Environmental Modifers', 'Carbon Input', 'Baseline Decomposition', 'Vertical Transport', 'Input Allocation')

p_soc = 
  ggplot() + 
  geom_bar(data = contribution_summary, aes(x = process_name, y = contribution_soc, fill = process_name), color = 'black', stat = 'identity', position = position_dodge(0.6), size = 1, width = 0.8) + 
  geom_errorbar(data = contribution_summary, aes(x = process_name, y = contribution_soc, ymin = contribution_soc - contribution_soc_std, ymax = contribution_soc + contribution_soc_std), stat = 'identity', size = 1, width = 0.4) +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  scale_x_discrete(name = '', labels = line_label) +
  geom_hline(yintercept = 0, color = 'black', size = 1, linetype = 'solid', alpha = 1) +
  # change the background to black and white
  theme_classic() +
  theme(legend.position = 'None') +
  # theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  # theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # add title
  labs(y = 'Contribution to SOC Stock (Pg)', x = '') +
  # modify the position of title
  # modify the font sizea
  theme(axis.title = element_text(size = 25)) +
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=25)) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) 

p_tau =
  ggplot() + 
  geom_bar(data = contribution_summary, aes(x = process_name, y = contribution_tau, fill = process_name), color = 'black', stat = 'identity', position = position_dodge(0.6), size = 1, width = 0.8) + 
  geom_errorbar(data = contribution_summary, aes(x = process_name, y = contribution_tau, ymin = contribution_tau - contribution_tau_std, ymax = contribution_tau + contribution_tau_std), stat = 'identity', size = 1, width = 0.4) +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  scale_x_discrete(name = '', labels = line_label) +
  geom_hline(yintercept = 0, color = 'black', size = 1, linetype = 'solid', alpha = 1) +
  # change the background to black and white
  theme_classic() +
  theme(legend.position = 'None') +
  # theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  # theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # add title
  labs(y = 'Contribution to Turnover Time (yr)', x = '') +
  # modify the position of title
  # modify the font sizea
  theme(axis.title = element_text(size = 25)) +
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=25)) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) 

proda_total_soc = mean(marginal_change_total[1, 1, 11, ])
proda_total_soc_std = sd(marginal_change_total[1, 1, 11, ])
proda_total_tau = mean(marginal_change_total[2, 1, 11, ])
proda_total_tau_std = sd(marginal_change_total[2, 1, 11, ])


jpeg('./Ensemble/bulk_process_contribution_proda_soc_CESM2.jpeg', width = 10, height = 10, units = 'in', res = 300)
p_soc + 
  scale_x_discrete(labels = c('', '', '', '', '', '')) +
  theme(legend.justification = c(1, 1), legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) + 
  theme(axis.title = element_text(size = 35)) +
  labs(y = 'Contribution to Global SOC Stock (Pg)', x = '+1% Change of Mechanisms') + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text.x = element_text(size=10)) +
  theme(axis.text.y = element_text(size=32)) + 
  theme(axis.line = element_line(size = 1), axis.ticks.y= element_line(size = 1), axis.ticks.x = element_line(size = 0))

dev.off()

####################################
marginal_change_total_mean = apply(marginal_change_total, c(1, 2, 3), mean, na.rm = TRUE)
marginal_change_total_std = apply(marginal_change_total, c(1, 2, 3), sd, na.rm = TRUE)

col_list = c('stock', 'tau', 'stock_std', 'tau_std', 'process', 'process_std', 'process_name')
marginal_change_summary = data.frame(array(NA, dim = c(length(process_list)*length(manage_proposal), length(col_list))))
colnames(marginal_change_summary) = col_list

counter = 1

for (iprocess in 1:length(process_list)) {
  
  for (imanage in 1:length(manage_proposal)) {
    marginal_change_summary[counter, 'process_name'] = process_label[iprocess]
    marginal_change_summary[counter, 'stock'] = marginal_change_total_mean[1, iprocess, imanage]
    marginal_change_summary[counter, 'stock_std'] = marginal_change_total_std[1, iprocess, imanage]
    
    marginal_change_summary[counter, 'tau'] = marginal_change_total_mean[2, iprocess, imanage]
    marginal_change_summary[counter, 'tau_std'] = marginal_change_total_std[2, iprocess, imanage]
    
    marginal_change_summary[counter, 'process'] = marginal_change_total_mean[3, iprocess, imanage]
    marginal_change_summary[counter, 'process_std'] = marginal_change_total_std[3, iprocess, imanage]
    
    counter = counter + 1
  }
}


jpeg(paste('./Ensemble/marginal_change_soc_stock_nolegend_proda.jpeg', sep = ''), width = 10, height = 10, units = 'in', res = 300)

line_label = c('Microbial CUE', 'Environmental Modifers', 'Carbon Input', 'Baseline Decomposition', 'Vertical Transport', 'Input Allocation')

ggplot(data = marginal_change_summary) +
  geom_ribbon(aes(x = process*100, ymin = stock - stock_std, ymax = stock + stock_std, fill = process_name), alpha = 0.15) +
  geom_line(aes(x = process*100, y = stock, color = process_name), alpha = 1, size = 2) +
  geom_hline(yintercept = marginal_change_summary$stock[11], color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  # geom_hline(yintercept = marginal_change_summary$stock[11]*(1+0.004)**30, color = 'red', size = 2, linetype = 'dashed', alpha = 1) +
  geom_hline(yintercept = marginal_change_summary$stock[11]*(1+0.1), color = 'red', size = 2, linetype = 'dashed', alpha = 1) +
  geom_vline(xintercept = 0, color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  scale_x_continuous(trans = 'identity', limits = c(-10, 10), n.breaks = 7) +
  scale_y_continuous(trans = 'identity', limits = c(1500, 3500), n.breaks = 7) +
  # annotate('text', x = -10, y = marginal_change_summary$stock[11]*(1+0.004)**30, label = '"4per1000" by 2050 target', hjust = 0, vjust = -0.5, size = 10, color = 'red') +
  annotate('text', x = -10, y = marginal_change_summary$stock[11]*(1+0.1), label = '10% increase of total stock', hjust = 0, vjust = -0.5, size = 10, color = 'red') +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') + 
  theme(legend.justification = 'right', legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # add title
  labs(y = 'Global SOC Stock (Pg C)', x = paste('Proportional Changes of Mechanisms (%)')) +
  # modify the position of title
  # modify the font sizea
  theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=35)) 

dev.off()


jpeg(paste('./Ensemble/marginal_change_soc_stock_proda.jpeg', sep = ''), width = 12, height = 7, units = 'in', res = 300)


line_label = c('Microbial CUE', 'Environmental Modifers', 'Carbon Input', 'Baseline Decomposition', 'Vertical Transport', 'Input Allocation')

ggplot(data = marginal_change_summary) +
  geom_ribbon(aes(x = process*100, ymin = stock - stock_std, ymax = stock + stock_std, fill = process_name), alpha = 0.15) +
  geom_line(aes(x = process*100, y = stock, color = process_name), alpha = 1, size = 2) +
  geom_hline(yintercept = marginal_change_summary$stock[11], color = 'grey4', size = 2, linetype = 'dashed', alpha = 0.7) +
  # geom_hline(yintercept = marginal_change_summary$stock[11]*(1+0.004*30), color = 'red', size = 2, linetype = 'dashed', alpha = 0.7) +
  geom_hline(yintercept = marginal_change_summary$stock[11]*(1+0.1), color = 'red', size = 2, linetype = 'dashed', alpha = 0.7) +
  scale_x_continuous(trans = 'identity', limits = c(-10, 10), n.breaks = 7) +
  scale_y_continuous(trans = 'identity', limits = c(1500, 3500)) +
  scale_color_manual(name = '', labels = line_label, values = color_scheme) +
  scale_fill_manual(name = '', labels = line_label, values = color_scheme) +
  # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') + 
  theme(legend.justification = 'right', legend.position = 'right', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # add title
  labs(y = 'Global SOC Stock (Pg C)', x = paste('Proportional Changes of Mechanisms (%)')) +
  # modify the position of title
  # modify the font sizea
  theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.15, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=25)) 

dev.off()

