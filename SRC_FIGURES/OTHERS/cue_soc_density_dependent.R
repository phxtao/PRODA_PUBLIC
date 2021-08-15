library(ggplot2)
library(cowplot)
library(jcolors)
library(gridExtra)
library(viridis)



dev.off()
##
rm(list = ls())

setwd('/Users/phoenix/Google_Drive/R')

## Jet colorbar function
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
diff.colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))


####################################
# density dependent microbial model
####################################
## two pool microbial model
input = 0.00016
k_b = 0.00028
k_mu = 250
v_max = 0.01


cue = seq(from = 0.1, to = 0.9, by = 0.01)

beta_max = 5
beta_seq = seq(from = 1, to = beta_max, by = 0.1)

current_data_micro = c()
for (beta in beta_seq) {
  carbon_biomass = (cue*input/((1 - cue)*k_b))**(1/beta)
  carbon_soc = k_mu*k_b*carbon_biomass**(beta - 1)/(cue*v_max - k_b*carbon_biomass**(beta - 1))
  
  current_data_micro = rbind(current_data_micro, cbind(cue, (carbon_soc/carbon_soc[1]-1)*100, (carbon_biomass/carbon_biomass[1]-1)*100, beta))
  
}

current_data_micro = data.frame(current_data_micro)
colnames(current_data_micro) = c('cue', 'soc', 'biomass', 'beta')

## three pool linear model
f = 0.94
f_s = 0.31
k_uptake = 0.005
k_d = 0.001
f_b = cue
f_bs = 0.5
f_d = 0.31
k_b = 0.00028
k_s = 5.6*10**(-6)

carbon_doc = input*((1 - f) + f*f_s)/(k_uptake + k_d + k_uptake*f_b*(f_bs - 1 - f_bs*f_s) - f_d*k_d*f_s)
carbon_biomass = k_uptake*carbon_doc/k_b
carbon_soc_linear = (f*input + carbon_doc*(f_d*k_d + k_uptake*f_b*f_bs))/k_s
carbon_soc_linear = (carbon_soc_linear/carbon_soc_linear[1] - 1)*100

p1 =
  ggplot() +
  geom_line(data = current_data_micro, aes(x = cue, y = soc, group = beta, color = beta), size = 1) + 
  scale_color_gradientn(name = expression(paste(beta, sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(1, beta_max), trans = 'identity', oob = scales::squish) +
    geom_line(aes(x = cue, y = carbon_soc_linear), color = 'black', size = 1.5) + 
    geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed', size = 1.5) + 
    scale_x_continuous(limits = c(0.1, 0.9)) + 
    scale_y_continuous(limits = c(-100, 100)) + 
    # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') + 
  theme(legend.justification = c(0., 1), legend.position = c(0., 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # add title
  labs(y = 'Relative changes of SOC content \n at steady state (%)', x = paste('CUE')) +
  # modify the position of title
  # modify the font sizea
  theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=30)) 

p2 =
  ggplot() +
  geom_line(data = current_data_micro, aes(x = cue, y = biomass, group = beta, color = beta), size = 1) + 
  scale_color_gradientn(name = expression(paste(beta, sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(1, beta_max), trans = 'identity', oob = scales::squish) +
  geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed', size = 1.5) + 
  
  scale_x_continuous(limits = c(0.1, 0.9)) + 
  scale_y_continuous(limits = c(0, 600)) + 
  
  # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') + 
  theme(legend.justification = c(0.9, 1), legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # add title
  labs(y = 'Relative changes of biomass carbon \n at steady state (%)', x = paste('CUE')) +
  # modify the position of title
  # modify the font sizea
  theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=30)) 


jpeg(paste('./Ensemble/cue_soc_density_dependent.jpeg', sep = ''), width = 20, height = 10, units = 'in', res = 300)
plot_grid(p1, p2, nrow = 1)
dev.off()


####################################
# logistic model
####################################
population_0 = 0.001
growth_rate = 0.005
carrying_capacity = 3.5

time_series = 4000

beta_max = 5

population_record = data.frame(array(NA, dim = c((time_series+1)*beta_max, 3)))
colnames(population_record) = c('time', 'population', 'beta')

population_record[(time_series+1)*0+1, ] = c(0, population_0, 1)
population_record[(time_series+1)*1+1, ] = c(0, population_0, 2)
population_record[(time_series+1)*2+1, ] = c(0, population_0, 3)
population_record[(time_series+1)*3+1, ] = c(0, population_0, 4)
population_record[(time_series+1)*4+1, ] = c(0, population_0, 5)

for (time in 1:time_series) {
  for (beta in 1:beta_max) {
    population_record[(time_series+1)*(beta - 1)+time+1, 1] = time+1
    current_polulation = population_record[time, 2]
    population_dynamic = current_polulation*growth_rate*(1 - current_polulation**(beta)/carrying_capacity)
    if (current_polulation + population_dynamic < 0) { 
      population_record[(time_series+1)*(beta - 1)+time+1, 2] = 0
    } else { 
      population_record[(time_series+1)*(beta - 1)+time+1, 2] = current_polulation + population_dynamic
      
    }
    population_record[(time_series+1)*(beta - 1)+time+1, 3] = beta
    
  }
  
}


current_data = population_record

jpeg(paste('./Ensemble/logistic_growth_density_dependent.jpeg', sep = ''), width = 15, height = 10, units = 'in', res = 300)
ggplot() + 
  geom_line(data = current_data, aes(x = time, y = population, group = beta, color = beta), size = 2) + 
  scale_color_gradientn(name = expression(paste(beta, sep = '')), colours = rev(viridis(15)), na.value="transparent", limits = c(1, beta_max), trans = 'identity', oob = scales::squish) +
  # scale_y_continuous(limits = c(0, 100)) + 
  # change the background to black and white
  theme_classic() +
  # theme(legend.position = 'None') + 
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # add title
  labs(y = 'Population', x = paste('Time')) +
  # modify the position of title
  # modify the font sizea
  theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=30)) 
dev.off()


