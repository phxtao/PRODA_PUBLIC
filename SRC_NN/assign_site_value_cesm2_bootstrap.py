from pandas import Series as se
from pandas import DataFrame as df
from scipy.io import loadmat
import scipy.stats
import pandas as pd
import numpy as np

import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Dropout, Flatten, MaxPooling2D
from tensorflow.keras.layers import Dropout
from tensorflow.keras.layers import Activation
from tensorflow.keras.optimizers import SGD
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.regularizers import l2
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import load_model

from matplotlib import pyplot as plt

import random

########################################################
# Environment Setting
########################################################
is_server = 1
model_name = 'cesm2_clm5_cen_vr_v2'
time_domain = 'whole_time' # 'whole_time', 'before_1985', 'after_1985', 'random_half_1', 'random_half_2'

is_median_scaled = 0
########################################################
# Data Path
########################################################
if is_server == 1:
	data_dir_output = '/GFPS8p/cess11/taof/datahub/ensemble/output_data/'
	data_dir_input = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/'
else:
	data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
	data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'


for bootstrap_num in np.arange(1, 201):
	
	print('processing bootstraping number ' + str(bootstrap_num))
	nn_training_name = 'exp_pc_cesm2_23' + '_bootstrap_' + str(bootstrap_num)
	# nn_training_name = 'exp_pc_cesm2_23' + '_cross_valid_0_' + str(bootstrap_num)
	########################################################
	# Import para info after bayesian method
	########################################################
	# laod para info after MCMC
	para_without_trans = loadmat(data_dir_output + 'mcmc_summary_' + model_name + '/' + model_name + '_' + time_domain + '_para_mean_cleansed.mat')
	para_without_trans = df(para_without_trans['para_mean'])
	
	if model_name[0:5] == 'cesm2':
		para_names = ['diffus', 'cryo', 'q10', 'efolding', 'taucwd', 'taul1', 'taul2', 'tau4s1', 'tau4s2', 'tau4s3', 'fl1s1', 'fl2s1', 'fl3s2', 'fs1s2', 'fs1s3', 'fs2s1', 'fs2s3', 'fs3s1', 'fcwdl2', 'w-scaling', 'beta']
	
	if model_name == 'clm_cen_vr_v3':
		para_names = ['diffus', 'cryo', 'q10', 'efolding', 'tau4s1', 'tau4s2', 'tau4s3', 'fl1s1', 'fl2s1', 'fl3s2', 'fs1s2', 'fs1s3', 'fs2s1', 'fs2s3', 'fs3s1', 'fcwdl2', 'beta', 'maxpsi']
	
	if model_name == 'clm_bgc_vr':
		para_names = ['diffus', 'cryo', 'q10', 'efolding', 'tau4s1', 'tau4s2', 'tau4s3', 'tau4s4', 'fl1s1', 'fl2s2', 'fl3s3', 'fs1s2', 'fs2s3', 'fs3s4', 'fcwdl2', 'beta', 'maxpsi']
	
	para_without_trans.columns = para_names
	
	# environmental info of soil profiles
	if is_median_scaled == 1:
		env_info = loadmat(data_dir_input + 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info_' + model_name + '_' + time_domain + '_median_scaled.mat')
	else:
		env_info = loadmat(data_dir_input + 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info_' + model_name  + '_' + time_domain + '_maxmin_scaled.mat')
	
	env_info = df(env_info['profile_env_info'])
	
	env_info_names = ['ProfileNum', 'ProfileID', 'LayerNum', 'Lon', 'Lat', 'Date', \
	'Rmean', 'Rmax', 'Rmin', \
	'ESA_Land_Cover', \
	'ET', \
	'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin', \
	'Veg_Cover', \
	'BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', \
	'Abs_Depth_to_Bedrock', \
	'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',\
	'CEC_0cm', 'CEC_30cm', 'CEC_100cm', \
	'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm', \
	'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', \
	'Depth_Bedrock_R', \
	'Garde_Acid', \
	'Occurrence_R_Horizon', \
	'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', \
	'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm', \
	'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', \
	'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', \
	'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', \
	'USDA_Suborder', \
	'WRB_Subgroup', \
	'Drought', \
	'Elevation', \
	'Max_Depth', \
	'Koppen_Climate_2018', \
	'cesm2_npp', 'cesm2_npp_std', \
	'cesm2_gpp', 'cesm2_gpp_std', \
	'cesm2_vegc', \
	'nbedrock', \
	'R_Squared']
	
	env_info.columns = env_info_names
	
	
	# variables used in training the NN
	var4nn = ['Lon', 'Lat', \
	'ESA_Land_Cover', \
	# 'IGBP', \
	# 'Climate', \
	# 'Soil_Type', \
	# 'NPPmean', 'NPPmax', 'NPPmin', \
	# 'Veg_Cover', \
	'BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', \
	'Abs_Depth_to_Bedrock', \
	'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',\
	'CEC_0cm', 'CEC_30cm', 'CEC_100cm', \
	'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm', \
	'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', \
	# 'Depth_Bedrock_R', \
	'Garde_Acid', \
	'Occurrence_R_Horizon', \
	'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', \
	'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm', \
	'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', \
	'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', \
	'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', \
	'USDA_Suborder', \
	'WRB_Subgroup', \
	# 'Drought', \
	'Elevation', \
	# 'Max_Depth', \
	'Koppen_Climate_2018', \
	'cesm2_npp', 'cesm2_npp_std', \
	# 'cesm2_gpp', 'cesm2_gpp_std', \
	'cesm2_vegc', \
	'nbedrock']
	
	# R2 of soil profiles in MCMC to eliminate some profiles
	r_squared = np.array(env_info.loc[:, 'R_Squared'])
	
	########################################################
	# specify the profiles 
	########################################################
	# delete nan value and R2 < 0
	valid_loc = env_info.loc[:, 'ProfileNum'] 
	currentdata_y = np.array(para_without_trans)
	currentdata_x = np.array(env_info.loc[:, var4nn])
	
	# eliminate nan value
	for ivar in range(len(var4nn)):
		r_squared = r_squared[(np.isnan(currentdata_x[:, ivar]) == False)]
		currentdata_y = currentdata_y[(np.isnan(currentdata_x[:, ivar]) == False), :]
		valid_loc = valid_loc[(np.isnan(currentdata_x[:, ivar]) == False)]
		currentdata_x = currentdata_x[(np.isnan(currentdata_x[:, ivar]) == False), :]
	
	train_loc = np.arange(0, len(currentdata_x[:, 0]))
	
	########################################################
	# Configuration NN
	########################################################
	
	# split into input and outputs
	train_x = currentdata_x[train_loc, :]
	train_y = currentdata_y[train_loc, :]
	
	print(train_x.shape, train_y.shape)
	
	# define a joint loss
	para_mse = 1000*0.5
	para_ratio = 50*0.5
	def joint_loss (y_true, y_pred):
		# mse
		mse_loss = K.mean(K.square(y_true - y_pred))
		# mean absolute ratio  error
		ratio_loss = K.mean(K.abs((y_true - y_pred)/y_true))
		# return the joint loss
		return para_mse*mse_loss * para_ratio*ratio_loss
	
	def ratio_loss (y_true, y_pred):
		# mean absolute ratio  error
		ratio_loss = K.mean(K.abs((y_true - y_pred)/y_true))
		return ratio_loss
		
	
	# load network
	model = tf.keras.models.load_model(data_dir_output + 'neural_networking/trained_model_' + model_name + '_' + time_domain + '_' + nn_training_name + '.h5', custom_objects={'joint_loss': joint_loss})
	
	# predict parameter values of each site
	site_predict = model.predict(train_x)
	
	######################################################
	# Output
	######################################################
	nn_site_loc =  np.asarray(valid_loc)[train_loc]
	np.savetxt(data_dir_output + 'neural_networking/nn_para_result_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', site_predict, delimiter = ',')
	np.savetxt(data_dir_output + 'neural_networking/nn_site_loc_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', nn_site_loc, delimiter = ',')
	
# end for bootstrap_num
