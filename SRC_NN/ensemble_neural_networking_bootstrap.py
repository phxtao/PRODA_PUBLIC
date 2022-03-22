import sys
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
random.seed(7)

bootstrap_num_start = int(sys.argv[1])
bootstrap_num_end = int(sys.argv[2])

print('bootstraping from ' + str(bootstrap_num_start) + ' to ' + str(bootstrap_num_end))

is_server = 1
model_name = 'cesm2_clm5_cen_vr_v2'
time_domain = 'whole_time' # 'whole_time', 'before_1985', 'after_1985', 'random_half_1', 'random_half_2'

for bootstrap_num in np.arange(bootstrap_num_start, (bootstrap_num_end+1)):
	print('processing bootstraping ' + str(bootstrap_num))
	nn_training_name = 'exp_pc_cesm2_23' + '_bootstrap_' + str(bootstrap_num)
	
	is_median_scaled = 0
	nn_loss = 'joint_loss' # 'mean_squared_error'
	nn_optimizer = 'adadelta'
	nn_batch_size = 32
	nn_epochs = 1200*5 # 1200*2
	early_stop_patience = 1200
	nn_layer_num = [256, 512, 512, 256]
	nn_drop_ratio = [0.0]*len(nn_layer_num) #[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	nn_l2_regularizer = [0.0]*len(nn_layer_num)
	
	nn_split_ratio = 0.1
	
	use_custom_activation = 0
	nn_activation = [None]*len(nn_layer_num)
	
	if use_custom_activation == 1:
		# define activation function
		def custom_activation(x):
			custom_activation = tf.keras.activations.relu(x, alpha = 0.1)
			return custom_activation
		
		for ilayer in range(len(nn_layer_num)):
			get_custom_objects().update({'custom_activation_'+str(ilayer): Activation(custom_activation)})
			nn_activation[ilayer] = 'custom_activation_'+str(ilayer)
	else:
		nn_activation = ['relu', 'relu', 'relu', 'relu']
	
	########################################################
	# Data Path
	########################################################
	if is_server == 1:
		data_dir_output = '/GFPS8p/cess11/taof/datahub/ensemble/output_data/'
		data_dir_input = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/'
	else:
		data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
		data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/'
	
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
	
	train_loc_bootstrap = loadmat(data_dir_input + 'data4nn/bootstrap_train_loc_' + model_name + '_' + time_domain + '.mat')
	train_loc_bootstrap = df(train_loc_bootstrap['train_loc'])
	train_loc = np.asarray(train_loc_bootstrap.iloc[:, (bootstrap_num-1)]) - 1
	
	# shuffle all the data
	site_list = list(range(len(currentdata_x)))
	  
	# for ishuffle in range(100):
	# 	print('Shuffling ' + str(ishuffle) + ' Time(s)')
	# 	random.shuffle(site_list)
	# # seperate all the data into train, test and validation set with the propotion of 8:1:1
	# train_loc = np.asarray(sorted(random.sample(site_list, round((1-nn_split_ratio)*len(site_list)))))
	# test_loc = np.setdiff1d(np.asarray(sorted(site_list)), train_loc)
	
	########################################################
	# Configuration NN
	########################################################
	
	# split into input and outputs
	train_x = currentdata_x[train_loc, :]
	train_y = currentdata_y[train_loc, :]
	 
	# test_x = currentdata_x[test_loc, :]
	# test_y = currentdata_y[test_loc, :]
	
	# validation_x = currentdata_x[validation_loc, :]
	# validation_y = currentdata_y[validation_loc, :]
	
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
	
	
	# design network
	
	model = Sequential()
	
	for ilayer in range(len(nn_layer_num)):
		if use_custom_activation == 1:
			if ilayer == 0:
				model.add(Dense(nn_layer_num[ilayer], input_dim = len(var4nn)))
				model.add(Activation(custom_activation, name = nn_activation[ilayer]))
				model.add(Dropout(nn_drop_ratio[ilayer]))
			else:
				model.add(Dense(nn_layer_num[ilayer]))
				model.add(Activation(custom_activation, name = nn_activation[ilayer]))
				model.add(Dropout(nn_drop_ratio[ilayer]))
		else:
			if ilayer == 0:
				model.add(Dense(nn_layer_num[ilayer], kernel_regularizer=l2(nn_l2_regularizer[ilayer]), input_dim = len(var4nn), activation=nn_activation[ilayer]))
				model.add(Dropout(nn_drop_ratio[ilayer]))
			else:
				model.add(Dense(nn_layer_num[ilayer], kernel_regularizer=l2(nn_l2_regularizer[ilayer]), activation=nn_activation[ilayer]))
				model.add(Dropout(nn_drop_ratio[ilayer]))
	
	model.add(Dense(len(para_names), activation='linear'))
	
	if nn_loss == 'joint_loss':
		model.compile(loss = joint_loss, optimizer = nn_optimizer, metrics = ['accuracy'])
	else:
		model.compile(loss = nn_loss, optimizer = nn_optimizer, metrics = ['accuracy'])
	
	model.summary()
	
	# early stopping
	early_stop = EarlyStopping(monitor = 'val_loss', mode = 'min', verbose = 2, patience = early_stop_patience)
	
	model_check = ModelCheckpoint(data_dir_output + 'neural_networking/trained_model_' + model_name + '_' + time_domain + '_' + nn_training_name + '.h5', monitor = 'val_loss', mode = 'min', verbose = 2, save_best_only = True)
	
	########################################################
	# NN Operation
	########################################################
	# fit network
	history = model.fit(x = train_x, y = train_y, epochs = nn_epochs, batch_size = nn_batch_size, validation_split = nn_split_ratio, verbose = 2, callbacks=[early_stop, model_check])
	# load best model
	best_model = load_model(data_dir_output + 'neural_networking/trained_model_' + model_name + '_' + time_domain + '_' + nn_training_name + '.h5', custom_objects={'joint_loss': joint_loss})
	
	# fig = plt.figure()    
	
	nn_predict = best_model.predict(train_x)
	
	plt.plot(history.history['loss'])    
	plt.plot(history.history['val_loss'])
	plt.yscale('log')
	plt.xscale('log')
	plt.legend(['train', 'validation'], loc = 'upper right')
	plt.savefig(data_dir_output + 'neural_networking/loss_' + model_name + '_' + time_domain + '_' + nn_training_name + '.pdf')
	plt.close()
	
	# reload_model = tf.keras.models.load_model('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DADL/fao_gsp/trained_model_test11.h5', custom_objects={'joint_loss': joint_loss})
	# model.save(data_dir_output + 'neural_networking/trained_model_' + model_name + '_' + time_domain + '_' + nn_training_name + '.h5')
	
	corr_para = [None]*len(para_names)
	for ipara in range(len(para_names)):
		corr_para[ipara] = np.corrcoef(train_y[:, ipara], nn_predict[:, ipara])[0, 1]
	
	
	
	print(corr_para)
	
	######################################################
	# Output
	######################################################
	nn_site_loc =  np.asarray(valid_loc)[train_loc]
	np.savetxt(data_dir_output + 'neural_networking/nn_para_result_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', nn_predict, delimiter = ',')
	np.savetxt(data_dir_output + 'neural_networking/nn_site_loc_' + model_name + '_' + time_domain + '_' + nn_training_name + '.csv', nn_site_loc, delimiter = ',')
# end for bootstrap in 

