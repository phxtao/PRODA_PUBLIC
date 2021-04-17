clear,
clc;
%%
model_name = 'cesm2_clm5_cen_vr_v2';
time_domain = 'whole_time'; % 'whole_time', 'before_1985', 'after_1985', 'random_half'
% model_name = 'direct_nn';

cross_validation_num = 10;

var_names = ...
    {'ProfileNum', 'ProfileID', 'LayerNum', 'Lon', 'Lat', 'Date',...
    'Rmean', 'Rmax', 'Rmin',...
    'ESA_Land_Cover',...
    'ET', ...
    'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin',...
    'Veg_Cover', ...
    'BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', ...
    'Abs_Depth_to_Bedrock',...
    'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',...
    'CEC_0cm', 'CEC_30cm', 'CEC_100cm',...
    'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',...
    'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', ...
    'Depth_Bedrock_R', ...
    'Garde_Acid', ...
    'Occurrence_R_Horizon', ...
    'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', ...
    'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm', ...
    'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', ...
    'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', ...
    'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', ...
    'USDA_Suborder', ...
    'WRB_Subgroup', ...
    'Drought', ...
    'Elevation' ...
    'Max_Depth', ...
    'Koppen_Climate_2018',...
    'cesm2_npp', 'cesm2_npp_std',...
    'cesm2_gpp', 'cesm2_gpp_std',...
    'cesm2_vegc',...
    'nbedrock',...
    'R_Squared'};


var_nn_list = {...
    'Lon', 'Lat', ...
    'ESA_Land_Cover',...
    'BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19', ...
    'Abs_Depth_to_Bedrock',...
    'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',...
    'CEC_0cm', 'CEC_30cm', 'CEC_100cm',...
    'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',...
    'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', ...
    'Garde_Acid', ...
    'Occurrence_R_Horizon', ...
    'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', ...
    'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm', ...
    'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', ...
    'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', ...
    'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', ...
    'USDA_Suborder', ...
    'WRB_Subgroup', ...
    'Elevation' ...
    'Koppen_Climate_2018',...
    'cesm2_npp', 'cesm2_npp_std',...
    'cesm2_vegc',...
    'nbedrock'};


var_valid_list = zeros(length(var_names), 1);
  
data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/';
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/';

load([data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info_', model_name, '_', time_domain, '_maxmin_scaled.mat'], 'profile_env_info');


for ivar = 1:length(var_nn_list)
    var_loc = find(strcmp(var_nn_list{ivar}, var_names) == 1);
    
    env_info_middle = profile_env_info(:, var_loc);
    
    profile_env_info = profile_env_info(isnan(env_info_middle) == 0, :);
end

%% 
site_loc = 1:length(profile_env_info(:, 1)) - 1;

candidate_site_loc = site_loc;

test_loc = nan(round(length(site_loc)/cross_validation_num), cross_validation_num);
for isample = 1:cross_validation_num
    test_loc(:, isample) = sort(datasample(candidate_site_loc, round(length(site_loc)/cross_validation_num), 'Replace', false));
    
    candidate_site_loc = setdiff(candidate_site_loc, test_loc);
end
% 
% reshaped_site_loc = sort(reshape(test_loc, [round(length(site_loc)/cross_validation_num)*cross_validation_num, 1]));
% 
% unique(reshaped_site_loc' - site_loc(1:length(reshaped_site_loc)))
% 

save([data_dir_input, 'data4nn/cross_valid_test_loc_', model_name, '_', time_domain, '.mat'], 'test_loc');
