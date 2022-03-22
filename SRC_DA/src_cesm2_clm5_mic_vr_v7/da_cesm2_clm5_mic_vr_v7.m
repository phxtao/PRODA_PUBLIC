% clear;
% clc;
%
% cesm2_case_name = 'sasu_f05_g16_checked_step4';
% start_year = 661;
% end_year = 680;
%
% model_name = 'cesm2_clm5_mic_vr_v1';
%
% start_id = 1; %36; 333; 17638; %17427; %17600;
% end_id = 5000; %36; 333; 17638; %17427; % 17600;
% is_resubmit = 0;
%
% workers = 3;
%
%============================above is only for pc use============================================

% disp(['Profiles start from ', num2str(start_id), ' to ', num2str(end_id)]);

% parallel setting
delete(gcp('nocreate'));
parpool(workers);

%
warning('off');
format long e;


%% paths
% mac
% data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/';
% cd(['/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/ENSEMBLE/SRC_DA/src_', model_name, '/']);

% server
data_path = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/';
cd(['/GFPS8p/cess11/taof/ensemble/src_da/src_', model_name, '/']);
%
%% set vertical soil pools
month_num = 12;
soil_cpool_num = 7;
soil_decom_num = 20;

%% load wosis data
env_info = load([data_path, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat'], 'EnvInfo'); % NPP in this file will be used
env_info = env_info.EnvInfo;
% layer_info: "profile_id, date, upper_depth, lower_depth, node_depth, soc_layer_weight, soc_stock, bulk_denstiy, is_pedo"
wosis_profile_info = ncread([data_path, 'wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'], 'soc_profile_info'); % wosis profile info
wosis_soc_info = ncread([data_path, 'wosis_2019_snap_shot/soc_data_integrate_wosis_2019_snapshot_hugelius_mishra.nc'], 'data_soc_integrate'); % wosis SOC info

%% soil depth information
% width between two interfaces
dz = [2.000000000000000E-002, 4.000000000000000E-002, 6.000000000000000E-002, ...
    8.000000000000000E-002, 0.120000000000000, 0.160000000000000, ...
    0.200000000000000, 0.240000000000000, 0.280000000000000, ...
    0.320000000000000, 0.360000000000000, 0.400000000000000, ...
    0.440000000000000, 0.540000000000000, 0.640000000000000, ...
    0.740000000000000, 0.840000000000000, 0.940000000000000, ...
    1.04000000000000, 1.14000000000000, 2.39000000000000, ...
    4.67553390593274, 7.63519052838329, 11.1400000000000, ...
    15.1154248593737]';

% depth of the interface
zisoi = [2.000000000000000E-002, 6.000000000000000E-002, ...
    0.120000000000000, 0.200000000000000, 0.320000000000000, ...
    0.480000000000000, 0.680000000000000, 0.920000000000000, ...
    1.20000000000000, 1.52000000000000, 1.88000000000000, ...
    2.28000000000000, 2.72000000000000, 3.26000000000000, ...
    3.90000000000000, 4.64000000000000, 5.48000000000000, ...
    6.42000000000000, 7.46000000000000, 8.60000000000000, ...
    10.9900000000000, 15.6655339059327, 23.3007244343160, ...
    34.4407244343160, 49.5561492936897]';

% depth of the node
zsoi = [1.000000000000000E-002, 4.000000000000000E-002, 9.000000000000000E-002, ...
    0.160000000000000, 0.260000000000000, 0.400000000000000, ...
    0.580000000000000, 0.800000000000000, 1.06000000000000, ...
    1.36000000000000, 1.70000000000000, 2.08000000000000, ...
    2.50000000000000, 2.99000000000000, 3.58000000000000, ...
    4.27000000000000, 5.06000000000000, 5.95000000000000, ...
    6.94000000000000, 8.03000000000000, 9.79500000000000, ...
    13.3277669529664, 19.4831291701244, 28.8707244343160, ...
    41.9984368640029]';

% depth between two node
dz_node = zsoi - [0; zsoi(1:end-1)];

%% input from cesm2
% cesm2 resolution
cesm2_resolution_lat = 180/384;
cesm2_resolution_lon = 360/576;
lon_grid = (-180 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 180 - cesm2_resolution_lon/2)';
lat_grid = (90 - cesm2_resolution_lat/2 : -cesm2_resolution_lat : -90 + cesm2_resolution_lat/2)';

% load cesm2 input
var_name_list = {'nbedrock', 'ALTMAX', 'ALTMAX_LASTYEAR', 'CELLSAND', 'NPP',...
    'SOILPSI', 'TSOI', ...
    'W_SCALAR', 'T_SCALAR', 'O_SCALAR', 'FPI_vr', ...
    'LITR1_INPUT_ACC_VECTOR', 'LITR2_INPUT_ACC_VECTOR', 'LITR3_INPUT_ACC_VECTOR', 'CWD_INPUT_ACC_VECTOR', ...
    'TOTSOMC'};

var_name_list_rename =  {'cesm2_simu_nbedrock', 'cesm2_simu_altmax', 'cesm2_simu_altmax_last_year', 'cesm2_simu_cellsand', 'cesm2_simu_npp',...
    'cesm2_simu_soil_water_potnetial', 'cesm2_simu_soil_temperature', ...
    'cesm2_simu_w_scalar', 'cesm2_simu_t_scalar', 'cesm2_simu_o_scalar', 'cesm2_simu_n_scalar', ...
    'cesm2_simu_input_vector_litter1', 'cesm2_simu_input_vector_litter2', 'cesm2_simu_input_vector_litter3', 'cesm2_simu_input_vector_cwd', ...
    'cesm2_simu_soc_stock'};

for ivar = 1:length(var_name_list)
    %% load simulation from CESM2
    load([data_path, 'cesm2_simu/spinup_ss/',...
        cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
        '_', var_name_list{ivar}, '.mat'], 'var_record_monthly_mean');
    eval([var_name_list_rename{ivar} ' = var_record_monthly_mean;']);
    
end

for ilayer = 1:soil_decom_num
    cesm2_simu_input_vector_litter1(:, :, ilayer, :) = cesm2_simu_input_vector_litter1(:, :, ilayer, :).*dz(ilayer);
    cesm2_simu_input_vector_litter2(:, :, ilayer, :) = cesm2_simu_input_vector_litter2(:, :, ilayer, :).*dz(ilayer);
    cesm2_simu_input_vector_litter3(:, :, ilayer, :) = cesm2_simu_input_vector_litter3(:, :, ilayer, :).*dz(ilayer);
    cesm2_simu_input_vector_cwd(:, :, ilayer, :) = cesm2_simu_input_vector_cwd(:, :, ilayer, :).*dz(ilayer);
end

cesm2_simu_input_sum_litter1 = reshape(sum(cesm2_simu_input_vector_litter1, 3), [384, 576, 12]);
cesm2_simu_input_sum_litter2 = reshape(sum(cesm2_simu_input_vector_litter2, 3), [384, 576, 12]);
cesm2_simu_input_sum_litter3 = reshape(sum(cesm2_simu_input_vector_litter3, 3), [384, 576, 12]);
cesm2_simu_input_sum_cwd = reshape(sum(cesm2_simu_input_vector_cwd, 3), [384, 576, 12]);

clearvars cesm2_simu_input_vector_litter1 cesm2_simu_input_vector_litter2 cesm2_simu_input_vector_litter3 cesm2_simu_input_vector_cwd

%% Shuffled Complex Evolution Method

% cleansed_profile_id = load([data_path, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info_cesm2_clm5_cen_vr_v2_whole_time_maxmin_scaled.mat']);
% cleansed_profile_id = cleansed_profile_id.profile_env_info(:, 1);
% 
% if is_resubmit == 0
%     profile_collection = cleansed_profile_id(cleansed_profile_id >= start_id & cleansed_profile_id <= end_id);
%     for iperm = 1:10
%         profile_collection = profile_collection(randperm(length(profile_collection)));
%     end
%     % profile_collection = (start_id:end_id);
% else
%     profile_collection = fun_resubmit_profiles(start_id, end_id, profile_num, model_name);
%     disp(['Resubmitting profile number: ', num2str(length(profile_collection))]);
% end

sample_profile_id = load([data_path, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_representative_profiles.mat']);
sample_profile_id = sample_profile_id.sample_profile_id;

profile_collection = reshape(sample_profile_id((2*(itask-1)+1):(2*(itask-1)+2), :), [100, 1]);

profile_range = (1:length(profile_collection));

% parfor iprofile_hat = profile_range
for iprofile_hat = profile_range
    iprofile = profile_collection(iprofile_hat);
    
    profile_id = wosis_profile_info(iprofile, 1);
    
    file_directory = [data_path, '../output_data/',  model_name, '/', model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.mat'];
    if exist(file_directory, 'file') == 0
        disp([datestr(now), ' Profile ', num2str(iprofile), ' started']);
    else
        da_result = load(file_directory);
        eval(['da_result = da_result.', model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ';'])
        
        para_mean_std_hist = mean(da_result.candidate_para_std_hist, 1, 'omitnan');
        shuffle_num = length(find(isnan(para_mean_std_hist) == 0)) + 1;
        
        if da_result.r2_opt >= 0.9 || para_mean_std_hist(shuffle_num) <= 0.005 || shuffle_num == 501
            disp([datestr(now), ' Skipped profile ', num2str(iprofile), ' due to satisfied cateria']);
            continue;
        else
            disp([datestr(now), ' Profile ', num2str(iprofile), ' started']);
        end
        
    end
    
    warning('off')
    
    try
        % find currently using profile
        loc_profile = find(wosis_soc_info(:, 1) == profile_id);
        % info of the node depth of profile, and change unit from cm to m
        wosis_layer_depth = wosis_soc_info(loc_profile, 5)/100;
        % observed C info (gC/m3)
        wosis_layer_obs = wosis_soc_info(loc_profile, 7);
        % specify the number of layers in studied profiles
        num_layers = length(wosis_layer_obs);
        
        % find the lon and lat info of soil profile
        lon_profile = wosis_profile_info(iprofile, 4);
        lat_profile = wosis_profile_info(iprofile, 5);
        lat_loc = find(abs(lat_profile - lat_grid) == min(abs(lat_profile - lat_grid)));
        lon_loc = find(abs(lon_profile - lon_grid) == min(abs(lon_profile - lon_grid)));
        
        if length(lon_loc) > 1
            lon_loc = lon_loc(1);
        end
        
        if length(lat_loc) > 1
            lat_loc = lat_loc(1);
        end
        
        % input vector
        % input_vector_cwd = reshape(cesm2_simu_input_vector_cwd(lat_loc, lon_loc, :, :), [soil_decom_num, month_num]);
        % input_vector_litter1 = reshape(cesm2_simu_input_vector_litter1(lat_loc, lon_loc, :, :), [soil_decom_num, month_num]);
        % input_vector_litter2 = reshape(cesm2_simu_input_vector_litter2(lat_loc, lon_loc, :, :), [soil_decom_num, month_num]);
        % input_vector_litter3 = reshape(cesm2_simu_input_vector_litter3(lat_loc, lon_loc, :, :), [soil_decom_num, month_num]);
        input_vector_cwd = reshape(cesm2_simu_input_sum_cwd(lat_loc, lon_loc, :), [1, month_num]);
        input_vector_litter1 = reshape(cesm2_simu_input_sum_litter1(lat_loc, lon_loc, :), [1, month_num]);
        input_vector_litter2 = reshape(cesm2_simu_input_sum_litter2(lat_loc, lon_loc, :), [1, month_num]);
        input_vector_litter3 = reshape(cesm2_simu_input_sum_litter3(lat_loc, lon_loc, :), [1, month_num]);
        
        % no input information
        if isnan(input_vector_cwd(1)) == 1 || isnan(input_vector_litter1(1)) == 1 ...
                || isnan(input_vector_litter2(1)) == 1 || isnan(input_vector_litter3(1)) == 1
            disp(['No input info ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        % npp from CESM2 simulatoin
        npp_mean = reshape(cesm2_simu_npp(lat_loc, lon_loc, :), [month_num, 1]);
        % NPP info of studied profile from soc_modIS NPP mean (year 2001-2016)
        % npp_mean = env_info(iprofile, 15);
        
        % altmax current and last year
        altmax_current_profile = reshape(cesm2_simu_altmax(lat_loc, lon_loc, :), [month_num, 1]);
        altmax_lastyear_profile = reshape(cesm2_simu_altmax_last_year(lat_loc, lon_loc, :), [month_num, 1]);
        
        % nbedrock
        nbedrock = reshape(cesm2_simu_nbedrock(lat_loc, lon_loc, :), [month_num, 1]);
        
        % oxygen scalar
        xio = reshape(cesm2_simu_o_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        % nitrogen
        xin = reshape(cesm2_simu_n_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        
        % sand content
        sand = reshape(cesm2_simu_cellsand(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        sand_vector = sand;
        
        % soil temperature and water potential
        soil_temp_profile = reshape(cesm2_simu_soil_temperature(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        % soil_water_profile = reshape(cesm2_simu_soil_water_potnetial(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        % soil_temp_profile = reshape(cesm2_simu_t_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        soil_water_profile = reshape(cesm2_simu_w_scalar(lat_loc, lon_loc, 1:soil_decom_num, :), [soil_decom_num, month_num, 1]);
        
        % deal with zero value in mcmc
        % if zero, cost function will always be inf
        if length(find(wosis_layer_obs == 0)) == length(wosis_layer_obs)
            disp(['All zero value ', ' Profile ', num2str(iprofile)]);
            continue
        elseif isempty(find(wosis_layer_obs == 0, 1)) == 0
            wosis_layer_depth = wosis_layer_depth(wosis_layer_obs > 0);
            wosis_layer_obs = wosis_layer_obs(wosis_layer_obs > 0);
        end
        
        % deal with nan values in obs or depth
        layer_valid_loc = find(isnan(wosis_layer_depth) == 0 & isnan(wosis_layer_obs) == 0);
        
        if isempty(layer_valid_loc) == 0
            wosis_layer_depth = wosis_layer_depth(layer_valid_loc);
            wosis_layer_obs = wosis_layer_obs(layer_valid_loc);
        else
            disp(['No valid obsrvations in ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        % eliminate profiles having single layer
        if length(wosis_layer_obs) == 1
            disp(['Single layer ', ' Profile ', num2str(iprofile)]);
            continue
        end
        
        % layer_weighting at different layers of soil profile
        
        % layer_weight = exp(-wosis_layer_depth);
        % layer_weight([1, end], 1) = 10;
        [~, soc_stock_rank] = sort(abs(wosis_layer_obs - mean(wosis_layer_obs)), 'descend');
        [~, soc_stock_rank] = sort(soc_stock_rank);
        [~, soc_depth_rank] = sort(abs(wosis_layer_depth - mean(wosis_layer_depth)), 'descend');
        [~, soc_depth_rank] = sort(soc_depth_rank);
        soc_rank = (soc_stock_rank + soc_depth_rank)/2;
        
        layer_weight = exp(-(soc_rank - 1))*10;
        layer_weight(layer_weight < 1) = 1;
        
        % Parameter names and their initial values in MCMC
        para_name = {'bio', 'cryo', 'q10', 'efolding', 'w_scaling', ...
            'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', ...
            'mm_const_assim', 'mm_const_decom', ...
            'fcwdl2', 'fl1s1', 'fl2s1', 'fl3s4', ...
            'mic_cue', 'pdeath2soc', ...
            'beta', ...
            'allo_slope_mic'};
        
        % number of paras
        npara = length(para_name);
        para_min = zeros(npara, 1);
        para_max = ones(npara, 1);
        % set initial values of parameters
        para0 = 0.3*ones(npara, 1);
        
        % scatter(wosis_layer_depth, wosis_layer_obs, 100, '*')
        % close
        %---------------------------------------------------------------
        % SCE main
        %---------------------------------------------------------------
        warning off
        % STEP 0: hyper parameters
        %==========================================================================
        shuffle_num = 500;
        complex_num = workers; % complex number
        point_num = 2*npara + 1; % points number in each complex
        sample_size = complex_num*point_num;
        % hyper para in cce
        parents_num = npara + 1;
        evolve_num = 1;
        generation_num = 2*npara + 1;
        %==========================================================================
        
        candidate_para_value_hist = nan(complex_num, point_num, npara, (shuffle_num+1));
        candidate_cost_value_hist = nan(complex_num, point_num, (shuffle_num+1));
        candidate_steady_state_hist = nan(complex_num, point_num, (shuffle_num+1));
        candidate_mod_soc_hist = nan(complex_num, point_num, length(wosis_layer_obs), (shuffle_num+1));
        candidate_r2_hist = nan(complex_num, point_num, (shuffle_num+1));
        candidate_para_std_hist = nan(npara, (shuffle_num+1));
        
        % STEP 1: generate sample (uniform distribution as prior)
        candidate_para_value = rand(complex_num, point_num, npara);
        candidate_mod_soc = nan(complex_num, point_num, length(wosis_layer_obs));
        candidate_soc_stock = nan(complex_num, point_num, 5);
        candidate_soc_layer = nan(complex_num, point_num, soil_decom_num);
        candidate_r2 = nan(complex_num, point_num);
        candidate_cost_value = nan(complex_num, point_num);
        candidate_steady_state = nan(complex_num, point_num);
        candidate_index = reshape(1:sample_size, [complex_num, point_num]);
        % calculate soc at steady state and corresponding cost function values
        parfor icomplex = 1:complex_num
            for ipoint = 1 : point_num
                para = reshape(candidate_para_value(icomplex, ipoint, :), [npara, 1]);
                
                [~, soc_stock_summary, soc_mod, ss_index] = ...
                    matrix_fun_seasonal(para, wosis_layer_obs, nbedrock, sand_vector, npp_mean, ...
                    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
                    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
                
                if isnan(sum(soc_mod(:, 5))) == 0
                    optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod(:, 5), wosis_layer_depth, 'pchip');
                else
                    optimize_profile_soc = -10^7*ones(length(wosis_layer_depth), 1);
                end
                
                candidate_cost_value(icomplex, ipoint) = cost_fun(layer_weight, optimize_profile_soc, wosis_layer_obs);
                candidate_steady_state(icomplex, ipoint) = ss_index;
                
                candidate_mod_soc(icomplex, ipoint, :) = optimize_profile_soc;
                candidate_soc_stock(icomplex, ipoint, :) = soc_stock_summary;
                candidate_soc_layer(icomplex, ipoint, :) = soc_mod(:, 5);
                
                mod_r2 = fun_r2(wosis_layer_obs, optimize_profile_soc);
                candidate_r2(icomplex, ipoint) = mod_r2;
                
                % disp([datestr(now), ' processing complex ', num2str(icomplex), ' point ', num2str(ipoint), ' steady state ', num2str(ss_index)]);
                
            end
        end
        
        max_ss_cost = max(candidate_cost_value(candidate_steady_state == 1));
        candidate_cost_value(candidate_steady_state == 0) = candidate_cost_value(candidate_steady_state == 0) + max_ss_cost;
        
        candidate_para_value_hist(:, :, :, 1) = candidate_para_value;
        candidate_cost_value_hist(:, :, 1) = candidate_cost_value;
        candidate_steady_state_hist(:, :, 1) = candidate_steady_state;
        candidate_mod_soc_hist(:, :, :, 1) = candidate_mod_soc;
        candidate_r2_hist(:, :, 1) = candidate_r2;
        
        % STEP 2 & 3: rank points and partition into complexes
        for ishuffle = 1 : shuffle_num
            disp([datestr(now), ' processing shuffle ', num2str(ishuffle)]);
            % rank cost and det the order index
            [ranked_cost_value_vector, order_index] = sort(reshape(candidate_cost_value, [sample_size, 1]), 'ascend');
            ranked_candidate_cost_value = reshape(ranked_cost_value_vector, [complex_num, point_num]);
            % rank candidate index
            candidate_index_vector = reshape(candidate_index, [sample_size, 1]);
            ranked_candidate_index = reshape(candidate_index_vector(order_index), [complex_num, point_num]);
            % rank steady state
            candidate_steady_state_vector = reshape(candidate_steady_state, [sample_size, 1]);
            ranked_candidate_steady_state = reshape(candidate_steady_state_vector(order_index), [complex_num, point_num]);
            % rank r2
            candidate_r2_vector = reshape(candidate_r2, [sample_size, 1]);
            ranked_candidate_r2 = reshape(candidate_r2_vector(order_index), [complex_num, point_num]);
            % rank para values of candidate
            ranked_candidate_para_value = nan(complex_num, point_num, npara);
            % rank simued soc corresponding to observations
            ranked_candidate_mod_soc = nan(complex_num, point_num, length(wosis_layer_obs));
            % rank simued soc stock and layer
            ranked_candidate_soc_stock = nan(complex_num, point_num, 5);
            ranked_candidate_soc_layer = nan(complex_num, point_num, soil_decom_num);
            
            for icomplex = 1 : complex_num
                for ipoint = 1 : point_num
                    [complex_loc, point_loc] = find(candidate_index == ranked_candidate_index(icomplex, ipoint));
                    ranked_candidate_para_value(icomplex, ipoint, :) = candidate_para_value(complex_loc, point_loc, :);
                    ranked_candidate_mod_soc(icomplex, ipoint, :) = candidate_mod_soc(complex_loc, point_loc, :);
                    
                    ranked_candidate_soc_stock(icomplex, ipoint, :) = candidate_soc_stock(complex_loc, point_loc, :);
                    ranked_candidate_soc_layer(icomplex, ipoint, :) = candidate_soc_layer(complex_loc, point_loc, :);
                end
            end
            
            
            % STEP 4: envolve complexes (competitive complex evolution)
            max_ss_cost_record = max_ss_cost*ones(complex_num, 1);
            
            for igeneration = 1:generation_num
                for ievolve = 1 : evolve_num
                    ranked_candidate_cost_value_middle = ranked_candidate_cost_value;
                    ranked_candidate_steady_state_middle = ranked_candidate_steady_state;
                    ranked_candidate_para_value_middle = ranked_candidate_para_value;
                    ranked_candidate_mod_soc_middle = ranked_candidate_mod_soc;
                    ranked_candidate_soc_stock_middle = ranked_candidate_soc_stock;
                    ranked_candidate_soc_layer_middle = ranked_candidate_soc_layer;
                    ranked_candidate_r2_middle = ranked_candidate_r2;
                    
                    max_ss_cost_record_middle = max_ss_cost_record;
                    parfor icomplex = 1 : complex_num
                        
                        [evolved_ranked_candidate_cost_value, evolved_ranked_candidate_para_value, evolved_ranked_candidate_steady_state, ...
                            evolved_ranked_candidate_mod_soc, evolved_ranked_candidate_soc_stock, evolved_ranked_candidate_soc_layer, ...
                            evolved_ranked_candidate_r2, max_ss_cost_record_update] = ...
                            fun_sce_cce(icomplex, npara, complex_num, point_num, parents_num, max_ss_cost_record_middle, ... % cce variablee
                            ranked_candidate_cost_value_middle, ranked_candidate_index, ranked_candidate_para_value_middle, ranked_candidate_steady_state_middle, ... % cce variablee
                            ranked_candidate_mod_soc_middle, ranked_candidate_soc_stock_middle, ranked_candidate_soc_layer_middle, ranked_candidate_r2_middle, ... % cce variablee
                            nbedrock, sand_vector, npp_mean, ... % mic model variable
                            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ... % mic model variable
                            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin, ... % mic model variable
                            zsoi, soil_decom_num, wosis_layer_depth, wosis_layer_obs, layer_weight); % other variable
                        
                        % disp([datestr(now), ' shuffle ', num2str(ishuffle), ' generation ', ...
                        %     num2str(igeneration), ' evolving ', num2str(ievolve), ' complex ', num2str(icomplex)])
                        
                        ranked_candidate_cost_value(icomplex, :) = evolved_ranked_candidate_cost_value(icomplex, :);
                        ranked_candidate_para_value(icomplex, :, :) = evolved_ranked_candidate_para_value(icomplex, :, :);
                        ranked_candidate_steady_state(icomplex, :) = evolved_ranked_candidate_steady_state(icomplex, :);
                        ranked_candidate_r2(icomplex, :) = evolved_ranked_candidate_r2(icomplex, :);
                        ranked_candidate_mod_soc(icomplex, :, :) = evolved_ranked_candidate_mod_soc(icomplex, :, :);
                        ranked_candidate_soc_stock(icomplex, :, :) = evolved_ranked_candidate_soc_stock(icomplex, :, :);
                        ranked_candidate_soc_layer(icomplex, :, :) = evolved_ranked_candidate_soc_layer(icomplex, :, :);
                        
                        max_ss_cost_record(icomplex) =  max_ss_cost_record_update(icomplex);
                    end
                end
            end
            
            % STEP 5: shuffle complexes
            candidate_para_value = ranked_candidate_para_value;
            candidate_cost_value = ranked_candidate_cost_value;
            candidate_steady_state = ranked_candidate_steady_state;
            candidate_mod_soc = ranked_candidate_mod_soc;
            candidate_soc_stock = ranked_candidate_soc_stock;
            candidate_soc_layer = ranked_candidate_soc_layer;
            
            candidate_r2 = ranked_candidate_r2;
            
            candidate_index = reshape(1:sample_size, [complex_num, point_num]);
            max_ss_cost = max(candidate_cost_value(candidate_steady_state == 1));
            % record into history
            candidate_para_value_hist(:, :, :, (ishuffle+1)) = candidate_para_value;
            candidate_cost_value_hist(:, :, (ishuffle+1)) = candidate_cost_value;
            candidate_steady_state_hist(:, :, (ishuffle+1)) = candidate_steady_state;
            candidate_r2_hist(:, :, (ishuffle+1)) = candidate_r2;
            candidate_mod_soc_hist(:, :, :, (ishuffle+1)) = candidate_mod_soc;
            
            % STEP 6: convergence test
            r2_converge_test = candidate_r2;
            r2_converge_test(candidate_steady_state == 0) = nan;
            para_value_converge_test = candidate_para_value;
            mod_soc_converge_test = candidate_mod_soc;
            soc_stock_converge_test = candidate_soc_stock;
            soc_layer_converge_test = candidate_soc_layer;
            
            for icomplex = 1:complex_num
                for ipoint = 1:point_num
                    if candidate_steady_state(icomplex, ipoint) == 0
                        para_value_converge_test(icomplex, ipoint, :) = nan;
                        mod_soc_converge_test(icomplex, ipoint, :) = nan;
                        soc_stock_converge_test(icomplex, ipoint, :) = nan;
                        soc_layer_converge_test(icomplex, ipoint, :) = nan;
                    end
                end
            end
            
            % population statistics of parameters
            pop_para_std = nan(npara, 1);
            para_mean = nan(npara, 1); % para_mean = reshape(mean(para_value_converge_test, [1, 2], 'omitnan'), [npara, 1]);
            for ipara = 1:npara
                pop_para_std(ipara) = std(reshape(para_value_converge_test(:, :, ipara), [sample_size, 1]), [], 'omitnan');
                para_mean(ipara) = mean(reshape(para_value_converge_test(:, :, ipara), [sample_size, 1]), 'omitnan');
            end
            
            candidate_para_std_hist(:, (ishuffle+1)) = pop_para_std;
            
            % best optimization
            soc_mod_opt = nan(length(wosis_layer_obs), 1); % reshape(mean(mod_soc_converge_test, [1, 2], 'omitnan'), [length(wosis_layer_obs), 1]);
            soc_std_opt = nan(length(wosis_layer_obs), 1);
            for ilayer = 1:length(wosis_layer_obs)
                soc_mod_opt(ilayer) = mean(reshape(mod_soc_converge_test(:, :, ilayer), [sample_size, 1]), 'omitnan');
                soc_std_opt(ilayer) = std(reshape(mod_soc_converge_test(:, :, ilayer), [sample_size, 1]), [], 'omitnan');
            end
            
            soc_layer_opt = nan(soil_decom_num, 1);
            for ilayer = 1:soil_decom_num
                soc_layer_opt(ilayer) = mean(reshape(soc_layer_converge_test(:, :, ilayer), [sample_size, 1]), 'omitnan');
            end
            
            soc_stock_opt = nan(5, 1);
            for icpool = 1:5
                soc_stock_opt(icpool) = mean(reshape(soc_stock_converge_test(:, :, icpool), [sample_size, 1]), 'omitnan');
            end
            
            r2_opt = mean(reshape(r2_converge_test, [sample_size, 1]), 'omitnan');
            
            if (mean(pop_para_std, 'omitnan') < 0.001 && mean(pop_para_std, 'omitnan') > 0 && length(find(candidate_steady_state == 1))/sample_size >= 0.9) || ...
                    (r2_opt >= 0.9 && mean(pop_para_std, 'omitnan') < 0.005 && mean(pop_para_std, 'omitnan') > 0 && length(find(candidate_steady_state == 1))/sample_size >= 0.9)
                terminate_opt = 1;
            else
                terminate_opt = 0;
            end
            
            % save hist record
            if mod(ishuffle, 5) == 0 || ishuffle == shuffle_num || terminate_opt == 1
                % hist record of para values
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.candidate_para_value = candidate_para_value;']);
                % hist record of cost value
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.candidate_cost_value = candidate_cost_value;']);
                % hist record of steady state value
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.candidate_steady_state = candidate_steady_state;']);
                % hist record of modelled soc
                % eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                %     '.candidate_mod_soc_hist = candidate_mod_soc_hist;']);
                % hist record of r2
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.candidate_r2_hist = candidate_r2_hist;']);
                % hist record of std
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.candidate_para_std_hist = candidate_para_std_hist;']);
                
                % obs
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.wosis_layer_obs = wosis_layer_obs;']);
                % mean value of optimised soc
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.soc_mod_opt = soc_mod_opt;']);
                % mean value of simued soc stock
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.soc_stock_opt = soc_stock_opt;']);
                % mean value of simued soc layer
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.soc_layer_opt = soc_layer_opt;']);
                % std value of optimised soc
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.soc_std_opt = soc_std_opt;']);
                % mean value of optimised para value
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.para_mean = para_mean;']);
                % mean value of r2
                eval([model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ...
                    '.r2_opt = r2_opt;']);
                
                save(['/GFPS8p/cess11/taof/datahub/ensemble/output_data/', model_name, '/',...
                    model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.mat'],...
                    [model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id)]);
            end
            
            if terminate_opt == 1
                disp(['convergence/good results reached at shuffle ', num2str(ishuffle)])
                break
            end
        end
        
    catch
        disp(['error happned in profile ', num2str(iprofile)])
    end
    
    disp([datestr(now), ' Profile ', num2str(iprofile), ' has been finished']);
    
end

disp([datestr(now), ' Program finished']);


