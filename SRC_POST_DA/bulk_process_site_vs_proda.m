clear;
clc;

cesm2_case_name = 'sasu_f05_g16_checked_step4';
start_year = 661;
end_year = 680;

model_name = 'cesm2_clm5_cen_vr_v2';

time_domain = 'whole_time'; % 'whole_time', 'before_1985', 'after_1985', 'random_half', 'random_half_2'

for icross_valid = 1:10
    
    exp_name = ['exp_pc_cesm2_23_cross_valid_', num2str(icross_valid)];
    %% paths
    % mac
    data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/';
    cd(['/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/ENSEMBLE/SRC_DA/src_', model_name, '/']);
    
    % server
    % data_path = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/';
    % cd(['/GFPS8p/cess11/taof/ensemble/src_da/src_', model_name, '/']);
    
    %% set vertical soil pools
    month_num = 12;
    soil_cpool_num = 7;
    soil_decom_num = 20;
    
    %% load wosis data
    env_info = load([data_path, 'INPUT_DATA/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat'], 'EnvInfo'); % NPP in this file will be used
    env_info = env_info.EnvInfo;
    % layer_info: "profile_id, date, upper_depth, lower_depth, node_depth, soc_layer_weight, soc_stock, bulk_denstiy, is_pedo"
    wosis_profile_info = ncread([data_path, 'INPUT_DATA/wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'], 'soc_profile_info'); % wosis profile info
    wosis_soc_info = ncread([data_path, 'INPUT_DATA/wosis_2019_snap_shot/soc_data_integrate_wosis_2019_snapshot_hugelius_mishra.nc'], 'data_soc_integrate'); % wosis SOC info
    
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
        load([data_path, 'INPUT_DATA/cesm2_simu/spinup_ss/',...
            cesm2_case_name, '_cesm2_ss_4da_', num2str(start_year), '_', num2str(end_year),...
            '_', var_name_list{ivar}, '.mat'], 'var_record_monthly_mean');
        eval([var_name_list_rename{ivar} ' = var_record_monthly_mean;']);
        
    end
    
    % for ilayer = 1:soil_decom_num
    % 	cesm2_simu_input_vector_litter1(:, :, ilayer, :) = cesm2_simu_input_vector_litter1(:, :, ilayer, :).*dz(ilayer);
    % 	cesm2_simu_input_vector_litter2(:, :, ilayer, :) = cesm2_simu_input_vector_litter2(:, :, ilayer, :).*dz(ilayer);
    % 	cesm2_simu_input_vector_litter3(:, :, ilayer, :) = cesm2_simu_input_vector_litter3(:, :, ilayer, :).*dz(ilayer);
    % 	cesm2_simu_input_vector_cwd(:, :, ilayer, :) = cesm2_simu_input_vector_cwd(:, :, ilayer, :).*dz(ilayer);
    % end
    %
    % cesm2_simu_input_sum_litter1 = reshape(sum(cesm2_simu_input_vector_litter1, 3), [384, 576, 12]);
    % cesm2_simu_input_sum_litter2 = reshape(sum(cesm2_simu_input_vector_litter2, 3), [384, 576, 12]);
    % cesm2_simu_input_sum_litter3 = reshape(sum(cesm2_simu_input_vector_litter3, 3), [384, 576, 12]);
    % cesm2_simu_input_sum_cwd = reshape(sum(cesm2_simu_input_vector_cwd, 3), [384, 576, 12]);
    %
    
    cesm2_simu_input_sum_litter1 = reshape(sum(cesm2_simu_input_vector_litter1, 4), [384, 576, soil_decom_num]);
    cesm2_simu_input_sum_litter2 = reshape(sum(cesm2_simu_input_vector_litter2, 4), [384, 576, soil_decom_num]);
    cesm2_simu_input_sum_litter3 = reshape(sum(cesm2_simu_input_vector_litter3, 4), [384, 576, soil_decom_num]);
    cesm2_simu_input_sum_cwd = reshape(sum(cesm2_simu_input_vector_cwd, 4), [384, 576, soil_decom_num]);
    
    clearvars cesm2_simu_input_vector_litter1 cesm2_simu_input_vector_litter2 cesm2_simu_input_vector_litter3 cesm2_simu_input_vector_cwd
    
    %% load results from neural networking
    nn_predict = csvread([data_path, 'OUTPUT_DATA/neural_networking/nn_para_result_', model_name, '_', time_domain, '_', exp_name, '.csv']);
    nn_site_loc = csvread([data_path, 'OUTPUT_DATA/neural_networking/nn_site_loc_', model_name, '_', time_domain, '_', exp_name, '.csv']);
    
    %% Difference between NN and MCMC results
    load([data_path, 'OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name '_para_mean.mat']);
    load([data_path, 'OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name, '_stat_r2.mat']);
    stat_r2 = max(stat_r2, [], 2);
    
    %% load MLE by site-by-site data assimilation
    para_global_mle = load([data_path, 'OUTPUT_DATA/mcmc_summary_', model_name, '/', model_name , '_mcmc_results.mat']);
    para_global_mle = para_global_mle.mcmc_results.para_mle;
    para_global_mle = real(para_global_mle);
    
    para_global_mle_mean = mean(para_global_mle, 1, 'omitnan')';
    
    %% Individual Simulation
    bulk_process_proda = nan(length(nn_site_loc), 5);
    bulk_process_da = nan(length(nn_site_loc), 5);

    
    for iprofile_hat = 1:length(nn_site_loc)
        
        iprofile = nn_site_loc(iprofile_hat);
        warning('off')
        
        disp([datestr(now), ' Processing profile ', num2str(iprofile_hat)]);
        
        profile_id = wosis_profile_info(iprofile, 1);
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
        
        % input_vector_cwd = reshape(cesm2_simu_input_sum_cwd(lat_loc, lon_loc, :), [1, month_num]);
        % input_vector_litter1 = reshape(cesm2_simu_input_sum_litter1(lat_loc, lon_loc, :), [1, month_num]);
        % input_vector_litter2 = reshape(cesm2_simu_input_sum_litter2(lat_loc, lon_loc, :), [1, month_num]);
        % input_vector_litter3 = reshape(cesm2_simu_input_sum_litter3(lat_loc, lon_loc, :), [1, month_num]);
        
        input_vector_cwd = reshape(cesm2_simu_input_sum_cwd(lat_loc, lon_loc, :), [soil_decom_num, 1]);
        input_vector_litter1 = reshape(cesm2_simu_input_sum_litter1(lat_loc, lon_loc, :), [soil_decom_num, 1]);
        input_vector_litter2 = reshape(cesm2_simu_input_sum_litter2(lat_loc, lon_loc, :), [soil_decom_num, 1]);
        input_vector_litter3 = reshape(cesm2_simu_input_sum_litter3(lat_loc, lon_loc, :), [soil_decom_num, 1]);
        
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
        
        %% data cleansing
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
        
        %% Parameterization
        para_name = {'diffus'; 'cryo';...
            'q10';...
            'efolding';...
            'taucwd'; 'taul1'; 'taul2';...
            'tau4s1'; 'tau4s2'; 'tau4s3';...
            'fl1s1'; 'fl2s1'; 'fl3s2'; 'fs1s2'; 'fs1s3'; 'fs2s1'; 'fs2s3'; 'fs3s1'; 'fcwdl2';...
            'w-scaling'; 'beta'};
        
        para_A_loc = (11:19);
        para_K_loc = (5:10);
        para_Xi_loc = [3, 4, 20];
        para_I_loc = 21;
        para_V_loc = [1, 2];
        
        %% Simulation
        % if using default parameterization
        is_default = 0;
        % parameters values
        para_proda = nn_predict(iprofile_hat, :);
        para_da = para_mean(iprofile, :);
        
        % clear warning info
        lastwarn('');
        [~, ~, ~, ~, ~, ~, ~, ~, ~, bulk_A_proda, bulk_K_proda, bulk_V_proda, bulk_Xi_proda, bulk_I_proda] ...
            = fun_bulk_process(para_proda, is_default, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % clear warning info
        lastwarn('');
        [~, ~, ~, ~, ~, ~, ~, ~, ~, bulk_A_da, bulk_K_da, bulk_V_da, bulk_Xi_da, bulk_I_da] ...
            = fun_bulk_process(para_da, is_default, nbedrock, sand_vector, npp_mean, ...
            input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
            altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            disp(['Assilimation_failed_in_Profile_', num2str(iprofile)]);
            % msgid = [];
            continue
        end
        
        % profiles
        if length(wosis_layer_obs) > 1
            bulk_process_proda(iprofile_hat, :) = [bulk_A_proda, bulk_K_proda, bulk_V_proda, bulk_Xi_proda, bulk_I_proda];
            bulk_process_da(iprofile_hat, :) = [bulk_A_da, bulk_K_da, bulk_V_da, bulk_Xi_da, bulk_I_da];
        end
        

    end
    
    bulk_process_site_da_vs_proda.bulk_process_proda = bulk_process_proda;
    bulk_process_site_da_vs_proda.bulk_process_da = bulk_process_da;
    
    save([data_path, 'OUTPUT_DATA/world_simulation_analyses/bulk_process_site_da_vs_proda_', model_name, '_', time_domain, '_', exp_name, '.mat'], 'bulk_process_site_da_vs_proda');

end

disp('Programme Finished');
