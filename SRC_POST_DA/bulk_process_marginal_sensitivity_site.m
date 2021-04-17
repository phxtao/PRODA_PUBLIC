clear;
clc;
process_id = 3;

process_manage_list = {'A', 'K', 'Xi', 'V', 'I', 'NPP'};

%%

cesm2_case_name = 'sasu_f05_g16_checked_step4';
start_year = 661;
end_year = 680;

model_name = 'cesm2_clm5_cen_vr_v2';

time_domain = 'whole_time'; % 'whole_time', 'before_1985', 'after_1985', 'random_half', 'random_half_2'

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
manage_proposal = -0.5:0.1:0.5;

bulk_process_A = nan(length(para_mean), length(manage_proposal));
bulk_process_I = nan(length(para_mean), length(manage_proposal));
bulk_process_K = nan(length(para_mean), length(manage_proposal));
bulk_process_V = nan(length(para_mean), length(manage_proposal));
bulk_process_Xi = nan(length(para_mean), length(manage_proposal));
bulk_process_NPP = nan(length(para_mean), length(manage_proposal));

Var_Decom_Mean_30cm = nan(length(para_mean), length(manage_proposal));
Var_Decom_Mean_100cm = nan(length(para_mean), length(manage_proposal));
Var_Decom_Mean_200cm = nan(length(para_mean), length(manage_proposal));
Var_Decom_Mean_Total = nan(length(para_mean), length(manage_proposal));

Var_Decom_Tau_30cm = nan(length(para_mean), length(manage_proposal));
Var_Decom_Tau_100cm = nan(length(para_mean), length(manage_proposal));
Var_Decom_Tau_200cm = nan(length(para_mean), length(manage_proposal));
Var_Decom_Tau_Total = nan(length(para_mean), length(manage_proposal));

Var_Decom_Tau_Total_Base = nan(length(para_mean), length(manage_proposal));

for iprofile_hat = 1:length(para_mean)
    
    iprofile = iprofile_hat;
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
    
    if process_id == 1
        para_manage_loc = para_A_loc;
    elseif process_id == 2
        para_manage_loc = para_K_loc;
    elseif process_id == 3
        para_manage_loc = para_Xi_loc;
    elseif process_id == 4
        para_manage_loc = para_V_loc;
    elseif process_id == 5
        para_manage_loc = para_I_loc;
    end
    
    %%
    for imanage = 1:length(manage_proposal)
        is_default = 0;
        
        Para = para_mean(iprofile_hat, :);
        
        if process_id ~= 6
            
            Para(para_manage_loc) = Para(para_manage_loc) + manage_proposal(imanage)*Para(para_manage_loc);
            
            if Para(21) > 0.98
                Para(21) = 0.98;
            end
            
            Para(Para > 1) = 1;
            Para(Para <= 0) = 10^(-4);
            
            [carbon_input, ~, ~, soc_layer, total_res_time, total_res_time_base, ~, ~, ~, bulk_A, bulk_K, bulk_V, bulk_Xi, bulk_I] ...
                = fun_bulk_process(Para, is_default, nbedrock, sand_vector, npp_mean, ...
                input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
                altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        else
            input_vector_litter1_manage = input_vector_litter1 + manage_proposal(imanage)*input_vector_litter1;
            input_vector_litter2_manage = input_vector_litter2 + manage_proposal(imanage)*input_vector_litter2;
            input_vector_litter3_manage = input_vector_litter3 + manage_proposal(imanage)*input_vector_litter3;
            input_vector_cwd_manage = input_vector_cwd + manage_proposal(imanage)*input_vector_cwd;
            npp_mean_manage = npp_mean + manage_proposal(imanage)*npp_mean;
            
            if Para(21) > 0.98
                Para(21) = 0.98;
            end
            
            Para(Para > 1) = 1;
            Para(Para <= 0) = 10^(-4);
            
            [carbon_input, ~, ~, soc_layer, total_res_time, total_res_time_base, ~, ~, ~, bulk_A, bulk_K, bulk_V, bulk_Xi, bulk_I] ...
                = fun_bulk_process(Para, is_default, nbedrock, sand_vector, npp_mean_manage, ...
                input_vector_cwd_manage, input_vector_litter1_manage, input_vector_litter2_manage, input_vector_litter3_manage, ...
                altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
        end
        
        var_loc = imanage;
        
        Var_Decom_Mean_30cm(iprofile_hat, var_loc) = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Mean_100cm(iprofile_hat, var_loc) = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Mean_200cm(iprofile_hat, var_loc) = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Mean_Total(iprofile_hat, var_loc) = sum(soc_layer(1:20).*dz(1:20));
        
        Var_Decom_Tau_30cm(iprofile_hat, var_loc) = sum(total_res_time(1:4).*dz(1:4)) + total_res_time(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
        Var_Decom_Tau_100cm(iprofile_hat, var_loc) = sum(total_res_time(1:8).*dz(1:8)) + total_res_time(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
        Var_Decom_Tau_200cm(iprofile_hat, var_loc) = sum(total_res_time(1:11).*dz(1:11)) + total_res_time(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11));
        Var_Decom_Tau_Total(iprofile_hat,  var_loc) = sum(total_res_time(1:20).*dz(1:20));
        
        Var_Decom_Tau_Total_Base(iprofile_hat,  var_loc) = sum(total_res_time_base.*dz(1:20));
        
        bulk_process_A(iprofile_hat, var_loc) = bulk_A;
        bulk_process_I(iprofile_hat, var_loc) = bulk_I;
        bulk_process_K(iprofile_hat, var_loc) = bulk_K;
        bulk_process_V(iprofile_hat, var_loc) = bulk_V;
        bulk_process_Xi(iprofile_hat, var_loc) = bulk_Xi;
        
        bulk_process_NPP(iprofile_hat, var_loc) = carbon_input;
    end
end

Var_Decom_Grid.GridInfo = GridInfo;

Var_Decom_Grid.Var_Decom_Mean_30cm = Var_Decom_Mean_30cm;
Var_Decom_Grid.Var_Decom_Mean_100cm = Var_Decom_Mean_100cm;
Var_Decom_Grid.Var_Decom_Mean_200cm = Var_Decom_Mean_200cm;
Var_Decom_Grid.Var_Decom_Mean_Total = Var_Decom_Mean_Total;

Var_Decom_Grid.Var_Decom_Tau_30cm = Var_Decom_Tau_30cm;
Var_Decom_Grid.Var_Decom_Tau_100cm = Var_Decom_Tau_100cm;
Var_Decom_Grid.Var_Decom_Tau_200cm = Var_Decom_Tau_200cm;
Var_Decom_Grid.Var_Decom_Tau_Total = Var_Decom_Tau_Total;
Var_Decom_Grid.Var_Decom_Tau_Total_Base = Var_Decom_Tau_Total_Base;

Var_Decom_Grid.bulk_process_A = bulk_process_A;
Var_Decom_Grid.bulk_process_K = bulk_process_K;
Var_Decom_Grid.bulk_process_Xi = bulk_process_Xi;
Var_Decom_Grid.bulk_process_V = bulk_process_V;
Var_Decom_Grid.bulk_process_I = bulk_process_I;
Var_Decom_Grid.bulk_process_NPP = bulk_process_NPP;

save([data_path, 'output_data/world_simulation_analyses/marginal_sensitivity_site_', process_manage_list{process_id}, '_', model_name, '_', nn_exp_name, '.mat'], 'Var_Decom_Grid');

disp('Programme Finished');
