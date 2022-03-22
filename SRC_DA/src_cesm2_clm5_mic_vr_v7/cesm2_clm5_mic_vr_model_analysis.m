clear;
clc;

cesm2_case_name = 'sasu_f05_g16_checked_step4';
start_year = 661;
end_year = 680;

model_name = 'cesm2_clm5_mic_vr_v1';

%%
start_id = 333; %36; 333; 17638; %17427; %17600;
end_id = 333; %36; 333; 17638; %17427; % 17600;
is_resubmit = 0;

disp(['Profiles start from ', num2str(start_id), ' to ', num2str(end_id)]);

% parallel setting
% delete(gcp('nocreate'));
% parpool(workers);

%
warning('off');
format long e;


%% paths
% mac
data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/';
cd(['/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/ENSEMBLE/SRC_DA/src_', model_name, '/']);

% server
% data_path = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/';
% cd(['/GFPS8p/cess11/taof/ensemble/src_da/src_', model_name, '/']);
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
if is_resubmit == 0
    profile_collection = (start_id:end_id);
else
    profile_collection = fun_resubmit_profiles(start_id, end_id, profile_num, model_name);
    disp(['Resubmitting profile number: ', num2str(length(profile_collection))]);
end

profile_range = (1:length(profile_collection));

% parfor iprofile_hat = profile_range
for iprofile_hat = profile_range
    
    iprofile = profile_collection(iprofile_hat);
    
    warning('off')
    
    disp([datestr(now), ' Profile ', num2str(iprofile), ' started']);
    
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
    
    % Parameter names and their initial values in MCMC
    para_name = {'bio', 'cryo', 'q10', 'efolding', 'w_scaling', ...
        'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2', 'tau4s3', 'tau4s4', ...
        'mm_const_assim', 'mm_const_decom', ...
        'fcwdl2', 'fl1s1', 'fl2s1', 'fl3s4', 'mic_cue', 'pcue2death', 'pdeath2soc', ...
        'beta'};
    
    % number of paras
    npara = length(para_name);
    para_min = zeros(npara, 1);
    para_max = ones(npara, 1);
    % set initial values of parameters
    para0 = 0.5*ones(npara, 1);
    % para0 = []';
    % clear warning info
    lastwarn('');
    
    allo_slope = [-0.5, -0.1, -0.05, 0,    0.05, 1,  2];
    tau_enz =    [0.01, 0.05, 0.2,   1,    5,    20, 25];
    mic_cue = 0.1:0.1:0.7;
    
    cue_soc_relation = nan(length(allo_slope), length(tau_enz), length(mic_cue));
    
    for iallo = 1: length(allo_slope)
        for itau = 1 : length(tau_enz)
            for icue = 1 : length(mic_cue)
                para0 = [allo_slope(iallo), tau_enz(itau), mic_cue(icue)];
                % original soc_modelled data
                [~, soc_stock, soc_mod] = ...
                    matrix_fun_model_analysis(para0, nbedrock, sand_vector, npp_mean, ...
                    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
                    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
                
                ss_year_loc = find(isnan(soc_stock(:, 12, 5)) == 0, 1, 'last');
                if isempty(ss_year_loc) == 0
                    cue_soc_relation(iallo, itau, icue) = soc_stock(ss_year_loc, 12, 5);
                end
            end
        end
    end
    save('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/mic_model_cue_soc_relationship/cue_soc_relation.mat', 'cue_soc_relation')
    
%     ax1 = axes('Color', 'none');
%     color_scheme_1 = autumn(6);
%     for iallo = 2:7
%         plot_data = reshape(cue_soc_relation(iallo, 1, :), [length(mic_cue), 1])/1000;
%         plot_data(plot_data < 0) = nan;
%         plot(ax1, mic_cue, plot_data, '*-', 'Color', color_scheme_1(iallo-1, :), 'LineWidth', 3, 'MarkerSize', 15)
%         hold on
%     end
%     xlabel('CUE')
%     ylabel('SOC stock (kgc/m2)')
%     ax1.FontSize = 25;
%     ax1.LineWidth = 2;
%     ax1.YScale = 'log';
% 
%     ax2 = axes('Color', 'none');
%     color_scheme_2 = winter(6);
%     for itau = 1:6
%         plot_data = reshape(cue_soc_relation(6, itau, :), [length(mic_cue), 1])/1000;
%         plot_data(plot_data < 0) = nan;
%         plot(ax2, mic_cue, plot_data, '^--', 'Color', color_scheme_2(itau, :), 'LineWidth', 3, 'MarkerSize', 15)
%         hold on
%     end
%     ax2.FontSize = 25;
%     ax2.LineWidth = 2;
%     ax2.YScale = 'log';
%     % Hide the top axes
%     ax2.Visible = 'off';
%     ax2.XTick = [];
%     ax2.YTick = [];
%     
%     % Link two axes together
%     linkprop([ax1,ax2],{'XLim', 'YLim', 'ZLim', 'CameraUpVector', 'CameraPosition', 'CameraTarget'});
%     
%     % Give each one its colormap
%     colormap(ax1, color_scheme_1)
%     colormap(ax2, color_scheme_2)
%     
%     colorbar(ax1, 'Position',[0.9 0.1 0.02 0.815])
%     colorbar(ax2, 'Position',[0.8 0.1 0.02 0.815])
%     fig = gcf;
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 12 10];
%     print('/Users/phoenix/Desktop/cue_soc_relation.jpeg', '-djpeg', '-r300')
%     close
% 
%     plot_data_1 = reshape(cue_soc_relation(3, 1, :), [length(mic_cue), 1])/1000;
%     plot_data_1(plot_data_1 < 0) = nan;
%     plot_data_2 = reshape(cue_soc_relation(6, 6, :), [length(mic_cue), 1])/1000;
%     plot_data_2(plot_data_2 < 0) = nan;
%     plot_data_3 = reshape(cue_soc_relation(5, 1, :), [length(mic_cue), 1])/1000;
%     plot_data_3(plot_data_3 < 0) = nan;
%     plot(mic_cue, plot_data_1, '*-b', 'LineWidth', 3, 'MarkerSize', 15)
%     hold on
%     plot(mic_cue, plot_data_2, '*-r', 'LineWidth', 3, 'MarkerSize', 15)
%     hold on
%     plot(mic_cue, plot_data_3, '*-g', 'LineWidth', 3, 'MarkerSize', 15)
%     xlabel('CUE')
%     ylabel('SOC stock (kgc/m2)')
%     axis = gca;
%     axis.FontSize = 25;
%     axis.LineWidth = 2;
%     axis.YScale = 'log';
%     fig = gcf;
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 12 10];
%     print('/Users/phoenix/Desktop/cue_soc_relation1.jpeg', '-djpeg', '-r300')
%     close
%     
%   
% 
%     close
%     subplot(2, 3, 1)
%     plot(soc_stock(:, 1, 1) - 0)
%     title('DOC')
%     xlabel('Simu Year')
%     ylabel('Stock (gc/m2)')
%     subplot(2, 3, 2)
%     plot(soc_stock(:, 1, 2) - 0)
%     title('MIC')
%     xlabel('Simu Year')
%     ylabel('Stock (gc/m2)')
%     subplot(2, 3, 3)
%     plot(soc_stock(:, 1, 3) - 0)
%     title('ENZ')
%     xlabel('Simu Year')
%     ylabel('Stock (gc/m2)')
%     subplot(2, 3, 4)
%     plot(soc_stock(:, 1, 4) - 0)
%     title('SOC')
%     xlabel('Simu Year')
%     ylabel('Stock (gc/m2)')
%     subplot(2, 3, 5)
%     plot(soc_stock(:, 1, 5) - 0)
%     title('TOTAL')
%     xlabel('Simu Year')
%     ylabel('Stock (gc/m2)')
end 