
clc;
clear;
model_name = 'cesm2_clm5_mic_vr_v7';
%%
disp([datestr(now), ' Programme started']);

data_dir_output = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/OUTPUT_DATA/';
data_dir_input = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/';

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

%%

sample_profile_num = load([data_dir_input, 'wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_representative_profiles.mat']);
sample_profile_num = reshape(sample_profile_num.sample_profile_id, [5000, 1]);

profile_env_info = load([data_dir_input, '/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat']);
profile_env_info = profile_env_info.EnvInfo;

para_name = {'bio', 'cryo', 'q10', 'efolding', 'w_scaling', ...
    'tau4cwd', 'tau4l1', 'tau4l2', 'tau4s1', 'tau4s2_death', 'tau4s2_enz', 'tau4s3', 'tau4s4', ...
    'mm_const_assim', 'mm_const_decom', ...
    'fcwdl2', 'fl1s1', 'fl2s1', 'fl3s4', ...
    'mic_cue', 'pdeath2soc', ...
    'beta', ...
    'allo_slope_mic'};

col_name_list = {'profile_num', 'cue', 'mm_ratio_assim', 'mm_ratio_decom', 'tau4s2_enz', 'allo_slope_mic', ...
    'soc_0_30cm', 'soc_30_100cm', 'soc_100_200cm', 'mic_0_30cm', 'mic_stock', 'enz_stock', 'doc_stock', 'poc_stock', 'soc_total', 'r2'};


da_summary_mic = nan(length(profile_env_info(:, 1)), length(col_name_list));
da_summary_mic_para = nan(length(profile_env_info(:, 1)), length(para_name));
non_mic_summary = nan(length(profile_env_info(:, 1)), 1);
%%

for iprofile_hat = 1:length(sample_profile_num)
    disp([datestr(now), ' processing profile hat num: ', num2str(iprofile_hat)]);
    
    iprofile =  sample_profile_num(iprofile_hat);
    profile_id = profile_env_info(iprofile, 2);
    file_directory = [data_dir_output, model_name, '/', model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), '.mat'];
    
    if isfile(file_directory) == 1
        da_result = load(file_directory);
        eval(['da_result = da_result.', model_name, '_result_profile_', num2str(iprofile), '_', num2str(profile_id), ';']);
        
        para_mean_std_hist = mean(da_result.candidate_para_std_hist, 1, 'omitnan');
        shuffle_num = length(find(isnan(para_mean_std_hist) == 0)) + 1;
        
        if para_mean_std_hist(shuffle_num) <= 0.005
            
            kp = da_result.para_mean;
            
            da_summary_mic_para(iprofile, 1) = kp(1)*(5*10^(-4) - 3*10^(-5)) + 3*10^(-5);
            da_summary_mic_para(iprofile, 2) = kp(2)*(16*10^(-4) - 3*10^(-5)) + 3*10^(-5);
            da_summary_mic_para(iprofile, 3) = kp(3)*(3 - 1.2) + 1.2;
            da_summary_mic_para(iprofile, 4) = kp(4)*(1 - 0) + 0;
            da_summary_mic_para(iprofile, 5) = kp(5)*(5 - 0) + 0;
            da_summary_mic_para(iprofile, 6) = kp(6)*(6 - 1) + 1;
            da_summary_mic_para(iprofile, 7) = kp(7)*(0.1 - 0) + 0;
            da_summary_mic_para(iprofile, 8) = kp(8)*(0.3 - 0.1) + 0.1;
            da_summary_mic_para(iprofile, 9) = 10^(kp(9)*((0) - (-3)) + (-3));
            da_summary_mic_para(iprofile, 10) = 10^(kp(10)*(0 - (-2)) + (-2));
            da_summary_mic_para(iprofile, 11) = 10^(kp(11)*(1.5 - (0)) + (0)); % kp(11)*(25 - 0.5) + 0.5;
            da_summary_mic_para(iprofile, 12) = kp(12)*(1.5 - 0) + 0;
            da_summary_mic_para(iprofile, 13) = 10^(kp(13)*((-2) - (-7)) + (-7));
            da_summary_mic_para(iprofile, 14) = 10^(kp(14)*(4 - 0) + 0);
            da_summary_mic_para(iprofile, 15) = 10^(kp(15)*(9 - 4) + 4);
            da_summary_mic_para(iprofile, 16) = kp(16)*(1 - 0.5) + 0.5;
            da_summary_mic_para(iprofile, 17) = kp(17)*(0.8 - 0.1) + 0.1;
            da_summary_mic_para(iprofile, 18) = kp(18)*(0.8 - 0.2) + 0.2;
            da_summary_mic_para(iprofile, 19) = kp(19)*(0.8 - 0.2) + 0.2;
            da_summary_mic_para(iprofile, 20) = kp(20)*(0.7 - 0.01) + 0.01;
            da_summary_mic_para(iprofile, 21) = kp(21)*(0.9 - 0.1) + 0.1;
            da_summary_mic_para(iprofile, 22) = kp(22)*(0.9999 - 0.5) + 0.5;
            da_summary_mic_para(iprofile, 23) = kp(23)*(1.5 - (0)) + (0);
            
            
            da_summary_mic(iprofile, strcmp(col_name_list, 'profile_num')) = iprofile;
            
            da_summary_mic(iprofile, strcmp(col_name_list, 'doc_stock')) = da_result.soc_stock_opt(1)/1000; % unit kgc/m2
            da_summary_mic(iprofile, strcmp(col_name_list, 'mic_stock')) = da_result.soc_stock_opt(2)/1000; % unit kgc/m2
            da_summary_mic(iprofile, strcmp(col_name_list, 'enz_stock')) = da_result.soc_stock_opt(3)/1000; % unit kgc/m2
            da_summary_mic(iprofile, strcmp(col_name_list, 'poc_stock')) = da_result.soc_stock_opt(4)/1000; % unit kgc/m2
            da_summary_mic(iprofile, strcmp(col_name_list, 'soc_total')) = da_result.soc_stock_opt(5)/1000; % unit kgc/m2
            
            soc_layer = da_result.soc_layer_opt; % unit gc/m3
            soc_0_30cm = sum(soc_layer(1:4).*dz(1:4)) + soc_layer(5)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4)); % unit gc/m2
            soc_0_100cm = sum(soc_layer(1:8).*dz(1:8)) + soc_layer(9)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8)); % unit gc/m2
            soc_0_200cm = sum(soc_layer(1:11).*dz(1:11)) + soc_layer(12)*dz(12)*(2 - zisoi(11))/(zisoi(12) - zisoi(11)); % unit gc/m2
            
            ratio_30cm = soc_0_30cm/da_result.soc_stock_opt(5);
            mic_30cm = da_result.soc_stock_opt(2)*ratio_30cm; % unit gc/m2
            
            doc_30cm = da_result.soc_stock_opt(1)*ratio_30cm/0.3; % unit gC/m3 for assimillation
            poc_30cm = da_result.soc_stock_opt(4)*ratio_30cm/0.3; % unit gC/m3 for decomposition
                        
            da_summary_mic(iprofile, strcmp(col_name_list, 'soc_0_30cm')) =  soc_0_30cm;
            da_summary_mic(iprofile, strcmp(col_name_list, 'mic_0_30cm')) =  mic_30cm;
            da_summary_mic(iprofile, strcmp(col_name_list, 'soc_30_100cm')) = soc_0_100cm - soc_0_30cm;
            da_summary_mic(iprofile, strcmp(col_name_list, 'soc_100_200cm')) = soc_0_200cm - soc_0_100cm;
            
            da_summary_mic(iprofile, strcmp(col_name_list, 'cue')) = da_result.para_mean(strcmp(para_name, 'mic_cue'))*(0.7 - 0.01) + 0.01;
            da_summary_mic(iprofile, strcmp(col_name_list, 'mm_ratio_assim')) = 10^(da_result.para_mean(strcmp(para_name, 'mm_const_assim'))*(4 - 0) + 0)/doc_30cm;
            da_summary_mic(iprofile, strcmp(col_name_list, 'mm_ratio_decom')) = 10^(da_result.para_mean(strcmp(para_name, 'mm_const_decom'))*(9 - 4) + 4)/poc_30cm;
            da_summary_mic(iprofile, strcmp(col_name_list, 'tau4s2_enz')) = 10^(da_result.para_mean(strcmp(para_name, 'tau4s2_enz'))*(1.5 - (0)) + (0));
            da_summary_mic(iprofile, strcmp(col_name_list, 'allo_slope_mic')) = da_result.para_mean(strcmp(para_name, 'allo_slope_mic'))*(1.5 - (0)) + (0);
            
            
            da_summary_mic(iprofile, strcmp(col_name_list, 'r2')) = da_result.r2_opt;
        end
        
    end
    
end

save([data_dir_output, '/da_summary_', model_name, '/', model_name, '_da_summary_mic_para.mat'], 'da_summary_mic_para');
save([data_dir_output, '/da_summary_', model_name, '/', model_name, '_da_summary_mic.mat'], 'da_summary_mic');
disp([datestr(now), ' Programme finished']);
