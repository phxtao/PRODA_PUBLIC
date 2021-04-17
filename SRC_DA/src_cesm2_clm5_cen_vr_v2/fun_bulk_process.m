function [carbon_input, cpools, cpools_layer, soc_layer, total_res_time, total_res_time_base, res_time_base_pools, t_scalar, w_scalar, bulk_A, bulk_K, bulk_V, bulk_Xi, bulk_I] = ...
    fun_bulk_process(kp, is_default, nbedrock, sand_vector, npp_mean, ...
    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin)
% define global parameters
% global diffus adv cryo q10 fq10 maxpsi minpsi efolding ...
%     tau4cwd tau4l1 tau4l2 tau4l3 tau4s1 tau4s2 tau4s3 fl1s1 fl2s1 fl3s2 fs2s1 fs3s1 ...
%     fs1s2 fs1s3 fs2s3 fcwdl2 fcwdl3 ins beta4cwd beta4l1 beta4l2 beta4l3 dmax p4cwd ...
%     p4ml p4cl p4ll
use_beta = 1;
month_num = 12;
normalize_q10_to_century_tfunc = false;

global kelvin_to_celsius
kelvin_to_celsius = 273.15;

global use_vertsoilc npool npool_vr n_soil_layer days_per_year secspday
use_vertsoilc = 1 ; % whether or not use vertical maxing part
npool = 7;  % number of pools if no vertical
npool_vr = 140; % number of pools if vertical
n_soil_layer = 20;  % number of soil layers
days_per_year = 365;
secspday = 24*60*60;
% dt = secspday*30;

global max_altdepth_cryoturbation max_depth_cryoturb
max_altdepth_cryoturbation = 2;
max_depth_cryoturb = 3;

global dz dz_node zisoi zsoi
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
zisoi_0 = 0;
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

% define parameter names
if is_default == 0
    % diffusion (bioturbation) 10^(-4) (m2/yr)
    bio = kp(1)*(5*10^(-4) - 3*10^(-5)) + 3*10^(-5);
    % cryoturbation 5*10^(-4) (m2/yr)
    cryo = kp(2)*(16*10^(-4) - 3*10^(-5)) + 3*10^(-5);
    % Q10 (unitless) 1.5
    q10 = kp(3)*(3 - 1.2) + 1.2;
    % Q10 when forzen (unitless) 1.5
    fq10 = q10;
    % parameters used in vertical discretization of carbon inputs 10 (metre)
    efolding = kp(4)*(1 - 0) + 0;
    % turnover time of CWD (yr) 3.3333
    tau4cwd = kp(5)*(6 - 1) + 1;
    % tau for metabolic litter (yr) 0.0541
    tau4l1 = kp(6)*(0.11 - 0) + 0;
    % tau for cellulose litter (yr) 0.2041
    tau4l2 = kp(7)*(0.3 - 0.1) + 0.1;
    % tau for lignin litter (yr)
    tau4l3 = tau4l2;
    % tau for fast SOC (yr) 0.1370
    tau4s1 = kp(8)*(1 - 0) + 0;
    % tau for slow SOC (yr) 5
    tau4s2 = kp(9)*(50 - 1) + 1;
    % tau for passive SOC (yr) 222.222
    tau4s3 = kp(10)*(1000 - 200) + 200;
    
    % fraction from l1 to s2, 0.45
    fl1s1 = kp(11)*(0.8 - 0.1) + 0.1;
    % fraction from l2 to s1, 0.5
    fl2s1 = kp(12)*(0.8 - 0.2) + 0.2;
    % fraction from l3 to s2, 0.5
    fl3s2 = kp(13)*(0.8 - 0.2) + 0.2;
    % fraction from s1 to s2, sand dependeted
    fs1s2 = kp(14)*(0.4 - 0) + 0;
    % fraction from s1 to s3, sand dependeted
    fs1s3 = kp(15)*(0.05 - 0) + 0;
    % fraction from s2 to s1, 0.42
    fs2s1 = kp(16)*(0.74 - 0.1) + 0.1;
    % fraction from s2 to s3, 0.03
    fs2s3 = kp(17)*(0.1 - 0) + 0;
    % fraction from s3 to s1, 0.45
    fs3s1 = kp(18)*(0.9 - 0) + 0;
    % fraction from cwd to l2, 0.76
    fcwdl2 = kp(19)*(1 - 0.5) + 0.5;
    
    % water scaling factor
    w_scaling = kp(20)*(5 - 0) + 0;
    
    % beta to describe the shape of vertical profile
    % beta = 0.95;
    beta = kp(21)*(0.9999 - 0.5) + 0.5;
else
    
    % default_para_dir = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/ENSEMBLE/INPUT_DATA/cesm2_parameters/clm5_params.c190829.nc';
    default_para_dir = '/GFPS8p/cess11/taof/datahub/ensemble/input_data/cesm2_parameters/clm5_params.c190829.nc';
	% define parameter names
    bio = ncread(default_para_dir, 'som_diffus')*3600*24*365; % diffusion (bioturbation?) (m2/yr)
    cryo = ncread(default_para_dir, 'cryoturb_diffusion_k')*3600*24*365; % exp(kp(2)*(-7 - (-12)) + (-12)); % cryoturbation (m2/yr)
    q10 = ncread(default_para_dir, 'q10_hr'); % Q10 (unitless)
    fq10 = ncread(default_para_dir, 'froz_q10'); % Q10 when forzen (unitless)
    efolding = ncread(default_para_dir, 'decomp_depth_efolding'); % parameters used in vertical discretization of carbon inputs (metre)
    
    tau4cwd = ncread(default_para_dir, 'tau_cwd'); % kp(5)*(10 - 1) + 1; % turnover time of CWD (yr) 4.1 in technote
    tau4l1 = ncread(default_para_dir, 'tau_l1'); % kp(6)*(0.2 - 0) + 0; % tau for metabolic litter (yr) 0.066 in technote
    tau4l2 = ncread(default_para_dir, 'tau_l2_l3'); %kp(7)*(1 - 0.2) + 0.2; % tau for cellulose litter (yr) 0.25 in technote
    tau4l3 = tau4l2; % tau for lignin litter (yr)
    
    tau4s1 = ncread(default_para_dir, 'tau_s1'); % tau for fast SOC (yr)
    tau4s2 = ncread(default_para_dir, 'tau_s2'); % tau for slow SOC (yr)
    tau4s3 = ncread(default_para_dir, 'tau_s3'); % tau for passive SOC (yr)
    
    fl1s1 = 1 - ncread(default_para_dir, 'rf_l1s1_bgc'); % 1 - 0.55;
    fl2s1 = 1 - ncread(default_para_dir, 'rf_l2s1_bgc'); % 1 - 0.5;
    fl3s2 = 1 - ncread(default_para_dir, 'rf_l3s2_bgc'); % 1 - 0.5;
    fs1s2 = nan;
    fs1s3 = nan;
    fs2s1 = (0.42/0.45)*(1 - ncread(default_para_dir, 'rf_s2s1_bgc')); % (0.42/0.45)*(1 - 0.55);
    fs2s3 = (0.03/0.45)*(1 - ncread(default_para_dir, 'rf_s2s3_bgc')); % (0.03/0.45)*(1 - 0.55);
    fs3s1 = 1 - ncread(default_para_dir, 'rf_s3s1_bgc'); % 1 - 0.55;
    fcwdl2 = 1 - ncread(default_para_dir, 'cwd_flig'); % 0.76;
    
    beta = 0.95;
end
% maximum and minimum water potential (MPa)
maxpsi= -0.0020;

minpsi= -2; % minimum water potential (MPa)

adv = 0; % parameter for advection (m/yr)

%% Environmental Scalar (Xi)
xit = nan(n_soil_layer, month_num);
xiw = nan(n_soil_layer, month_num);

for imonth = 1:month_num
    % temperature related function xit
    % calculate rate constant scalar for soil temperature
    % assuming that the base rate constants are assigned for non-moisture
    % limiting conditions at 25 C.
    for ilayer = 1 : n_soil_layer
        if soil_temp_profile(ilayer, imonth) >= 0 + kelvin_to_celsius
            xit(ilayer, imonth) = q10^((soil_temp_profile(ilayer, imonth) - (kelvin_to_celsius + 25))/10);
        else
            xit(ilayer, imonth) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer, imonth) - (0 + kelvin_to_celsius))/10));
        end
    end
    
    catanf_30 = catanf(30);
    normalization_tref = 15;
    if normalize_q10_to_century_tfunc == true
        % scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
        normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10^((normalization_tref-25)/10));
        xit = xit*normalization_factor;
    end
    
    % water related function xiw
    % calculate the rate constant scalar for soil water content.
    % Uses the log relationship with water potential given in
    % Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
    % a comparison of models. Ecology, 68(5):1190-1200.
    % and supported by data in
    % Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
    % and soil moisture. Soil Biol. Biochem., 15(4):447-453.
    
    %     for ilayer = 1 : n_soil_layer
    %         psi = min(maxpsi, soil_water_profile(ilayer, imonth));
    %
    %         if psi > minpsi
    %             xiw(ilayer, imonth) = (log(minpsi/psi))/(log(minpsi/maxpsi));
    %         else
    %             xiw(ilayer, imonth) = 0;
    %         end
    %         xiw(ilayer, imonth) = max(xiw(ilayer, imonth), 0.05);
    %     end
    
end

% xit = soil_temp_profile;

if is_default == 0
    xiw = soil_water_profile*w_scaling;
    xiw(xiw > 1) = 1;
else
    xiw = soil_water_profile;
end

%% Triangle Matrix, A Matrix and K Matrix
sand_vector = mean(sand_vector, 2, 'omitnan');
% Allocation matrix
a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector);

kk_ma_middle = nan(npool_vr, npool_vr, month_num);
tri_ma_middle = nan(npool_vr, npool_vr, month_num);
for imonth = 1:month_num
    % decomposition matrix
    monthly_xit = xit(:, imonth);
    monthly_xiw = xiw(:, imonth);
    monthly_xio = xio(:, imonth);
    monthly_xin = xin(:, imonth);
    kk_ma_middle(:, :, imonth) = kk_matrix(monthly_xit, monthly_xiw, monthly_xio, monthly_xin, efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
    % tri matrix
    monthly_nbedrock = nbedrock(imonth);
    monthly_altmax_current_profile = altmax_current_profile(imonth);
    monthly_altmax_lastyear_profile = altmax_lastyear_profile(imonth);
    tri_ma_middle(:, :, imonth) = tri_matrix(monthly_nbedrock, monthly_altmax_current_profile, monthly_altmax_lastyear_profile, bio, adv, cryo);
end
tri_ma = mean(tri_ma_middle, 3, 'omitnan');
kk_ma = mean(kk_ma_middle, 3, 'omitnan');


%% Vertical Profile
% in the original beta model in Jackson et al 1996, the unit for the depth
% of the soil is cm (dmax*100)

m_to_cm = 100;
vertical_prof = nan(n_soil_layer, 1);

if altmax_lastyear_profile > 0
    for j = 1:n_soil_layer
        if j == 1
            vertical_prof(j) = (beta^((zisoi_0)*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
        else
            vertical_prof(j) = (beta^((zisoi(j - 1))*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
        end
    end
else
    vertical_prof(1) = 1/dz(1);
    vertical_prof(2:end) = 0;
end

vertical_input = dz(1:n_soil_layer).*vertical_prof/sum(vertical_prof.*dz(1:n_soil_layer));

%% Analytical Solution
matrix_in=nan(140,1);

% % calculate annual carbon input, sum monthly input (gc/m3/month) as annual input (gc/m3/year)
% input_vector_cwd = sum(input_vector_cwd, 2, 'omitnan');
% input_vector_litter1 = sum(input_vector_litter1, 2, 'omitnan');
% input_vector_litter2 = sum(input_vector_litter2, 2, 'omitnan');
% input_vector_litter3 = sum(input_vector_litter3, 2, 'omitnan');

% calculate annual carbon input, (gc/m3/year)

if use_beta == 1
    input_tot_cwd = sum(input_vector_cwd.*dz(1:n_soil_layer)); % (gc/m2/year)
    input_tot_litter1 = sum(input_vector_litter1.*dz(1:n_soil_layer));
    input_tot_litter2 = sum(input_vector_litter2.*dz(1:n_soil_layer));
    input_tot_litter3 = sum(input_vector_litter3.*dz(1:n_soil_layer));
    
    if is_default == 0
        matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_year); % litter input gc/m3/day
        matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_year);
        matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_year);
        matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_year);
        matrix_in(81:140,1) = 0;
    else
        matrix_in(1:20,1) = input_vector_cwd./(days_per_year); % litter input gc/m3/day
        matrix_in(21:40,1) = input_vector_litter1./(days_per_year);
        matrix_in(41:60,1) = input_vector_litter2./(days_per_year);
        matrix_in(61:80,1) = input_vector_litter3./(days_per_year);
        matrix_in(81:140,1) = 0;
    end
    carbon_input = input_tot_cwd + input_tot_litter1 + input_tot_litter2 + input_tot_litter3;
else
    matrix_in(1:20,1) = input_vector_cwd./days_per_year; % litter input gc/m3/day
    matrix_in(21:40,1) = input_vector_litter1./days_per_year;
    matrix_in(41:60,1) = input_vector_litter2./days_per_year;
    matrix_in(61:80,1) = input_vector_litter3./days_per_year;
    matrix_in(81:140,1) = 0;
    
    carbon_input = sum((input_vector_cwd + input_vector_litter1 + input_vector_litter2 + input_vector_litter3).*dz(1:n_soil_layer));
end

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]';
site_input_ratio = carbon_input/(sum(npp_mean.*days_in_month*24*3600));

cpools = (a_ma*kk_ma-tri_ma)\(-matrix_in);
cpools(cpools < 0) = 0;
cpools_layer = sum([cpools(1:20), cpools(21:40), cpools(41:60), cpools(61:80), cpools(81:100), cpools(101:120), cpools(121:140)], 2);
soc_layer = sum([cpools(81:100), cpools(101:120), cpools(121:140)], 2);


input_matrix = nan(npool_vr,1);

if is_default == 0
    input_matrix(1:20,1) = input_tot_cwd/carbon_input.*vertical_input./dz(1:n_soil_layer);
    input_matrix(21:40,1) = input_tot_litter1/carbon_input.*vertical_input./dz(1:n_soil_layer);
    input_matrix(41:60,1) = input_tot_litter2/carbon_input.*vertical_input./dz(1:n_soil_layer);
    input_matrix(61:80,1) = input_tot_litter3/carbon_input.*vertical_input./dz(1:n_soil_layer);
    input_matrix(81:140,1) = 0;
else
	vertical_input = (input_vector_cwd + input_vector_litter1 + input_vector_litter2 + input_vector_litter3).*dz(1:n_soil_layer)./(input_tot_cwd + input_tot_litter1 + input_tot_litter2 + input_tot_litter3);
    input_matrix(1:20,1) = input_tot_cwd/carbon_input.*input_vector_cwd.*dz(1:n_soil_layer)./input_tot_cwd./dz(1:n_soil_layer);
    input_matrix(21:40,1) = input_tot_litter1/carbon_input.*input_vector_litter1.*dz(1:n_soil_layer)./input_tot_litter1./dz(1:n_soil_layer);
    input_matrix(41:60,1) = input_tot_litter2/carbon_input.*input_vector_litter2.*dz(1:n_soil_layer)./input_tot_litter2./dz(1:n_soil_layer);
    input_matrix(61:80,1) = input_tot_litter3/carbon_input.*input_vector_litter3.*dz(1:n_soil_layer)./input_tot_litter3./dz(1:n_soil_layer);
    input_matrix(81:140,1) = 0;
end

diag_scalar_middle = nan(n_soil_layer, month_num);
for imonth = 1:month_num
    diag_scalar_middle(:, imonth) = 1./(xiw(:, imonth).*xit(:, imonth).*xio(:, imonth).*xin(:, imonth));
end
t_scalar = mean(xit, 2, 'omitnan');
w_scalar = mean(xiw, 2, 'omitnan');

diag_scalar = mean(diag_scalar_middle, 2, 'omitnan');
diag_scalar = [diag_scalar; diag_scalar; diag_scalar; diag_scalar; diag_scalar; diag_scalar; diag_scalar];
diag_scalar = diag(diag_scalar);


residence_time = ((a_ma*kk_ma - tri_ma)\(-input_matrix))/days_per_year;
residence_time_baseline = (((a_ma*kk_ma - tri_ma)*diag_scalar)\(-input_matrix))/days_per_year;
% residence_time_baseline = ((a_ma*kk_ma - tri_ma))\(-input_matrix)./diag_scalar/days_per_year;

total_res_time =  sum([residence_time(81:100), residence_time(101:120), residence_time(121:140)], 2);


% sum(total_res_time.*dz(1:20)) - sum(soc_layer(1:20).*dz(1:20))/carbon_input


total_res_time_base = sum([residence_time_baseline(81:100), residence_time_baseline(101:120), residence_time_baseline(121:140)], 2);

res_time_base_pools = residence_time_baseline;


%% bulk process
% I
cum_fraction_input = cumsum(vertical_input);
cum_fraction_input(cum_fraction_input >= 1) = nan;
bulk_I = mean(exp(log(1 - cum_fraction_input)./(zsoi(1:n_soil_layer)*100)), 'omitnan');

% K
% cpools_total = ... % soc stock unit gc/m2
%     sum([cpools(1:20).*dz(1:n_soil_layer), ...
%     cpools(21:40).*dz(1:n_soil_layer), ...
%     cpools(41:60).*dz(1:n_soil_layer), ...
%     cpools(61:80).*dz(1:n_soil_layer), ...
%     cpools(81:100).*dz(1:n_soil_layer), ...
%     cpools(101:120).*dz(1:n_soil_layer), ...
%     cpools(121:140).*dz(1:n_soil_layer)], 1);
% decom_cpool = 1./[tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3]; % unit yr-1
% bulk_K = sum(decom_cpool.*(cpools_total/(sum(cpools_total)))); % unit yr-1

cpools_total = ... % soc stock unit gc/m2
    sum([cpools(81:100).*dz(1:n_soil_layer), ...
    cpools(101:120).*dz(1:n_soil_layer), ...
    cpools(121:140).*dz(1:n_soil_layer)], 1);
decom_cpool = 1./[tau4s1, tau4s2, tau4s3]; % unit yr-1
bulk_K = sum(decom_cpool.*(cpools_total/(sum(cpools_total)))); % unit yr-1

% xi
depth_scalar = exp(-zsoi(1:n_soil_layer)/efolding);
bulk_Xi_layer = t_scalar.*w_scalar.*depth_scalar;

bulk_Xi = sum(bulk_Xi_layer.*(soc_layer.*dz(1:n_soil_layer))/sum((soc_layer.*dz(1:n_soil_layer))));
% bulk_Xi = mean(bulk_Xi_layer.^(1/3), 'omitnan');

if isreal(bulk_Xi) == 0
    bulk_Xi = nan;
end

% A
if is_default == 1
    sand_vector(sand_vector < 0) = 0;
    sand_vector(sand_vector > 100) = 100;
    
    t = 0.85 - 0.68*0.01*(100 - sand_vector);
    f_s1s2 = 1 - 0.004./(1 - t);
    f_s1s3 = 0.004./(1 - t);
    
    rf_s1s2 = t;
    rf_s1s3 = t;
    
    t = f_s1s2;
    
    fs1s2_vector = (1 - rf_s1s2).*t; % fs1s2
    fs1s3_vector = (1 - rf_s1s3).*(1 - t); % fs1s3
    
    fs1s2 = mean(fs1s2_vector);
    fs1s3 = mean(fs1s3_vector);
end

% donor_pool_layer = [cpools(1:20), cpools(21:40), cpools(41:60), cpools(61:80), cpools(81:100), cpools(101:120), cpools(121:140)];
%
% cue_cpool = [fl1s1, fl2s1, fl3s2, ...
%     fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, ...
%     fcwdl2, 1 - fcwdl2];
% donor_pool_size = [donor_pool_layer(:, 2), donor_pool_layer(:, 3), donor_pool_layer(:, 4), ...
%     donor_pool_layer(:, 5), donor_pool_layer(:, 5), donor_pool_layer(:, 6), donor_pool_layer(:, 6), donor_pool_layer(:, 7), ...
%     donor_pool_layer(:, 1), donor_pool_layer(:, 1)];
% donor_decomp = [decom_cpool(2), decom_cpool(3), decom_cpool(4), ...
%     decom_cpool(5), decom_cpool(5), decom_cpool(6), decom_cpool(6), decom_cpool(7), ...
%     decom_cpool(1), decom_cpool(1)];


donor_pool_layer = [cpools(1:20), cpools(21:40), cpools(41:60), cpools(61:80), cpools(81:100), cpools(101:120), cpools(121:140)];

cue_cpool = [fl1s1, fl2s1, fl3s2, ...
    fs1s2, fs1s3, fs2s1, fs2s3, fs3s1];
donor_pool_size = [donor_pool_layer(:, 2), donor_pool_layer(:, 3), donor_pool_layer(:, 4), ...
    donor_pool_layer(:, 5), donor_pool_layer(:, 5), donor_pool_layer(:, 6), donor_pool_layer(:, 6), donor_pool_layer(:, 7)];
decom_cpool = 1./[tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3]; % unit yr-1
donor_decomp = [decom_cpool(2), decom_cpool(3), decom_cpool(4), ...
    decom_cpool(5), decom_cpool(5), decom_cpool(6), decom_cpool(6), decom_cpool(7)];

total_doner_flow = nan(1, length(cue_cpool));
for idoner = 1:length(cue_cpool)
    total_doner_flow(idoner) = sum(donor_pool_size(:, idoner)*donor_decomp(idoner).*bulk_Xi_layer.*dz(1:n_soil_layer));
end

bulk_A = sum(cue_cpool.*total_doner_flow)/sum(total_doner_flow([1, 2, 3, 4, 6, 8]));


% V
bulk_V_monthly = nan(month_num, 1);
for imonth = 1:month_num
    bulk_V_middle = diag(tri_ma_middle(21:40, 21:40, imonth))*365; % unit from day-1 to yr-1
    bulk_V_monthly(imonth) = sum(bulk_V_middle.*(soc_layer.*dz(1:n_soil_layer))/sum((soc_layer.*dz(1:n_soil_layer))));
end

% for imonth = 1:month_num
%     % tri matrix
%     monthly_nbedrock = nbedrock(imonth);
%     altmax = altmax_current_profile(imonth);
%     altmax_lastyear = altmax_lastyear_profile(imonth);
%     som_diffus = bio/days_per_year;
%     som_adv_flux = adv/days_per_year;
%     cryoturb_diffusion_k = cryo/days_per_year;
%
%     %------ first get diffusivity / advection terms -------%
%     som_adv_coef = zeros(n_soil_layer+1, 1);             	% SOM advective flux (m/day)
%     som_diffus_coef = zeros(n_soil_layer+1, 1);               % SOM diffusivity due to bio/cryo-turbation (m2/day)
%
%     % use different mixing rates for bioturbation and cryoturbation, with fixed bioturbation and cryoturbation set to a maximum depth
%     if  (( max(altmax, altmax_lastyear) <= max_altdepth_cryoturbation ) && ...
%             ( max(altmax, altmax_lastyear) > 0.) )
%         % use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
%         for j = 1 : n_soil_layer+1
%             if ( j <= monthly_nbedrock+1 )
%                 if ( zisoi(j) < max(altmax, altmax_lastyear) )
%                     som_diffus_coef(j) = cryoturb_diffusion_k;
%                     som_adv_coef(j) = 0.;
%                 else
%                     som_diffus_coef(j) = max(cryoturb_diffusion_k * ...
%                         ( 1. - ( zisoi(j) - max(altmax, altmax_lastyear) ) / ...
%                         ( min(max_depth_cryoturb, zisoi(monthly_nbedrock+1)) - max(altmax, altmax_lastyear) ) ), 0.);  % go linearly to zero between ALT and max_depth_cryoturb
%                     som_adv_coef(j) = 0.;
%                 end
%             else
%                 som_adv_coef(j) = 0.;
%                 som_diffus_coef(j) = 0.;
%             end
%         end
%     elseif (  max(altmax, altmax_lastyear) > 0. )
%         % constant advection, constant diffusion
%         for j = 1 : n_soil_layer+1
%             if ( j <= monthly_nbedrock+1 )
%                 som_adv_coef(j) = som_adv_flux;
%                 som_diffus_coef(j) = som_diffus;
%             else
%                 som_adv_coef(j) = 0.;
%                 som_diffus_coef(j) = 0.;
%             end
%         end
%     else
%         % completely frozen soils--no mixing
%         for j = 1 : n_soil_layer+1
%             som_adv_coef(j) = 0.;
%             som_diffus_coef(j) = 0.;
%         end
%     end
%
%     bulk_V_monthly(imonth) = sum(som_diffus_coef(1:n_soil_layer).*(dz(1:n_soil_layer)./(sum(dz(1:n_soil_layer)))));
%
% end
%
bulk_V = mean(bulk_V_monthly, 'omitnan');

end

% ----- CENTURY T response function
function catanf_results = catanf(t1)
catanf_results = 11.75 +(29.7 / pi) * atan( pi * 0.031  * ( t1 - 15.4 ));
end
