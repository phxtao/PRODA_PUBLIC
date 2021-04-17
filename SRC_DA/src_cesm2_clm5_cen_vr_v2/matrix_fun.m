function [cpools, soc_stock, soc_layer] = matrix_fun(kp, nbedrock, sand_vector, npp_mean, ...
    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin)
% define global parameters
% global diffus adv cryo q10 fq10 maxpsi minpsi efolding ...
%     tau4cwd tau4l1 tau4l2 tau4l3 tau4s1 tau4s2 tau4s3 fl1s1 fl2s1 fl3s2 fs2s1 fs3s1 ...
%     fs1s2 fs1s3 fs2s3 fcwdl2 fcwdl3 ins beta4cwd beta4l1 beta4l2 beta4l3 dmax p4cwd ...
%     p4ml p4cl p4ll
is_default = 0; 

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

xiw = soil_water_profile*w_scaling;
xiw(xiw > 1) = 1;

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

% calculate annual carbon input, sum monthly input (gc/m3/month) as annual input (gc/m3/year)
input_vector_cwd = sum(input_vector_cwd, 2, 'omitnan');
input_vector_litter1 = sum(input_vector_litter1, 2, 'omitnan');
input_vector_litter2 = sum(input_vector_litter2, 2, 'omitnan');
input_vector_litter3 = sum(input_vector_litter3, 2, 'omitnan');

if use_beta == 1
    input_tot_cwd = input_vector_cwd;
    input_tot_litter1 = input_vector_litter1;
    input_tot_litter2 = input_vector_litter2;
    input_tot_litter3 = input_vector_litter3;
    
    matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_year); % litter input gc/m3/day
    matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_year);
    matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_year);
    matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_year);
    matrix_in(81:140,1) = 0;

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
cpools = [cpools(1:20), cpools(21:40), cpools(41:60), cpools(61:80), cpools(81:100), cpools(101:120), cpools(121:140)];
cpools_layer = sum(cpools, 2);
soc_layer = cpools(:, 5:7);
soc_layer = [soc_layer, sum(soc_layer, 2)]; % unit gC/m3

soc_stock = sum(soc_layer(:, 4).*dz(1:n_soil_layer)); % unit gC/m2
end

% ----- CENTURY T response function
function catanf_results = catanf(t1)
catanf_results = 11.75 +(29.7 / pi) * atan( pi * 0.031  * ( t1 - 15.4 ));
end



