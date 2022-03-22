function [total_res_time, total_res_time_baseline, soc_res_time, soc_res_time_baseline, t_scalar, w_scalar] = fun_tau_and_scalar(kp, npp_mean, ins_profile,...
    fraction_cwd_profile, fraction_litter1_profile, fraction_litter2_profile, fraction_litter3_profile,...
    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin)
% define global parameters
% global diffus adv cryo q10 fq10 maxpsi minpsi efolding ...
%     tau4cwd tau4l1 tau4l2 tau4l3 tau4s1 tau4s2 tau4s3 fl1s1 fl2s1 fl3s2 fs2s1 fs3s1 ...
%     fs1s2 fs1s3 fs2s3 fcwdl2 fcwdl3 ins beta4cwd beta4l1 beta4l2 beta4l3 dmax p4cwd ...
%     p4ml p4cl p4ll
global kelvin_to_celsius
kelvin_to_celsius = 273.15;
global use_vertsoilc nspools nspools_vr n_soil_layer days_per_year secspday
use_vertsoilc = 1 ; % whether or not use vertical maxing part
nspools = 7;  % number of pools if no vertical
nspools_vr = 70; % number of pools if vertical
n_soil_layer = 10;  % number of soil layers
days_per_year = 365;
secspday = 24*60*60;
% dt = secspday*30;
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
bio = kp(1)*(5*10^(-4) - 3*10^(-5)) + 3*10^(-5); % exp(kp(1)*(-7 - (-12)) + (-12)); % diffusion (bioturbation?) (m2/yr)
cryo = kp(2)*(16*10^(-4) - 3*10^(-5)) + 3*10^(-5); % exp(kp(2)*(-7 - (-12)) + (-12)); % cryoturbation (m2/yr)
q10 = kp(3)*(3 - 1.2) + 1.2; % Q10 (unitless)
fq10 = q10; % Q10 when forzen (unitless)
efolding = kp(4)*(1 - 0) + 0; % parameters used in vertical discretization of carbon inputs (metre)

tau4cwd = 1/0.3; % kp(5)*(10 - 1) + 1; % turnover time of CWD (yr) 4.1 in technote
tau4l1 = 1/18.5; % kp(6)*(0.2 - 0) + 0; % tau for metabolic litter (yr) 0.066 in technote
tau4l2 = 1/4.9; %kp(7)*(1 - 0.2) + 0.2; % tau for cellulose litter (yr) 0.25 in technote
tau4l3 = tau4l2; % tau for lignin litter (yr)

tau4s1 = kp(5)*(1 - 0) + 0; % tau for fast SOC (yr)
tau4s2 = kp(6)*(50 - 1) + 1; % tau for slow SOC (yr)
tau4s3 = kp(7)*(1000 - 200) + 200; % tau for passive SOC (yr)

fl1s1 = kp(8)*(0.8 - 0.2) + 0.2;
fl2s1 = kp(9)*(0.8 - 0.2) + 0.2;
fl3s2 = kp(10)*(0.8 - 0.2) + 0.2;
fs1s2 = kp(11)*(0.4 - 0) + 0;
fs1s3 = kp(12)*(0.05 - 0) + 0;
fs2s1 = kp(13)*(0.6 - 0.1) + 0.1; 
fs2s3 = kp(14)*(0.1 - 0) + 0;
fs3s1 = kp(15)*(1 - 0) + 0;
fcwdl2 = kp(16)*(1 - 0.5) + 0.5; 

% vegetation root distribution parameter (beta function) (unitless)
beta = kp(17)*(0.9999 - 0.5) + 0.5;
% maximum water potential (MPa)
maxpsi= kp(18)*(0 - (-0.015)) + (-0.015);

% maximum depth of root (the maximum depth that the NPP and be distributed) (metre)
dmax = 3.8;
% fraction of NPP input (unitless) for NPP (gC/m2/yr)
ins = ins_profile;
p4ll = fraction_litter3_profile; % allocation ratio from NPP to lignin litter (unitless)
p4ml = fraction_litter1_profile; % allocation ratio from NPP to metabolic litter (unitless)
p4cl = fraction_litter2_profile; % allocation ratio from NPP to cellulose litter (unitless)
minpsi= -10; % minimum water potential (MPa)
adv = 0; % parameter for advection (m/yr)

p4cwd = fraction_cwd_profile; % allocation ratio from NPP to lignin litter (unitless)
%% Environmental Scalar (Xi)
% temperature related function xit
% calculate rate constant scalar for soil temperature
% assuming that the base rate constants are assigned for non-moisture
% limiting conditions at 25 C.
xit = nan(n_soil_layer, 1);
for ilayer = 1 : n_soil_layer
    if soil_temp_profile(ilayer) >= 0 + kelvin_to_celsius
        xit(ilayer) = q10^((soil_temp_profile(ilayer) - (kelvin_to_celsius + 25))/10);
    else
        xit(ilayer) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer) - (0 + kelvin_to_celsius))/10));
    end
    xit(ilayer) = min(xit(ilayer), 1);
end

% water related function xiw
% calculate the rate constant scalar for soil water content.
% Uses the log relationship with water potential given in
% Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
% a comparison of models. Ecology, 68(5):1190-1200.
% and supported by data in
% Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
% and soil moisture. Soil Biol. Biochem., 15(4):447-453.
xiw = nan(n_soil_layer, 1);
% maxpsi = -(9.8*10^(-5))*exp((1.54 - 0.0095*sand + 0.0063*(100 - sand - clay))*log(10));
for ilayer = 1 : n_soil_layer
    psi = min(maxpsi, soil_water_profile(ilayer));
    if psi > maxpsi
        xiw(ilayer) = 1;
    elseif psi < minpsi
        xiw(ilayer) = 0.05;
    else
        xiw(ilayer)=(log(minpsi/psi))/(log(minpsi/maxpsi));
    end
    % adapted from matrix ORCHIDEE
    xiw(ilayer) = min(max(0.05, xiw(ilayer)), 1);
end

%% Triangle Matrix, A Matrix and K Matrix
% Allocation matrix
a_ma = a_matrix(fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2);
% decomposition matrix
kk_ma = kk_matrix(xit, xiw, xio, xin, efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
% tri matrix
tri_ma = tri_matrix(altmax_current_profile, altmax_lastyear_profile, bio, adv, cryo);

%% B Matrix
% in the original beta model in Jackson et al 1996, the unit for the depth
% of the soil is cm (dmax*100)

vertical_input = nan(n_soil_layer, 1);
if altmax_lastyear_profile > 0
    for j = 1:n_soil_layer
        if zisoi(j) < dmax && (zisoi(j)-dz(j)) < dmax
            vertical_input(j, 1) = (beta^((zisoi(j)-dz(j))*100) - beta^(zisoi(j)*100))/dz(j);
        elseif zisoi(j) > dmax && (zisoi(j)-dz(j)) > dmax
            vertical_input(j, 1) = 0;
        else
            vertical_input(j, 1) = (beta^((zisoi(j)-dz(j))*100) - beta^(dmax*100))/dz(j);
        end
    end    
else
    vertical_input(1) = 1/dz(1);
    vertical_input(2:end) = 0;
end
vertical_input = dz(1:n_soil_layer).*vertical_input/sum(vertical_input.*dz(1:n_soil_layer));


input_matrix = NaN(70,1);

input_matrix(1:10,1) = p4cwd.*vertical_input;
input_matrix(11:20,1) = p4ml.*vertical_input;
input_matrix(21:30,1) = p4cl.*vertical_input;
input_matrix(31:40,1) = p4ll.*vertical_input;
input_matrix(41:70,1) = 0;

t_scalar = xit;
w_scalar = xiw;

diag_scalar = 1./(xiw.*xit);
diag_scalar = [diag_scalar; diag_scalar; diag_scalar; diag_scalar; diag_scalar; diag_scalar; diag_scalar];
diag_scalar = diag(diag_scalar);

kk_ma_without_scalar = kk_ma*diag_scalar;

residence_time = ((a_ma*kk_ma - tri_ma)\(-input_matrix))/days_per_year;
residence_time_baseline = ((a_ma*kk_ma_without_scalar - tri_ma)\(-input_matrix))/days_per_year;

soc_res_time = sum([residence_time(41:50), residence_time(51:60), residence_time(61:70)], 2);
soc_res_time_baseline = sum([residence_time_baseline(41:50), residence_time_baseline(51:60), residence_time_baseline(61:70)], 2);

total_res_time =  sum([residence_time(1:10), residence_time(11:20), residence_time(21:30), residence_time(31:40),...
    residence_time(41:50), residence_time(51:60), residence_time(61:70)], 2);
total_res_time_baseline =  sum([residence_time_baseline(1:10), residence_time_baseline(11:20), residence_time_baseline(21:30), residence_time_baseline(31:40),...
    residence_time_baseline(41:50), residence_time_baseline(51:60), residence_time_baseline(61:70)], 2);

end