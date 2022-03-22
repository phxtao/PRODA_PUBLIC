function [cpool_current, soc_stock, soc_layer, ss_index] = matrix_fun(kp, nbedrock, sand_vector, npp_mean, ...
    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin)
warning off
% define global parameters
% global diffus adv cryo q10 fq10 maxpsi minpsi efolding ...
%     tau4cwd tau4l1 tau4l2 tau4l3 tau4s1 tau4s2 tau4s3 fl1s1 fl2s1 fl3s2 fs2s1 fs3s1 ...
%     fs1s2 fs1s3 fs2s3 fcwdl2 fcwdl3 ins beta4cwd beta4l1 beta4l2 beta4l3 dmax p4cwd ...
%     p4ml p4cl p4ll

use_beta = 1;
month_num = 12;
normalize_q10_to_century_tfunc = false;
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]';

global kelvin_to_celsius
kelvin_to_celsius = 273.15;

global use_vertsoilc npool npool_vr n_soil_layer days_per_year secspday
use_vertsoilc = 1 ; % whether or not use vertical maxing part
npool = 8;  % number of pools if no vertical
n_soil_layer = 20;  % number of soil layers
npool_vr = n_soil_layer*npool; % number of pools if vertical
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

%-------------------------------------------
% define parameters
%-------------------------------------------
% ------------vertical transport
% define parameter names
% diffusion (bioturbation) 10^(-4) (m2/yr)
bio = kp(1)*(5*10^(-4) - 3*10^(-5)) + 3*10^(-5);
% cryoturbation 5*10^(-4) (m2/yr)
cryo = kp(2)*(16*10^(-4) - 3*10^(-5)) + 3*10^(-5);
% ------------env. scalar
% Q10 (unitless) 1.5
q10 = kp(3)*(3 - 1.2) + 1.2;
% Q10 when forzen (unitless) 1.5
fq10 = q10;
% parameters used in vertical discretization of carbon inputs 10 (metre) 0.5
efolding = kp(4)*(1 - 0) + 0;
% water scaling factor 1
w_scaling = kp(5)*(5 - 0) + 0;
% ------------max decomposition
% turnover time of CWD (yr) 3.3333
tau4cwd = kp(6)*(6 - 1) + 1;
% tau for metabolic litter (yr) 0.0541
tau4l1 = kp(7)*(0.1 - 0) + 0;
% tau for cellulose litter (yr) 0.2041
tau4l2 = kp(8)*(0.3 - 0.1) + 0.1;
% tau for lignin litter (yr)
tau4l3 = tau4l2;
% tau for DOC (yr) 1.14*10^(-12) ~ exp(-27)
tau4s1 = exp(-(kp(9)*(50 - 1) + 1));
% tau for MIC (yr) 0.57 for death
tau4s2_death = kp(10)*(1 - 0.1) + 0.1;
% tau for MIC (yr) 22 for enz production
tau4s2_enz = kp(11)*(30 - 0) + 0;
% tau for ENZ (yr) 0.11
tau4s3 = kp(12)*(0.5 - 0) + 0;
% tau for SOC (yr)  exp(-6.5) for allo slope > 0, exp(-22) for allo slope < 0
tau4s4 = exp(-(kp(13)*(10 - (-6)) + (-6)));
% ------------michaelis-menten
% (michaelis-menten) concentration DOC (yr) for half max assimlation reaction from DOC to MIC, unit: gc/m3, 1*10^(2) ~ exp(4.6)
mm_const_assim = exp(kp(14)*(10 - 0) + 0);
% (michaelis-menten) concentration DOC (yr) for half max decomposition reaction from SOC to MIC, unit: gc.m3, 5*10^(5) ~ exp(13)
mm_const_decom = exp(kp(15)*(20 - 10) + 10);
% ------------transfer fraction intra litter
% fraction from cwd to l2 0.75
fcwdl2 = kp(16)*(1 - 0.5) + 0.5;
% ------------transfer fraction litter -> soil
% fraction from l1 to s2 0.45
fl1s1 = kp(17)*(0.8 - 0.1) + 0.1;
% fraction from l2 to s1 0.5
fl2s1 = kp(18)*(0.8 - 0.2) + 0.2;
% fraction from l3 to s2 0.5
fl3s4 = kp(19)*(0.8 - 0.2) + 0.2;
% ------------transfer fraction intra soil
% cue
mic_cue = kp(20)*(0.8 - 0.1) + 0.1;
% fraction from s1 to s2
fs1s2 = mic_cue;
% fraction of cue that leads to death (doc + soc) 0.5
pdeath2soc = kp(21)*(1 - 0) + 0;
% fraction from enz to doc
fs3s1 = 1;
% fraction from soc to doc
fs4s1 = 1;
% ------------input allocation
% beta to describe the shape of vertical profile 0.95
beta = kp(22)*(0.9999 - 0.5) + 0.5;
% allometric slope for microbial vs enzyme production relationship 1
allo_slope_mic = kp(23)*(1.5 - (-0.5)) + (-0.5);
% maximum and minimum water potential (MPa)
maxpsi= -0.0020;

minpsi= -2;
% parameter for advection (m/yr)
adv = 0;

%-------------------------------------------
% vertical allocation
%-------------------------------------------
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

%-------------------------------------------
% analytical solution for litter pools at steady state
%-------------------------------------------

%-------------env nscalar----------------
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

% Triangle Matrix, A Matrix and K Matrix
sand_vector = mean(sand_vector, 2, 'omitnan');
% Allocation matrix
a_ma = a_matrix(fcwdl2, fl1s1, fl2s1, fl3s4, nan, nan(n_soil_layer, 1), nan(n_soil_layer, 1), nan(n_soil_layer, 1), nan, nan);

kk_ma_middle = nan(npool_vr, npool_vr, month_num);
tri_ma_middle = nan(npool_vr, npool_vr, month_num);
for imonth = 1:month_num
    % decomposition matrix
    monthly_xit = xit(:, imonth);
    monthly_xiw = xiw(:, imonth);
    monthly_xio = xio(:, imonth);
    monthly_xin = xin(:, imonth);
    kk_ma_middle(:, :, imonth) = kk_matrix(monthly_xit, monthly_xiw, monthly_xio, monthly_xin, efolding, ...
        tau4cwd, tau4l1, tau4l2, tau4l3, nan, nan, nan, nan, nan, ...
        nan, nan, nan(npool_vr, 1), nan);
    % tri matrix
    monthly_nbedrock = nbedrock(imonth);
    monthly_altmax_current_profile = altmax_current_profile(imonth);
    monthly_altmax_lastyear_profile = altmax_lastyear_profile(imonth);
    tri_ma_middle(:, :, imonth) = tri_matrix(monthly_nbedrock, monthly_altmax_current_profile, monthly_altmax_lastyear_profile, bio, adv, cryo);
end
tri_ma = mean(tri_ma_middle, 3, 'omitnan');
kk_ma = mean(kk_ma_middle, 3, 'omitnan');

% Analytical Solution
matrix_in=nan(140,1);

% calculate annual carbon input, sum monthly input (gc/m3/month) as annual input (gc/m3/year)
input_tot_cwd = sum(input_vector_cwd, 2, 'omitnan');
input_tot_litter1 = sum(input_vector_litter1, 2, 'omitnan');
input_tot_litter2 = sum(input_vector_litter2, 2, 'omitnan');
input_tot_litter3 = sum(input_vector_litter3, 2, 'omitnan');
matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_year); % litter input gc/m3/day
matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(81:npool_vr,1) = 0;

litter_cpools_initial = (a_ma(:, 1:80)*kk_ma(1:80, 1:80)-tri_ma(:, 1:80))\(-matrix_in(1:end));
litter_cpools_initial(litter_cpools_initial < 0) = 0;

%% set initial values
soc_initial_max = 10^5;
mic_initial_max = 5*10^2;
doc_initial_max = 5*10^2;
enz_initial_max = 1*10^1;

cpool_initial = [litter_cpools_initial; ...
    doc_initial_max*vertical_input; ...
    mic_initial_max*vertical_input; ...
    enz_initial_max*vertical_input; ...
    soc_initial_max*vertical_input]; % 1*ones(npool_vr, 1); % zeros(npool_vr, 1);

%% forward simulation to steady state
% disp([datestr(now), ' simulation started'])

simu_year_size = 3000;

soc_stock = nan(simu_year_size, month_num, 5);

for iyear = 1:simu_year_size
    for imonth = 1:month_num
        if imonth == 1 && iyear == 1
            cpool_current = cpool_initial;
        end
        %-------------------------------------------
        % env. scalar
        %-------------------------------------------
        % temperature related function xit
        % calculate rate constant scalar for soil temperature
        % assuming that the base rate constants are assigned for non-moisture
        % limiting conditions at 25 C.
        xit = nan(n_soil_layer, 1);
        for ilayer = 1 : n_soil_layer
            if soil_temp_profile(ilayer, imonth) >= 0 + kelvin_to_celsius
                xit(ilayer) = q10^((soil_temp_profile(ilayer, imonth) - (kelvin_to_celsius + 25))/10);
            else
                xit(ilayer) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer, imonth) - (0 + kelvin_to_celsius))/10));
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
        
        xiw = soil_water_profile(:, imonth)*w_scaling;
        xiw(xiw > 1) = 1;
        
        %-------------------------------------------
        % K matrix
        %-------------------------------------------
        [kk_ma, p2death] = kk_matrix(xit, xiw, xio(:, imonth), xin(:, imonth), efolding, ...
            tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2_death, tau4s2_enz, tau4s3, tau4s4, ...
            mm_const_assim, mm_const_decom, cpool_current, allo_slope_mic);
        
        %-------------------------------------------
        % A matrix
        %-------------------------------------------
        % fraction from mic to doc
        fs2s1 = p2death*(1-pdeath2soc);
        % fraction from mic to enz
        fs2s3 = 1 - p2death;
        % fraction from mic to soc
        fs2s4 = p2death*pdeath2soc;
        
        a_ma = a_matrix(fcwdl2, fl1s1, fl2s1, fl3s4, fs1s2, fs2s1, fs2s3, fs2s4, fs3s1, fs4s1);
        
        %-------------------------------------------
        % triangle matrix
        %-------------------------------------------
        monthly_nbedrock = nbedrock(imonth);
        monthly_altmax_current_profile = altmax_current_profile(imonth);
        monthly_altmax_lastyear_profile = altmax_lastyear_profile(imonth);
        tri_ma = tri_matrix(monthly_nbedrock, monthly_altmax_current_profile, monthly_altmax_lastyear_profile, bio, adv, cryo);
        
        %-------------------------------------------
        % input vector
        %-------------------------------------------
        % monthly input (gc/m2/month)
        input_tot_cwd = input_vector_cwd(imonth);
        input_tot_litter1 = input_vector_litter1(imonth);
        input_tot_litter2 = input_vector_litter2(imonth);
        input_tot_litter3 = input_vector_litter3(imonth);
        
        matrix_in=nan(npool_vr,1);
        matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_month(imonth)); % litter input gc/m3/day
        matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_month(imonth));
        matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_month(imonth));
        matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_month(imonth));
        matrix_in(81:npool_vr,1) = 0;
        
        
        % decomposition size
        decom_operator = (kk_ma*(cpool_current))*days_per_month(imonth);
        decom_operator((decom_operator - cpool_current) >= 0) = cpool_current((decom_operator - cpool_current) >= 0);
        % calculate carbon transfer by A
        transfer_operator = a_ma*decom_operator + cpool_current;
        % calculate vertical transport
        cpool_next = matrix_in*days_per_month(imonth) - tri_ma*transfer_operator + transfer_operator;
        
        %         cpool_next(cpool_next < 0) = 0;
        %         plot(cpool_next)
        
        cpool_current = cpool_next;
        
        soc_layer = [cpool_current(81:100), cpool_current(101:120), cpool_current(121:140), cpool_current(141:160)];
        soc_layer = [soc_layer, sum(soc_layer, 2)]; % unit gC/m3
        
        soc_stock(iyear, imonth, :) = [sum(soc_layer(:, 1).*dz(1:n_soil_layer)), ...
            sum(soc_layer(:, 2).*dz(1:n_soil_layer)), ...
            sum(soc_layer(:, 3).*dz(1:n_soil_layer)), ...
            sum(soc_layer(:, 4).*dz(1:n_soil_layer)), ...
            sum(soc_layer(:, 5).*dz(1:n_soil_layer))]; % unit gC/m2
        
    end
    if mod(iyear, 20) == 0
        delta_stock = (soc_stock(iyear, month_num, 5) - soc_stock(iyear-19, 1, 5))/20; % unit gC/m2/year
        if delta_stock <= 1
            break
        end
    end
    
    if iyear == simu_year_size
        delta_stock = (soc_stock(iyear, month_num, 5) - soc_stock(iyear-19, 1, 5))/20; % unit gC/m2/year
        if delta_stock > 1
            % soc_stock(:, :, :) = nan;
            soc_stock = -soc_stock;
            soc_layer = -soc_layer;
        end
    end
end

if delta_stock > 1
    ss_index = 0;
    disp([datestr(now), ' simulation finished without steady state'])
else
    ss_index = 1;
    disp([datestr(now), ' simulation finished with steady state'])
end
end

% ----- CENTURY T response function
function catanf_results = catanf(t1)
catanf_results = 11.75 +(29.7 / pi) * atan( pi * 0.031  * ( t1 - 15.4 ));
end



