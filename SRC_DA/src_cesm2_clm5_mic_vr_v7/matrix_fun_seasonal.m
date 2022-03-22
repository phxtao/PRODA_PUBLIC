function [cpool_current, soc_stock_final, soc_layer, ss_index] = matrix_fun_seasonal(kp, wosis_layer_obs, nbedrock, sand_vector, npp_mean, ...
    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin)
warning off

%-------------------------------------------
% define global parameters
%-------------------------------------------
use_beta = 1;
month_num = 12;
season_num = 4;

normalize_q10_to_century_tfunc = false;
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]';
days_per_season = [90, 91, 92, 92]';

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
% convert monthly input into seasonal
%-------------------------------------------
nbedrock = mean(reshape(nbedrock, [month_num/season_num, season_num]), 1);
npp_mean = mean(reshape(npp_mean, [month_num/season_num, season_num]), 1);
sand_vector = reshape(mean(reshape(sand_vector, [n_soil_layer, 3, 4]), [2]), [20, 4]);
input_vector_cwd = sum(reshape(input_vector_cwd, [3, 4]), 1);
input_vector_litter1 = sum(reshape(input_vector_litter1, [3, 4]), 1);
input_vector_litter2 = sum(reshape(input_vector_litter2, [3, 4]), 1);
input_vector_litter3 = sum(reshape(input_vector_litter3, [3, 4]), 1);
altmax_current_profile = mean(reshape(altmax_current_profile, [3, 4]), 1);
altmax_lastyear_profile = mean(reshape(altmax_lastyear_profile, [3, 4]), 1);
soil_temp_profile = reshape(mean(reshape(soil_temp_profile, [n_soil_layer, 3, 4]), 2), [n_soil_layer, 4]);
soil_water_profile = reshape(mean(reshape(soil_water_profile, [n_soil_layer, 3, 4]), 2), [n_soil_layer, 4]);
xio = reshape(mean(reshape(xio, [n_soil_layer, 3, 4]), 2), [n_soil_layer, 4]);
xin = reshape(mean(reshape(xin, [n_soil_layer, 3, 4]), 2), [n_soil_layer, 4]);
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
% tau for DOC (yr) 1.1*10-2
tau4s1 = 10^(kp(9)*((0) - (-3)) + (-3));
% tau for MIC (yr) 0.57 for death
tau4s2_death = 10^(kp(10)*(0 - (-2)) + (-2));
% tau for MIC (yr) 22 for enz production
tau4s2_enz = 10^(kp(11)*(1.5 - (0)) + (0)); % kp(11)*(25 - 0.5) + 0.5;
% tau for ENZ (yr) 0.11
tau4s3 = kp(12)*(1.5 - 0) + 0;
% tau for SOC (yr)  default: 4.6*10^-5 -- 1.1*10^-4
tau4s4 = 10^(kp(13)*((-2) - (-7)) + (-7));
% ------------michaelis-menten
% (michaelis-menten) concentration DOC (yr) for half max assimlation reaction from DOC to MIC, unit: gc/m3, 4*10^2
mm_const_assim = 10^(kp(14)*(4 - 0) + 0);
% (michaelis-menten) concentration DOC (yr) for half max decomposition reaction from SOC to MIC, unit: gc.m3, 5*10^(4) -- 6*10^5
mm_const_decom = 10^(kp(15)*(9 - 4) + 4);
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
mic_cue = kp(20)*(0.7 - 0.01) + 0.01;
% fraction from s1 to s2
fs1s2 = mic_cue;
% fraction of cue that leads to death (doc + soc) 0.5
pdeath2soc = kp(21)*(0.9 - 0.1) + 0.1;
% fraction from enz to doc
fs3s1 = 1;
% fraction from soc to doc
fs4s1 = 1;
% ------------input allocation
% beta to describe the shape of vertical profile 0.95
beta = kp(22)*(0.9999 - 0.5) + 0.5;
% allometric slope for microbial vs enzyme production relationship 1
allo_slope_mic = kp(23)*(1.5 - (0)) + (0);
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
xit = nan(n_soil_layer, season_num);
xiw = nan(n_soil_layer, season_num);

for iseason = 1:season_num
    % temperature related function xit
    % calculate rate constant scalar for soil temperature
    % assuming that the base rate constants are assigned for non-moisture
    % limiting conditions at 25 C.
    for ilayer = 1 : n_soil_layer
        if soil_temp_profile(ilayer, iseason) >= 0 + kelvin_to_celsius
            xit(ilayer, iseason) = q10^((soil_temp_profile(ilayer, iseason) - (kelvin_to_celsius + 25))/10);
        else
            xit(ilayer, iseason) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer, iseason) - (0 + kelvin_to_celsius))/10));
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

kk_ma_middle = nan(npool_vr, npool_vr, season_num);
tri_ma_middle = nan(npool_vr, npool_vr, season_num);
for iseason = 1:season_num
    % decomposition matrix
    seasonal_xit = xit(:, iseason);
    seasonal_xiw = xiw(:, iseason);
    seasonal_xio = xio(:, iseason);
    seasonal_xin = xin(:, iseason);
    kk_ma_middle(:, :, iseason) = kk_matrix(seasonal_xit, seasonal_xiw, seasonal_xio, seasonal_xin, efolding, ...
        tau4cwd, tau4l1, tau4l2, tau4l3, nan, nan, nan, nan, nan, ...
        nan, nan, nan(npool_vr, 1), nan);
    % tri matrix
    seasonal_nbedrock = nbedrock(iseason);
    seasonal_altmax_current_profile = altmax_current_profile(iseason);
    seasonal_altmax_lastyear_profile = altmax_lastyear_profile(iseason);
    tri_ma_middle(:, :, iseason) = tri_matrix(seasonal_nbedrock, seasonal_altmax_current_profile, seasonal_altmax_lastyear_profile, bio, adv, cryo);
end
tri_ma = mean(tri_ma_middle, 3, 'omitnan');
kk_ma = mean(kk_ma_middle, 3, 'omitnan');

% Analytical Solution
matrix_in=nan(140,1);

% calculate annual carbon input, sum seasonal input (gc/m3/season) as annual input (gc/m3/year)
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
soc_initial_max = max(wosis_layer_obs, [], 'omitnan'); % 10^5;
mic_initial_max = 0.005*soc_initial_max; % 5*10^2;
doc_initial_max = 0.005*soc_initial_max; % 5*10^2;
enz_initial_max = 0.02*mic_initial_max; % 1*10^1;

cpool_initial = [litter_cpools_initial; ...
    doc_initial_max*vertical_input; ...
    mic_initial_max*vertical_input; ...
    enz_initial_max*vertical_input; ...
    soc_initial_max*vertical_input]; % 1*ones(npool_vr, 1); % zeros(npool_vr, 1);

%% forward simulation to steady state
% disp([datestr(now), ' simulation started'])

simu_year_size = 10000;

soc_stock = nan(simu_year_size, season_num, 5);

for iyear = 1:simu_year_size
    for iseason = 1:season_num
        if iseason == 1 && iyear == 1
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
            if soil_temp_profile(ilayer, iseason) >= 0 + kelvin_to_celsius
                xit(ilayer) = q10^((soil_temp_profile(ilayer, iseason) - (kelvin_to_celsius + 25))/10);
            else
                xit(ilayer) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer, iseason) - (0 + kelvin_to_celsius))/10));
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
        
        xiw = soil_water_profile(:, iseason)*w_scaling;
        xiw(xiw > 1) = 1;
        
        %-------------------------------------------
        % K matrix
        %-------------------------------------------
        [kk_ma, p2death] = kk_matrix(xit, xiw, xio(:, iseason), xin(:, iseason), efolding, ...
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
        seasonal_nbedrock = nbedrock(iseason);
        seasonal_altmax_current_profile = altmax_current_profile(iseason);
        seasonal_altmax_lastyear_profile = altmax_lastyear_profile(iseason);
        tri_ma = tri_matrix(seasonal_nbedrock, seasonal_altmax_current_profile, seasonal_altmax_lastyear_profile, bio, adv, cryo);
        
        %-------------------------------------------
        % input vector
        %-------------------------------------------
        % seasonal input (gc/m2/season)
        input_tot_cwd = input_vector_cwd(iseason);
        input_tot_litter1 = input_vector_litter1(iseason);
        input_tot_litter2 = input_vector_litter2(iseason);
        input_tot_litter3 = input_vector_litter3(iseason);
        
        matrix_in=nan(npool_vr,1);
        matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason)); % litter input gc/m3/day
        matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
        matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
        matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
        matrix_in(81:npool_vr,1) = 0;
        
        
        % decomposition size
        decom_operator = (kk_ma*(cpool_current))*days_per_season(iseason);
        decom_operator((decom_operator - cpool_current) >= 0) = cpool_current((decom_operator - cpool_current) >= 0);
        % calculate carbon transfer by A
        transfer_operator = a_ma*decom_operator + cpool_current;
        % calculate vertical transport
        cpool_next = matrix_in*days_per_season(iseason) - tri_ma*transfer_operator + transfer_operator;
        
        %         cpool_next(cpool_next < 0) = 0;
        %         plot(cpool_next)
        
        cpool_current = cpool_next;
        
        soc_layer = [cpool_current(81:100), cpool_current(101:120), cpool_current(121:140), cpool_current(141:160)];
        soc_layer = [soc_layer, sum(soc_layer, 2)]; % unit gC/m3
        
        soc_stock(iyear, iseason, :) = [sum(soc_layer(:, 1).*dz(1:n_soil_layer)), ...
            sum(soc_layer(:, 2).*dz(1:n_soil_layer)), ...
            sum(soc_layer(:, 3).*dz(1:n_soil_layer)), ...
            sum(soc_layer(:, 4).*dz(1:n_soil_layer)), ...
            sum(soc_layer(:, 5).*dz(1:n_soil_layer))]; % unit gC/m2
        
    end
	% constraints to model simulation
    if mod(iyear, 20) == 0
        delta_stock = (soc_stock(iyear, season_num, 5) - soc_stock(iyear-19, 1, 5))/20; % unit gC/m2/year
        if delta_stock <= 1
            soil_15cm_loc = 4;
			if soc_layer(soil_15cm_loc, 4)/soc_layer(soil_15cm_loc, 2) >= 15 && ...
				soc_layer(soil_15cm_loc, 4)/soc_layer(soil_15cm_loc, 1) >= 100 && ...
				soc_stock(iyear, season_num, 4)/soc_stock(iyear, season_num, 2) >= 50 && ... % soc >> mic (total)
				soc_stock(iyear, season_num, 4)/soc_stock(iyear, season_num, 1) >= 10 && ... % soc > doc	
				soc_stock(iyear, season_num, 2)/soc_stock(iyear, season_num, 3) >= 1 % mic > enz
			% if soc_stock(iyear, season_num, 4)/soc_stock(iyear, season_num, 2) >= 100 && ... % soc >> mic
			% 		soc_stock(iyear, season_num, 4)/soc_stock(iyear, season_num, 1) >= 1 && ... % soc > doc
			% 		soc_stock(iyear, season_num, 2)/soc_stock(iyear, season_num, 3) >= 1 % mic > enz
				
				judge_pool_size = 1;
			else
				judge_pool_size = 0;
			end

			break
        end
    end
    
    if iyear == simu_year_size
        delta_stock = (soc_stock(iyear, season_num, 5) - soc_stock(iyear-19, 1, 5))/20; % unit gC/m2/year
    end
end

if delta_stock > 1 || judge_pool_size == 0
    ss_index = 0;
    soc_stock = -soc_stock;
	soc_layer = -soc_layer;
	% disp([datestr(now), ' simulation finished without steady state'])
else
    ss_index = 1;
    % disp([datestr(now), ' simulation finished with steady state'])
end

soc_stock_final = reshape(soc_stock(iyear, season_num, :), [5, 1]);
end

% ----- CENTURY T response function
function catanf_results = catanf(t1)
catanf_results = 11.75 +(29.7 / pi) * atan( pi * 0.031  * ( t1 - 15.4 ));
end



