%clear all;
%close all;
function kk_out = kk_matrix(xit, xiw, xio, xin, efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3)
global use_vertsoilc n_soil_layer days_per_year secspday npool_vr
global zsoi
%use_vertsoilc = 1;
% load soildepth.mat;

% note scalars change with time
n_scalar = xin(1:n_soil_layer, 1);  % nitrogen, fpi
t_scalar = xit;  % temperature
w_scalar = xiw;  % water
o_scalar = xio(1:n_soil_layer, 1);  % oxgen

% turnover rate and time. Copied from SoilBiogeochemDecompCascadeBGCMod
% the belowground parameters from century
% tau_l1 = 1./18.5;
% tau_l2_l3 = 1./4.9;
% tau_s1 = 1./7.3;
% tau_s2 = 1./0.2;
% tau_s3 = 1./.0045;
% century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1
% tau_cwd  = 1./0.3;
% translate to per-day time constant

% tau4cwd = kp(9);               % turnover time of CWD
% tau4l1 = kp(10);               % tau for metabolic litter
% tau4l2 = kp(11);               % tau for cellulose litter
% tau4l3 = kp(12);               % tau for lignin litter
% tau4s1 = kp(13);               % tau for fast SOC
% tau4s2 = kp(14);               % tau for slow SOC
% tau4s3 = kp(15);               % tau for passive SOC

kl1 = 1/(days_per_year * tau4l1);
kl2 = 1/(days_per_year  * tau4l2);
kl3 = 1/(days_per_year  * tau4l3);
ks1 = 1/(days_per_year  * tau4s1);
ks2 = 1/(days_per_year  * tau4s2);
ks3 = 1/(days_per_year  * tau4s3);
kcwd = 1/(days_per_year  * tau4cwd);


decomp_depth_efolding = efolding;

depth_scalar = exp(-zsoi/decomp_depth_efolding);
xi_tw = t_scalar.*w_scalar.*o_scalar;
xi_tw = xi_tw(1:n_soil_layer, 1);

if (use_vertsoilc)
    kk_ma_vr = zeros(npool_vr);   % kk_matrix, decay matrix * scalar matrix
    for j = 1:n_soil_layer
        % CWD exists only on the surface of land
        kk_ma_vr(j,j) = kcwd * xi_tw(j) * depth_scalar(j);
        % other litters and SOC and be decomposed at each layer of the soil
        % profile
        kk_ma_vr(1*n_soil_layer+j, 1*n_soil_layer+j) = kl1 * xi_tw(j) * depth_scalar(j) * n_scalar(j) ;
        kk_ma_vr(2*n_soil_layer+j, 2*n_soil_layer+j) = kl2 * xi_tw(j) * depth_scalar(j) * n_scalar(j);
        kk_ma_vr(3*n_soil_layer+j, 3*n_soil_layer+j) = kl3 * xi_tw(j) * depth_scalar(j) * n_scalar(j);
        kk_ma_vr(4*n_soil_layer+j, 4*n_soil_layer+j) = ks1 * xi_tw(j) * depth_scalar(j);
        kk_ma_vr(5*n_soil_layer+j, 5*n_soil_layer+j) = ks2 * xi_tw(j) * depth_scalar(j);
        kk_ma_vr(6*n_soil_layer+j, 6*n_soil_layer+j) = ks3 * xi_tw(j) * depth_scalar(j);
    end
    kk_out = kk_ma_vr;
else
    kk_ma = zeros(nspools);
    xi_tw(1) = t_scalar(1)*w_scalar(1)*o_scalar(1);
    kk_ma(1,1) = kcwd * xi_tw(1);
    kk_ma(2,2) = kl1 * xi_tw(1)* n_scalar(1);
    kk_ma(3,3) = kl2 * xi_tw(1)* n_scalar(1) ;
    kk_ma(4,4) = kl2 * xi_tw(1)* n_scalar(1) ;
    kk_ma(5,5) = ks1 * xi_tw(1);
    kk_ma(6,6) = ks2 * xi_tw(1);
    kk_ma(7,7) = ks3 * xi_tw(1);
    kk_out = kk_ma;
end
end
