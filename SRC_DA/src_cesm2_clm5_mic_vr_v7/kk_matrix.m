function [kk_out, p2death] = kk_matrix(xit, xiw, xio, xin, efolding, ...
    tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2_death, tau4s2_enz, tau4s3, tau4s4, ...
    mm_const_assim, mm_const_decom, cpool_current, allo_slope_mic)

global use_vertsoilc n_soil_layer days_per_year secspday npool_vr
global zsoi

doc_layer = cpool_current(81:100);
mic_layer = cpool_current(101:120);
enz_layer = cpool_current(121:140);
soc_layer = cpool_current(141:160);

% note scalars change with time
n_scalar = xin(1:n_soil_layer, 1);  % nitrogen, fpi
t_scalar = xit;  % temperature
w_scalar = xiw;  % water
o_scalar = xio(1:n_soil_layer, 1);  % oxgen

% turnover rate and time. Copied from SoilBiogeochemDecompCascadeBGCMod
% the belowground parameters from century
kcwd = 1/(days_per_year * tau4cwd);
kl1 = 1/(days_per_year * tau4l1);
kl2 = 1/(days_per_year * tau4l2);
kl3 = 1/(days_per_year * tau4l3);
ks1 = 1/(days_per_year * tau4s1);


ks2 = 1/(days_per_year * tau4s2_death) + 1/(days_per_year * tau4s2_enz)*mic_layer.^(allo_slope_mic - 1);
ks2(ks2 > 10^10) = 10^10;

ks3 = 1/(days_per_year * tau4s3);
ks4 = 1/(days_per_year * tau4s4);

p2death = 1/(days_per_year * tau4s2_death)./ks2;

decomp_depth_efolding = efolding;

depth_scalar = exp(-zsoi/decomp_depth_efolding);
xi_tw = t_scalar.*w_scalar.*o_scalar;
xi_tw = xi_tw(1:n_soil_layer, 1);


decom_vector_cwd = kcwd .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer);

decom_vector_l1 = kl1 .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer) .* n_scalar(1:n_soil_layer);
decom_vector_l2 = kl2 .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer) .* n_scalar(1:n_soil_layer);
decom_vector_l3 = kl3 .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer) .* n_scalar(1:n_soil_layer);

decom_vector_mic = ks2 .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer);
decom_vector_enz = ks3 .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer);

mm_const_assim_layer = mm_const_assim .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer);
mm_const_decom_layer = mm_const_decom .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer);
decom_vector_doc = ks1 .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer) .* (mic_layer./(mm_const_assim_layer + doc_layer));
decom_vector_soc = ks4 .* xi_tw(1:n_soil_layer) .* depth_scalar(1:n_soil_layer) .* (enz_layer./(mm_const_decom_layer + soc_layer));

% kk_matrix, decay matrix * scalar matrix
kk_ma_vr = diag([decom_vector_cwd; decom_vector_l1; decom_vector_l2; decom_vector_l3; ...
    decom_vector_doc; decom_vector_mic; decom_vector_enz; decom_vector_soc]);

kk_out = kk_ma_vr;
end
