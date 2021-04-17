%clear all;
%close all;
function a_out = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector)
global use_vertsoilc n_soil_layer npool npool_vr

nlevdecomp = n_soil_layer;
nspools = npool;
nspools_vr = npool_vr;

%use_vertsoilc = 1;
% creat diagnal matrics
a_ma_vr = diag(-ones(nspools_vr, 1));
a_ma = diag(-ones(nspools, 1));

% fs2s1 = 1 - fs2s3;
fcwdl3 = 1 - fcwdl2;
% % set path fractions

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
end

if (use_vertsoilc)
    for j = 1:nlevdecomp
        % respiration fraction
        % 1,1,1,fs1s2,fs1s3,fs2s1,fs2s3,1, cwd_fcel,cwd_flig
        % litter 1 to 3 have to go out after respiration, same for S3
%         pathfrac_decomp_cascade =[1, 1, 1, fs1s2, fs1s3, fs2s1, fs2s3, 1, fcwdl2, fcwdl3];
        % transfer fraction  l1s1, l2s1, l3s2, s1s2, s1s3, s2s1, s2s3, s3s1, cwdl2, cwdl3
%         rf_decomp_cascade = [rfl1, rfl2, rfl3, rfs1, rfs1, rfs2, rfs2, rfs3, rfcwd, rfcwd];
        % rf_decomp_cascade = [0.55, 0.5, 0.5, rfs1s2, rfs1s3, 0.55, 0.55, 0.55, 0, 0];
        transfer_fraction = [fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, fcwdl3];
        a_ma_vr((3-1)*nlevdecomp+j,(1-1)*nlevdecomp+j) = transfer_fraction(9);
        a_ma_vr((4-1)*nlevdecomp+j,(1-1)*nlevdecomp+j) = transfer_fraction(10);
        a_ma_vr((5-1)*nlevdecomp+j,(2-1)*nlevdecomp+j) = transfer_fraction(1);
        a_ma_vr((5-1)*nlevdecomp+j,(3-1)*nlevdecomp+j) = transfer_fraction(2);
        a_ma_vr((5-1)*nlevdecomp+j,(6-1)*nlevdecomp+j) = transfer_fraction(6);
        a_ma_vr((5-1)*nlevdecomp+j,(7-1)*nlevdecomp+j) = transfer_fraction(8);
        a_ma_vr((6-1)*nlevdecomp+j,(4-1)*nlevdecomp+j) = transfer_fraction(3);
        if is_default == 1
            a_ma_vr((6-1)*nlevdecomp+j,(5-1)*nlevdecomp+j) = fs1s2_vector(j);
            a_ma_vr((7-1)*nlevdecomp+j,(5-1)*nlevdecomp+j) = fs1s3_vector(j);
            
        else
            a_ma_vr((6-1)*nlevdecomp+j,(5-1)*nlevdecomp+j) = transfer_fraction(4);
            a_ma_vr((7-1)*nlevdecomp+j,(5-1)*nlevdecomp+j) = transfer_fraction(5);
        end
        a_ma_vr((7-1)*nlevdecomp+j,(6-1)*nlevdecomp+j) = transfer_fraction(7);
    end
    a_out = a_ma_vr;
else
    pathfrac_decomp_cascade =[1, 1, 1, fs1s2, fs1s3, fs2s1, fs2s3, 1, fcwdl2, fcwdl3];
    % transfer fraction  l1s1, l2s1, l3s2, s1s2, s1s3, s2s1, s2s3, s3s1, cwdl2, cwdl3
    rf_decomp_cascade = [rfl1, rfl2, rfl3, rfs1, rfs1, rfs2, rfs2, rfs3, rfcwd, rfcwd];
    % rf_decomp_cascade = [0.55, 0.5, 0.5, rfs1s2, rfs1s3, 0.55, 0.55, 0.55, 0, 0];
    a_ma(3,1) = (1.0-rf_decomp_cascade(9))*pathfrac_decomp_cascade(9);
    a_ma(4,1) = (1.0-rf_decomp_cascade(10))*pathfrac_decomp_cascade(10);
    a_ma(5,2) = (1.0-rf_decomp_cascade(1))*pathfrac_decomp_cascade(1);
    a_ma(5,3) = (1.0-rf_decomp_cascade(2))*pathfrac_decomp_cascade(2);
    a_ma(5,6) = (1.0-rf_decomp_cascade(6))*pathfrac_decomp_cascade(6);
    a_ma(5,7) = (1.0-rf_decomp_cascade(8))*pathfrac_decomp_cascade(8);
    a_ma(6,4) = (1.0-rf_decomp_cascade(3))*pathfrac_decomp_cascade(3);
    a_ma(6,5) = (1.0-rf_decomp_cascade(4))*pathfrac_decomp_cascade(4);
    a_ma(7,5) = (1.0-rf_decomp_cascade(5))*pathfrac_decomp_cascade(5);
    a_ma(7,6) = (1.0-rf_decomp_cascade(7))*pathfrac_decomp_cascade(7) ;
    a_out = a_ma;
end %!!!!!!!!!  end of tansfer matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end
