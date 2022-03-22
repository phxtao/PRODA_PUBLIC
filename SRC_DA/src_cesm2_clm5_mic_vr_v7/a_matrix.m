%clear all;
%close all;
function a_out = a_matrix(fcwdl2, fl1s1, fl2s1, fl3s4, fs1s2, fs2s1, fs2s3, fs2s4, fs3s1, fs4s1)
global use_vertsoilc n_soil_layer npool npool_vr

nlevdecomp = n_soil_layer;
nspools = npool;
nspools_vr = npool_vr;

%use_vertsoilc = 1;
% creat diagnal matrics
a_ma_vr = diag(-ones(nspools_vr, 1));

fcwdl3 = 1 - fcwdl2;

if (use_vertsoilc)
    for j = 1:nlevdecomp
        a_ma_vr((3-1)*nlevdecomp+j,(1-1)*nlevdecomp+j) = fcwdl2;
        a_ma_vr((4-1)*nlevdecomp+j,(1-1)*nlevdecomp+j) = fcwdl3;
        a_ma_vr((5-1)*nlevdecomp+j,(2-1)*nlevdecomp+j) = fl1s1;
        a_ma_vr((5-1)*nlevdecomp+j,(3-1)*nlevdecomp+j) = fl2s1;
        a_ma_vr((8-1)*nlevdecomp+j,(4-1)*nlevdecomp+j) = fl3s4;
        a_ma_vr((6-1)*nlevdecomp+j,(5-1)*nlevdecomp+j) = fs1s2;
        a_ma_vr((5-1)*nlevdecomp+j,(6-1)*nlevdecomp+j) = fs2s1(j);
        a_ma_vr((7-1)*nlevdecomp+j,(6-1)*nlevdecomp+j) = fs2s3(j);
        a_ma_vr((8-1)*nlevdecomp+j,(6-1)*nlevdecomp+j) = fs2s4(j);
        a_ma_vr((5-1)*nlevdecomp+j,(7-1)*nlevdecomp+j) = fs3s1;
        a_ma_vr((5-1)*nlevdecomp+j,(8-1)*nlevdecomp+j) = fs4s1;
    end
    a_out = a_ma_vr;
end %!!!!!!!!!  end of tansfer matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end
