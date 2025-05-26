function Stab_Der = getCybdot(modelo,alpha,Stab_Der, AC_CONFIGURATION, Geo_tier, Aero)
b_w       = modelo.ala.b;
LAMc4_w   = modelo.ala.LAMc4;
TR_w      = modelo.ala.TR;
diedro_w  = modelo.ala.diedro;
Zca_w     = -Geo_tier.z_w1_LE;

%diedro_h  = modelo.horizontal.diedro;
Zca_v     = modelo.vertical.Zca;
Xca_v     = modelo.vertical.Xca;

if AC_CONFIGURATION.VTP == 1
    CLa_v     = Aero.CL_alpha_VTP_CR; 
elseif AC_CONFIGURATION.Vee == 1
    CLa_v = Aero.CYbeta_vee; %% Revisar si es mejor CYbeta_twin
end

twin_VTP = AC_CONFIGURATION.twin_VTP;
if twin_VTP ==1
    S_v       = 2*modelo.vertical.S;
else
    S_v       = modelo.vertical.S;
end

Minf      = modelo.general.Minf;
Xcg       = modelo.general.Xcg;
Sref      = modelo.general.Sref;

Wmax_fus  = modelo.fuselaje.W;


%% Cy_betaDot
% DATCOM 7.4.4.4 (2825 PDF)
sigma_beta_alpha    = sigma_beta_alpha_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, TR_w);
sigma_beta_diedro   = sigma_beta_diedro_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, Minf);
sigma_beta_WB       = sign(Zca_w)*sigma_beta_wb_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, TR_w, Wmax_fus/b_w);

sigma_beta          = sigma_beta_alpha*alpha*180/pi + sigma_beta_diedro*diedro_w + sigma_beta_WB;
Cy_betaDot          = 2*CLa_v*sigma_beta*S_v/Sref/b_w*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha)); 

if isnan(Cy_betaDot)
    Cy_betaDot = 0;
end

Stab_Der.Cybpunto = Cy_betaDot;
end
