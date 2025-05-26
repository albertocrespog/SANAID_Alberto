function [Stab_Der_parts Stab_Der] = getCybdot_v2(AC_CONFIGURATION,modelo,alpha,Stab_Der, Geo_tier, Aero,Trim_ITER,Stab_Der_parts)


W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

b_w       = modelo.ala.b;
LAMc4_w   = modelo.ala.LAMc4;
TR_w      = modelo.ala.TR;
diedro_w  = modelo.ala.diedro;
Zca_w     = -Geo_tier.z_w1_LE;

%diedro_h  = modelo.horizontal.diedro;
Zca_v     = modelo.vertical.Zca;
Xca_v     = modelo.vertical.Xca;

Minf      = modelo.general.Minf;
Xcg       = modelo.general.Xcg;
Sref      = modelo.general.Sref;

Wmax_fus  = modelo.fuselaje.W;

%% Cy_betaDot
if AC_CONFIGURATION.VTP == 1
    CL_w = Trim_ITER.CL_w1;
    CLa_we      = Stab_Der_parts.CLalpha_w1_e_pw;
    AR_w    = Geo_tier.AR_w1;
    AR_we       =  Geo_tier.AR_w1_e;
    b_w     = Geo_tier.b_w1;
    LAMLE_w = Geo_tier.Lambda_LE_w1;
    LAMc2_w = Geo_tier.Lambda_c2_w1;
    LAMc4_w = Geo_tier.Lambda_c4_w1;
    TR_w    = Geo_tier.lambda_w1;
    diedro_w  = Geo_tier.dihedral_w1;
   
    Zca_v       = -Geo_tier.z_VTP_LE; % positive if the wing is below thefuselage centerline.
    Xca_v     = Geo_tier.x_xbar_VTP;
    CLa_v     = Aero.CL_alpha_VTP_CR;
    S_v       =  Geo_tier.S_VTP;

    twin_VTP = AC_CONFIGURATION.twin_VTP;
    if twin_VTP ==1
        S_v       = 2*S_v;
    else
        S_v       = 1*S_v;
    end

    % DATCOM 7.4.4.4 (2825 PDF)
    sigma_beta_alpha    = sigma_beta_alpha_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, TR_w);
    sigma_beta_diedro   = sigma_beta_diedro_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, Minf);
    sigma_beta_WB       = sign(Zca_w)*sigma_beta_wb_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, TR_w, Wmax_fus/b_w);
    sigma_beta          = sigma_beta_alpha*alpha*180/pi + sigma_beta_diedro*diedro_w + sigma_beta_WB;
    Cy_betaDot_VTP          = 2*CLa_v*sigma_beta*S_v/Sref/b_w*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha));

    if isnan(Cy_betaDot_VTP)
        Cy_betaDot_VTP = 0;
    end
else
    Cy_betaDot_VTP = 0;
end

if AC_CONFIGURATION.Vee == 1
    Zca_v       = -Geo_tier.z_vee_LE; % positive if the wing is below thefuselage centerline.
    Xca_v     = Geo_tier.x_xbar_vee;
    CLa_v = Aero.CYbeta_vee; %% Revisar si es mejor CYbeta_twin
    S_v       =  Geo_tier.S_vee;

    % DATCOM 7.4.4.4 (2825 PDF)
    sigma_beta_alpha    = sigma_beta_alpha_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, TR_w);
    sigma_beta_diedro   = sigma_beta_diedro_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, Minf);
    sigma_beta_WB       = sign(Zca_w)*sigma_beta_wb_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, TR_w, Wmax_fus/b_w);
    sigma_beta          = sigma_beta_alpha*alpha*180/pi + sigma_beta_diedro*diedro_w + sigma_beta_WB;
    Cy_betaDot_vee          = 2*CLa_v*sigma_beta*S_v/Sref/b_w*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha));

    if isnan(Cy_betaDot_vee)
        Cy_betaDot_vee = 0;
    end
else
    Cy_betaDot_vee = 0;
end

if AC_CONFIGURATION.Vee2 == 1
    %diedro_h  = modelo.horizontal.diedro;
    % Zca_v2     = modelo.vertical2.Zca;
    % Xca_v2    = modelo.vertical2.Xca;
    Zca_v2       = -Geo_tier.z_vee2_LE; % positive if the wing is below thefuselage centerline.
    Xca_v2     = Geo_tier.x_xbar_vee2;
    CLa_v2 = Aero.CYbeta_vee2; %% Revisar si es mejor CYbeta_twin
    S_v2       =  Geo_tier.S_vee2;

    % S_v2       = modelo.vertical2.S;
    % S_v2       = modelo.vee2.S;

    % DATCOM 7.4.4.4 (2825 PDF)
    sigma_beta_alpha    = sigma_beta_alpha_calc(2*Zca_v2/b_w, 180/pi*LAMc4_w, TR_w);
    sigma_beta_diedro   = sigma_beta_diedro_calc(2*Zca_v2/b_w, 180/pi*LAMc4_w, Minf);
    sigma_beta_WB       = sign(Zca_w)*sigma_beta_wb_calc(2*Zca_v2/b_w, 180/pi*LAMc4_w, TR_w, Wmax_fus/b_w);
    sigma_beta          = sigma_beta_alpha*alpha*180/pi + sigma_beta_diedro*diedro_w + sigma_beta_WB;
    Cy_betaDot_vee2          = 2*CLa_v2*sigma_beta*S_v2/Sref/b_w*((Xca_v2 - Xcg)*cos(alpha) + Zca_v2*sin(alpha));

    if isnan(Cy_betaDot_vee2)
        Cy_betaDot_vee2 = 0;
    end

else
        Cy_betaDot_vee2 = 0;
end

Cy_betaDot_w1 = 0;
Cy_betaDot_can = 0;
Cy_betaDot_HTP = 0;
Cy_betaDot_nac = 0;

Cybpunto = Cy_betaDot_VTP + Cy_betaDot_vee + Cy_betaDot_vee2 + Cy_betaDot_HTP + Cy_betaDot_nac + Cy_betaDot_can + Cy_betaDot_w1;
Stab_Der.Cybpunto = Cybpunto;
Stab_Der.Cybpunto_VTP = Cy_betaDot_VTP;
Stab_Der.Cybpunto_vee = Cy_betaDot_vee;
Stab_Der.Cybpunto_vee2 = Cy_betaDot_vee2;
Stab_Der.Cybpunto_w1 = Cy_betaDot_w1;
Stab_Der.Cybpunto_can = Cy_betaDot_can;
Stab_Der.Cybpunto_nac = Cy_betaDot_nac;
Stab_Der.Cybpunto_HTP = Cy_betaDot_HTP;

% Stores Derivatives per parts
Stab_Der_parts.Cybpunto = Cybpunto;
Stab_Der_parts.Cy_betaDot_VTP = Cy_betaDot_VTP;
Stab_Der_parts.Cy_betaDot_vee = Cy_betaDot_vee;
Stab_Der_parts.Cy_betaDot_vee2 = Cy_betaDot_vee2;
Stab_Der_parts.Cy_betaDot_w1 = Cy_betaDot_w1;
Stab_Der_parts.Cy_betaDot_can = Cy_betaDot_can;
Stab_Der_parts.Cy_betaDot_nac = Cy_betaDot_nac;
Stab_Der_parts.Cy_betaDot_HTP = Cy_betaDot_HTP;
end
