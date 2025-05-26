function [Stab_Der_parts Stab_Der] = getClr_v2(AC_CONFIGURATION,modelo,alpha,Stab_Der, Stab_Der_parts, Trim_ITER, Geo_tier, OUTPUT_read_XLSX)


W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

AR_we     = modelo.ala.ARwe;
LAMc4_w   = modelo.ala.LAMc4;
TR_we     = modelo.ala.TR_we;
diedro_w  = modelo.ala.diedro;
CL_w1     = Trim_ITER.CL_w1;
b_w       = modelo.ala.b;

Minf      = modelo.general.Minf;
Xcg       = modelo.general.Xcg;

Xca_v     = modelo.vertical.Xca;
Zca_v     = modelo.vertical.Zca;

Cy_beta_VTP = Stab_Der_parts.Cy_beta_VTP;
Cy_beta_vee = Stab_Der_parts.Cy_beta_vee;
Cy_beta_vee2 = Stab_Der_parts.Cy_beta_vee2;

beta = sqrt(1 - Minf^(2)*cos(LAMc4_w)^2);
%% Cl_r

%% Ala
if W1 == 1

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
    Zca_w       = -Geo_tier.z_w1_LE; % positive if the wing is below thefuselage centerline.

    LAMc4_w = Geo_tier.Lambda_c4_w1;
    AR_w    = Geo_tier.AR_w1;
    TR_w    = Geo_tier.lambda_w1;
    TR_we    = Geo_tier.lambda_w1_e;
    b_w     = Geo_tier.b_w1;
    diedro_w  = Geo_tier.dihedral_w1;
    Zca_w       = -Geo_tier.z_w1_LE; % positive if the wing is below thefuselage centerline.
    t_c = Geo_tier.t_c_flap;

    Num_clr     = 1 + (AR_we*(1 - beta^2))/(2*beta*(AR_we*beta + 2*cos(LAMc4_w))) + ...
        ((AR_we*beta + 2*cos(LAMc4_w))/(AR_we*beta + 4*cos(LAMc4_w)))*...
        ((tan(LAMc4_w))^2)/8;
    Den_clr     = 1 + ((AR_we + 2*cos(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)))*...
        ((tan(LAMc4_w))^2)/8;
    Clr_CL_M0   = Clr_CL_M0_calc(AR_we, TR_we, LAMc4_w); % Pamadi Fig. 4.28
    Clr_CL_CL0  = Num_clr/Den_clr*Clr_CL_M0; % Pamadi Ec. 4.614

    DClr_diedro = 1/12*(pi*AR_we*sin(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)); % Pamadi Ec. 4.617

    Cl_r_w1      = CL_w*Clr_CL_CL0 + DClr_diedro*diedro_w; % Pamadi Ec. 4.613
    if isnan(Cl_r_w1)
        Cl_r_w1 = 0;
    end

else
    Cl_r_w1 = 0;
end


%% Can
if Can == 1

    CL_w = Trim_ITER.CL_can;
    CLa_we      = Stab_Der_parts.CLalpha_can_e_pw;
    AR_w    = Geo_tier.AR_can;
    AR_we       =  Geo_tier.AR_can_e;
    b_w     = Geo_tier.b_can;
    LAMLE_w = Geo_tier.Lambda_LE_can;
    LAMc2_w = Geo_tier.Lambda_c2_can;
    LAMc4_w = Geo_tier.Lambda_c4_can;
    TR_w    = Geo_tier.lambda_can;
    diedro_w  = Geo_tier.dihedral_can;
    Zca_w       = -Geo_tier.z_can_LE; % positive if the wing is below thefuselage centerline.

    LAMc4_w = Geo_tier.Lambda_c4_can;
    AR_w    = Geo_tier.AR_can;
    TR_w    = Geo_tier.lambda_can;
    TR_we    = Geo_tier.lambda_can_e;
    b_w     = Geo_tier.b_can;
    diedro_w  = Geo_tier.dihedral_can;
    Zca_w       = -Geo_tier.z_can_LE; % positive if the wing is below thefuselage centerline.
    t_c = Geo_tier.t_c_flap;

    Num_clr     = 1 + (AR_we*(1 - beta^2))/(2*beta*(AR_we*beta + 2*cos(LAMc4_w))) + ...
        ((AR_we*beta + 2*cos(LAMc4_w))/(AR_we*beta + 4*cos(LAMc4_w)))*...
        ((tan(LAMc4_w))^2)/8;
    Den_clr     = 1 + ((AR_we + 2*cos(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)))*...
        ((tan(LAMc4_w))^2)/8;
    Clr_CL_M0   = Clr_CL_M0_calc(AR_we, TR_we, LAMc4_w); % Pamadi Fig. 4.28
    Clr_CL_CL0  = Num_clr/Den_clr*Clr_CL_M0; % Pamadi Ec. 4.614

    DClr_diedro = 1/12*(pi*AR_we*sin(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)); % Pamadi Ec. 4.617

    Cl_r_can      = CL_w*Clr_CL_CL0 + DClr_diedro*diedro_w; % Pamadi Ec. 4.613
    if isnan(Cl_r_can)
        Cl_r_can = 0;
    end

else
    Cl_r_can = 0;
end


%% HTP
if HTP == 1

    CL_w = Trim_ITER.CL_HTP;
    CLa_we      = Stab_Der_parts.CLalpha_HTP_e_pw;
    AR_w    = Geo_tier.AR_HTP;
    AR_we       =  Geo_tier.AR_HTP_e;
    b_w     = Geo_tier.b_HTP;
    LAMLE_w = Geo_tier.Lambda_LE_HTP;
    LAMc2_w = Geo_tier.Lambda_c2_HTP;
    LAMc4_w = Geo_tier.Lambda_c4_HTP;
    TR_w    = Geo_tier.lambda_HTP;
    diedro_w  = Geo_tier.dihedral_HTP;
    Zca_w       = -Geo_tier.z_HTP_LE; % positive if the wing is below thefuselage centerline.

    LAMc4_w = Geo_tier.Lambda_c4_HTP;
    AR_w    = Geo_tier.AR_HTP;
    TR_w    = Geo_tier.lambda_HTP;
    TR_we    = Geo_tier.lambda_HTP_e;
    b_w     = Geo_tier.b_HTP;
    diedro_w  = Geo_tier.dihedral_HTP;
    Zca_w       = -Geo_tier.z_HTP_LE; % positive if the wing is below thefuselage centerline.
    t_c = Geo_tier.t_c_flap;

    Num_clr     = 1 + (AR_we*(1 - beta^2))/(2*beta*(AR_we*beta + 2*cos(LAMc4_w))) + ...
        ((AR_we*beta + 2*cos(LAMc4_w))/(AR_we*beta + 4*cos(LAMc4_w)))*...
        ((tan(LAMc4_w))^2)/8;
    Den_clr     = 1 + ((AR_we + 2*cos(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)))*...
        ((tan(LAMc4_w))^2)/8;
    Clr_CL_M0   = Clr_CL_M0_calc(AR_we, TR_we, LAMc4_w); % Pamadi Fig. 4.28
    Clr_CL_CL0  = Num_clr/Den_clr*Clr_CL_M0; % Pamadi Ec. 4.614

    DClr_diedro = 1/12*(pi*AR_we*sin(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)); % Pamadi Ec. 4.617

    Cl_r_HTP      = CL_w*Clr_CL_CL0 + DClr_diedro*diedro_w; % Pamadi Ec. 4.613
    if isnan(Cl_r_HTP)
        Cl_r_HTP = 0;
    end

else
    Cl_r_HTP = 0;
end


if VTP == 1
    % Vertical
    % Pamadi 4.618
    Cl_r_VTP   = -2/b_w^2*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))*(Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))*Cy_beta_VTP;

    if isnan(Cl_r_VTP)
        Cl_r_VTP = 0;
    end

    % Twin Vertical Tail configuration
    if AC_CONFIGURATION.twin_VTP == 1
        Cl_r_VTP = 2*Cl_r_VTP;
    end
else
    Cl_r_VTP = 0;
end

if Vee == 1
    % Vertical
    % Pamadi 4.618
    Cl_r_vee   = -2/b_w^2*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))*(Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))*Cy_beta_vee;

    if isnan(Cl_r_vee)
        Cl_r_vee = 0;
    end
else
    Cl_r_vee = 0;
end

if Vee2 == 1
    Zca_v2      = modelo.vertical2.Zca;
    Xca_v2      = modelo.vertical2.Xca;
    % Vertical
    % Vertical
    % Pamadi 4.618
    Cl_r_vee2   = -2/b_w^2*((Xca_v2 - Xcg)*cos(alpha) + Zca_v2*sin(alpha))*(Zca_v2*cos(alpha) - (Xca_v2 - Xcg)*sin(alpha))*Cy_beta_vee2;

    if isnan(Cl_r_vee2)
        Cl_r_vee2 = 0;
    end

else
    Cl_r_vee2 = 0;
end


if Nac == 1
    Cl_r_nac = 0;
else
    Cl_r_nac = 0;
end

% DERIVADA TOTAL
Cl_r        = Cl_r_w1 + Cl_r_VTP + Cl_r_vee + Cl_r_vee2 + Cl_r_HTP + Cl_r_can + Cl_r_nac;

Stab_Der.Clr_w1 = Cl_r_w1;
Stab_Der.Clr_VTP = Cl_r_VTP;
Stab_Der.Clr_vee = Cl_r_vee;
Stab_Der.Clr_vee2 = Cl_r_vee2;
Stab_Der.Clr_HTP = Cl_r_HTP;
Stab_Der.Clr_can = Cl_r_can;
Stab_Der.Clr_nac = Cl_r_nac;
Stab_Der.Clr = Cl_r;

% Stores Derivatives per parts
Stab_Der_parts.Cl_r = Cl_r;
Stab_Der_parts.Cl_r_w1 = Cl_r_w1;
Stab_Der_parts.Cl_r_VTP = Cl_r_VTP;
Stab_Der_parts.Cl_r_vee = Cl_r_vee;
Stab_Der_parts.Cl_r_vee2 = Cl_r_vee2;
Stab_Der_parts.Cl_r_nac = Cl_r_nac;
Stab_Der_parts.Cl_r_can = Cl_r_can;
Stab_Der_parts.Cl_r_HTP = Cl_r_HTP;
end