function [Stab_Der_parts Stab_Der] = getCnr(AC_CONFIGURATION,modelo,trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER,Geo_tier, OUTPUT_read_XLSX)
W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

%% Cn_r
CD0      = modelo.ala.CD0;
k1       = modelo.ala.CD1;
k2       = modelo.ala.CD2;
CL_w1    = Trim_ITER.CL_w1;
b_w      = modelo.ala.b;

Xcg      = modelo.general.Xcg;
Zcg      = modelo.general.Zcg;

Xca_v    = modelo.vertical.Xca;
Zca_v    = modelo.vertical.Zca;
Cy_beta_VTP = Stab_Der_parts.Cy_beta_VTP;
Cy_beta_vee = Stab_Der_parts.Cy_beta_vee;
Cy_beta_vee2 = Stab_Der_parts.Cy_beta_vee2;


%% Cl_r

if W1 == 1
    % Ala
% pagina 419
Cn_r_w1      = -(1/3)*(CD0 + k1*CL_w1 + k2*CL_w1^2);

    if isnan(Cn_r_w1)
        Cn_r_w1 = 0;
    end

else
    Cn_r_w1 = 0;
end


if VTP == 1
    % Vertical
    % Pamadi 4.618
    Cn_r_VTP   = (2/b_w^2)*(((Xca_v - Xcg)*cos(trim_alpha) + (Zca_v-Zcg)*sin(trim_alpha))^2)*Cy_beta_VTP;

    if isnan(Cn_r_VTP)
        Cn_r_VTP = 0;
    end

    % Twin Vertical Tail configuration
    if AC_CONFIGURATION.twin_VTP == 1
        Cn_r_VTP = 2*Cn_r_VTP;
    end

    if isnan(Cn_r_VTP)
        Cn_r_VTP = 0;
    end
else
    Cn_r_VTP = 0;
end

if Vee == 1
    % Vertical
    % Pamadi 4.618
    Cn_r_vee   = (2/b_w^2)*(((Xca_v - Xcg)*cos(trim_alpha) + (Zca_v-Zcg)*sin(trim_alpha))^2)*Cy_beta_vee;

    if isnan(Cn_r_vee)
        Cn_r_vee = 0;
    end
else
    Cn_r_vee = 0;
end

if Vee2 == 1
    Zca_v2      = modelo.vertical2.Zca;
    Xca_v2      = modelo.vertical2.Xca;
    % Vertical
    % Vertical
    % Pamadi 4.618
    Cn_r_vee2   = (2/b_w^2)*(((Xca_v - Xcg)*cos(trim_alpha) + (Zca_v-Zcg)*sin(trim_alpha))^2)*Cy_beta_vee;

    if isnan(Cn_r_vee2)
        Cn_r_vee2 = 0;
    end

else
    Cn_r_vee2 = 0;
end

Cn_r_can = 0;
Cn_r_HTP = 0;
Cn_r_nac = 0;



% DERIVADA TOTAL
Cn_r = Cn_r_w1 + Cn_r_VTP + Cn_r_vee + Cn_r_vee2 + Cn_r_can + Cn_r_HTP + Cn_r_nac;

Stab_Der.Cnr_w1 = Cn_r_w1;
Stab_Der.Cnr_can = Cn_r_can;
Stab_Der.Cnr_vee = Cn_r_vee;
Stab_Der.Cnr_vee2 = Cn_r_vee2;
Stab_Der.Cnr_HTP = Cn_r_HTP;
Stab_Der.Cnr_VTP = Cn_r_VTP;
Stab_Der.Cnr_nac = Cn_r_nac;
Stab_Der.Cnr = Cn_r;

% Stores Derivatives per parts
Stab_Der_parts.Cn_r = Cn_r;
Stab_Der_parts.Cn_r_w1 = Cn_r_w1;
Stab_Der_parts.Cn_r_can = Cn_r_can;
Stab_Der_parts.Cn_r_vee = Cn_r_vee;
Stab_Der_parts.Cn_r_vee2 = Cn_r_vee2;
Stab_Der_parts.Cn_r_HTP = Cn_r_HTP;
Stab_Der_parts.Cn_r_VTP = Cn_r_VTP;
Stab_Der_parts.Cn_r_nac = Cn_r_nac;
end