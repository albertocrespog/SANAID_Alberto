function [Stab_Der_parts Stab_Der] = getCyr (AC_CONFIGURATION,modelo,alpha,Stab_Der,Stab_Der_parts)


W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

b_w     = modelo.ala.b;

Xcg     = modelo.general.Xcg;
% Vee and VTP
Xca_v   = modelo.vertical.Xca;
Zca_v   = modelo.vertical.Zca;
Cy_beta_VTP = Stab_Der_parts.Cy_beta_VTP;
Cy_beta_vee = Stab_Der_parts.Cy_beta_vee;
Cy_beta_vee2 = Stab_Der_parts.Cy_beta_vee2;

%% Cy_r
Cy_r_w1      = 0;
Cy_r_nac      = 0;
Cy_r_can      = 0;
Cy_r_HTP      = 0;


if VTP == 1
    % Vertical
    Cy_r_VTP      = -(2/b_w)*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))*Cy_beta_VTP;

    if isnan(Cy_r_VTP)
        Cy_r_v = 0;
    end
    % Twin Vertical Tail configuration
    if AC_CONFIGURATION.twin_VTP == 1
        Cy_r_VTP = 2*Cy_r_VTP;
    end
else
    Cy_r_VTP = 0;;
end

if Vee == 1
    % Vertical
    Cy_r_vee      = -(2/b_w)*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))*Cy_beta_vee;

    if isnan(Cy_r_vee)
        Cy_r_vee = 0;
    end

else
    Cy_r_vee = 0;
end

if Vee2 == 1
    Zca_v2      = modelo.vertical2.Zca;
    Xca_v2      = modelo.vertical2.Xca;
    % Vertical
    Cy_r_vee2      = -(2/b_w)*((Xca_v2 - Xcg)*cos(alpha) + Zca_v2*sin(alpha))*Cy_beta_vee2;

    if isnan(Cy_r_vee2)
        Cy_r_vee2 = 0;
    end

else
    Cy_r_vee2 = 0;
end

% Total
Cy_r        = Cy_r_w1 + Cy_r_VTP + Cy_r_vee + Cy_r_vee2 + Cy_r_can + Cy_r_nac + Cy_r_HTP;

Stab_Der.Cyr_w1 = Cy_r_w1;
Stab_Der.Cyr_VTP = Cy_r_VTP;
Stab_Der.Cyr_vee = Cy_r_vee;
Stab_Der.Cyr_vee2 = Cy_r_vee2;
Stab_Der.Cyr_nac = Cy_r_nac;
Stab_Der.Cyr_can = Cy_r_can;
Stab_Der.Cyr_HTP = Cy_r_HTP;
Stab_Der.Cyr = Cy_r;

% Stores Derivatives per parts
Stab_Der_parts.Cy_r = Cy_r;
Stab_Der_parts.Cy_r_w1 = Cy_r_w1;
Stab_Der_parts.Cy_r_VTP = Cy_r_VTP;
Stab_Der_parts.Cy_r_vee = Cy_r_vee;
Stab_Der_parts.Cy_r_vee2 = Cy_r_vee2;
Stab_Der_parts.Cy_r_nac = Cy_r_nac;
Stab_Der_parts.Cy_r_can = Cy_r_can;
Stab_Der_parts.Cy_r_HTP = Cy_r_HTP;

end