function [Stab_Der_parts Stab_Der] = getClbdot_v2(AC_CONFIGURATION,modelo,alpha,Stab_Der,Geo_tier, Aero,Trim_ITER,Stab_Der_parts)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

b_w       = modelo.ala.b;
Zca_v     = modelo.vertical.Zca;
Xca_v     = modelo.vertical.Xca;
Xcg       = modelo.general.Xcg;

% Cy_betaDot   = Stab_Der.Cybpunto;
Cy_betaDot_VTP = Stab_Der.Cybpunto_VTP;
Cy_betaDot_vee = Stab_Der.Cybpunto_vee;
Cy_betaDot_vee2 = Stab_Der.Cybpunto_vee2;


%% Cl_betaDot
if AC_CONFIGURATION.VTP == 1
    Zca_v       = -Geo_tier.z_VTP_LE; % positive if the wing is below thefuselage centerline.
    Xca_v     = Geo_tier.x_xbar_VTP;
    Cl_betaDot_VTP          = Cy_betaDot_VTP*(Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w;

    if isnan(Cl_betaDot_VTP)
        Cl_betaDot_VTP = 0;
    end
else
    Cl_betaDot_VTP = 0;
end

if AC_CONFIGURATION.Vee == 1
    Zca_v       = -Geo_tier.z_vee_LE; % positive if the wing is below thefuselage centerline.
    Xca_v     = Geo_tier.x_xbar_vee;
    Cl_betaDot_vee          = Cy_betaDot_vee*(Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w;

    if isnan(Cl_betaDot_vee)
        Cl_betaDot_vee = 0;
    end
else
    Cl_betaDot_vee = 0;
end


if AC_CONFIGURATION.Vee2 == 1
    Zca_v2       = -Geo_tier.z_vee2_LE; % positive if the wing is below thefuselage centerline.
    Xca_v2     = Geo_tier.x_xbar_vee2;

    Cl_betaDot_vee2          = Cy_betaDot_vee2*(Zca_v2*cos(alpha) - (Xca_v2 - Xcg)*sin(alpha))/b_w;

    if isnan(Cl_betaDot_vee2)
        Cl_betaDot_vee2 = 0;
    end

else
    Cl_betaDot_vee2 = 0;
end

Cl_betaDot_w1 = 0;
Cl_betaDot_can = 0;
Cl_betaDot_HTP = 0;
Cl_betaDot_nac = 0;

Clbpunto = Cl_betaDot_VTP + Cl_betaDot_vee + Cl_betaDot_vee2 + Cl_betaDot_HTP + Cl_betaDot_nac + Cl_betaDot_can + Cl_betaDot_w1;

Stab_Der.Clbpunto = Clbpunto;
Stab_Der.Clbpunto_VTP = Cl_betaDot_VTP;
Stab_Der.Clbpunto_vee = Cl_betaDot_vee;
Stab_Der.Clbpunto_vee2 = Cl_betaDot_vee2;
Stab_Der.Clbpunto_w1 = Cl_betaDot_w1;
Stab_Der.Clbpunto_can = Cl_betaDot_can;
Stab_Der.Clbpunto_nac = Cl_betaDot_nac;
Stab_Der.Clbpunto_HTP = Cl_betaDot_HTP;

% Stores Derivatives per parts
Stab_Der_parts.Clbpunto = Clbpunto;
Stab_Der_parts.Cl_betaDot_VTP = Cl_betaDot_VTP;
Stab_Der_parts.Cl_betaDot_vee = Cl_betaDot_vee;
Stab_Der_parts.Cl_betaDot_vee2 = Cl_betaDot_vee2;
Stab_Der_parts.Cl_betaDot_w1 = Cl_betaDot_w1;
Stab_Der_parts.Cl_betaDot_can = Cl_betaDot_can;
Stab_Der_parts.Cl_betaDot_nac = Cl_betaDot_nac;
Stab_Der_parts.Cl_betaDot_HTP = Cl_betaDot_HTP;
end