function [Stab_Der_parts Stab_Der] = getCnbdot_v2(AC_CONFIGURATION,modelo,alpha,Stab_Der,Geo_tier, Aero,Trim_ITER,Stab_Der_parts)

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


%% Cn_betaDot

if AC_CONFIGURATION.VTP == 1
    Zca_v       = -Geo_tier.z_VTP_LE; % positive if the wing is below thefuselage centerline.
    Xca_v     = Geo_tier.x_xbar_VTP;

    Cn_betaDot_VTP          = -Cy_betaDot_VTP*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))/b_w;

    if isnan(Cn_betaDot_VTP)
        Cn_betaDot_VTP = 0;
    end
else
    Cn_betaDot_VTP = 0;
end

if AC_CONFIGURATION.Vee == 1
    Zca_v       = -Geo_tier.z_vee_LE; % positive if the wing is below thefuselage centerline.
    Xca_v     = Geo_tier.x_xbar_vee;

    Cn_betaDot_vee          = -Cy_betaDot_vee*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))/b_w;

    if isnan(Cn_betaDot_vee)
        Cn_betaDot_vee = 0;
    end
else
    Cn_betaDot_vee = 0;
end


if AC_CONFIGURATION.Vee2 == 1
    Zca_v2       = -Geo_tier.z_vee2_LE; % positive if the wing is below thefuselage centerline.
    Xca_v2     = Geo_tier.x_xbar_vee2;

    Cn_betaDot_vee2          = -Cy_betaDot_vee2*((Xca_v2 - Xcg)*cos(alpha) + Zca_v2*sin(alpha))/b_w;

    if isnan(Cn_betaDot_vee2)
        Cn_betaDot_vee2 = 0;
    end

else
    Cn_betaDot_vee2 = 0;
end

Cn_betaDot_w1 = 0;
Cn_betaDot_can = 0;
Cn_betaDot_HTP = 0;
Cn_betaDot_nac = 0;

Cnbpunto = Cn_betaDot_VTP + Cn_betaDot_vee + Cn_betaDot_vee2 + Cn_betaDot_HTP + Cn_betaDot_nac + Cn_betaDot_can + Cn_betaDot_w1;

Stab_Der.Cnbpunto = Cnbpunto;
Stab_Der.Cnbpunto_VTP = Cn_betaDot_VTP;
Stab_Der.Cnbpunto_vee = Cn_betaDot_vee;
Stab_Der.Cnbpunto_vee2 = Cn_betaDot_vee2;
Stab_Der.Cnbpunto_w1 = Cn_betaDot_w1;
Stab_Der.Cnbpunto_can = Cn_betaDot_can;
Stab_Der.Cnbpunto_nac = Cn_betaDot_nac;
Stab_Der.Cnbpunto_HTP = Cn_betaDot_HTP;

% Stores Derivatives per parts
Stab_Der_parts.Cnbpunto = Cnbpunto;
Stab_Der_parts.Cn_betaDot_VTP = Cn_betaDot_VTP;
Stab_Der_parts.Cn_betaDot_vee = Cn_betaDot_vee;
Stab_Der_parts.Cn_betaDot_vee2 = Cn_betaDot_vee2;
Stab_Der_parts.Cn_betaDot_w1 = Cn_betaDot_w1;
Stab_Der_parts.Cn_betaDot_can = Cn_betaDot_can;
Stab_Der_parts.Cn_betaDot_nac = Cn_betaDot_nac;
Stab_Der_parts.Cn_betaDot_HTP = Cn_betaDot_HTP;

end