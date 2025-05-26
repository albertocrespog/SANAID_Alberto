function Stab_Der = getCnr(AC_CONFIGURATION,modelo,trim_alpha,Stab_Der, Stab_Der_parts, Trim_ITER)
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

Cy_beta_vert = Stab_Der_parts.Cy_beta_vert;

% Ala
% pagina 419
Cn_r_w      = -(1/3)*(CD0 + k1*CL_w1 + k2*CL_w1^2);

% Vertical
Cn_r_vert   = (2/b_w^2)*(((Xca_v - Xcg)*cos(trim_alpha) + (Zca_v-Zcg)*sin(trim_alpha))^2)*Cy_beta_vert;

% Twin Vertical Tail configuration
if AC_CONFIGURATION.twin_VTP == 1
    Cn_r_vert = 2*Cn_r_vert;
end

if isnan(Cn_r_vert)
    Cn_r_vert = 0;
end

% DERIVADA TOTAL
Cn_r = Cn_r_vert + Cn_r_w;
Stab_Der.Cnr_w = Cn_r_w;
Stab_Der.Cnr_v = Cn_r_vert;
Stab_Der.Cnr = Cn_r;
end