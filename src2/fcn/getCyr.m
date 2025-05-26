function Stab_Der = getCyr (AC_CONFIGURATION,modelo,alpha,Stab_Der,Stab_Der_parts)

b_w     = modelo.ala.b;

Xcg     = modelo.general.Xcg;

Xca_v   = modelo.vertical.Xca;
Zca_v   = modelo.vertical.Zca;
Cy_beta_vert = Stab_Der_parts.Cy_beta_vert;
%% Cy_r
Cy_r_w      = 0;
Cy_r_v      = -(2/b_w)*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))*Cy_beta_vert;

% Twin Vertical Tail configuration
if AC_CONFIGURATION.twin_VTP == 1
    Cy_r_v = 2*Cy_r_v;
end

if isnan(Cy_r_v)
    Cy_r_v = 0;
end

Cy_r        = Cy_r_w + Cy_r_v;

Stab_Der.Cyr_w = Cy_r_w;
Stab_Der.Cyr_v = Cy_r_v;
Stab_Der.Cyr = Cy_r;

end