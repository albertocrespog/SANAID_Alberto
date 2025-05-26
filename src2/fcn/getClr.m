function Stab_Der = getClr(AC_CONFIGURATION,modelo,alpha,Stab_Der, Stab_Der_parts, Trim_ITER)
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

Cy_beta_vert   = Stab_Der_parts.Cy_beta_vert;

beta = sqrt(1 - Minf^(2)*cos(LAMc4_w)^2);
%% Cl_r

% Ala
Num_clr     = 1 + (AR_we*(1 - beta^2))/(2*beta*(AR_we*beta + 2*cos(LAMc4_w))) + ...
            ((AR_we*beta + 2*cos(LAMc4_w))/(AR_we*beta + 4*cos(LAMc4_w)))*...
            ((tan(LAMc4_w))^2)/8;
Den_clr     = 1 + ((AR_we + 2*cos(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)))*...
            ((tan(LAMc4_w))^2)/8;
Clr_CL_M0   = Clr_CL_M0_calc(AR_we, TR_we, LAMc4_w); % Pamadi Fig. 4.28
Clr_CL_CL0  = Num_clr/Den_clr*Clr_CL_M0; % Pamadi Ec. 4.614

DClr_diedro = 1/12*(pi*AR_we*sin(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)); % Pamadi Ec. 4.617

Cl_r_w      = CL_w1*Clr_CL_CL0 + DClr_diedro*diedro_w; % Pamadi Ec. 4.613



% Vertical
% Pamadi 4.618
Cl_r_vert   = -2/b_w^2*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))*(Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))*Cy_beta_vert;

% Twin Vertical Tail configuration
if AC_CONFIGURATION.twin_VTP == 1
    Cl_r_vert = 2*Cl_r_vert;
end

if isnan(Cl_r_vert)
    Cl_r_vert = 0;
end
% DERIVADA TOTAL
Cl_r        = Cl_r_w + Cl_r_vert;

Stab_Der.Clr_w = Cl_r_w;
Stab_Der.Clr_v = Cl_r_vert;
Stab_Der.Clr = Cl_r;
end