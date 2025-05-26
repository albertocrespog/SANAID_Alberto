function Stab_Der = getClp(AC_CONFIGURATION,modelo,alpha,Stab_Der,Stab_Der_parts, Geo_tier, OUTPUT_read_XLSX)

Zca_v     = modelo.vertical.Zca;
Xca_v     = modelo.vertical.Xca;

b_w       = modelo.ala.b;
Zca_w     = -Geo_tier.z_w1_LE;
diedro_w  = modelo.ala.diedro;
TR_w      = modelo.ala.TR;
AR_w      = modelo.ala.AR;
LAMc4_w   = modelo.ala.LAMc4;

Xcg       = modelo.general.Xcg;
Minf      = modelo.general.Minf;

% Cy_beta_vert   = Stab_Der.Cyb_v;
% Cy_p_vert      = Stab_Der.Cyp_v;

beta = sqrt(1 - Minf^2);
%% Cl_p
% 
Cy_beta_vert = Stab_Der_parts.Cy_beta_vert;
% Pamadi Ecuacion 4.579, Pagina 412

Cl_p_vert           = abs(2*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w)*(((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w))*Cy_beta_vert;

if isnan(Cl_p_vert)
    Cl_p_vert = 0;
end

% Twin Vertical Tail configuration
if AC_CONFIGURATION.twin_VTP == 1
    Cl_p_vert = 2*Cl_p_vert;
end

% Ala
% Pamadi Ecuacion 4.576, Pagina 412
% Effect of drag on rolling moment is ignored

t_c = Geo_tier.t_c_flap;
%% Perfil revisado para el Ala
airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
[cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
k_Clp       = clalpha_w2_M0/2/pi;   
betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
Clp_ClpDiedro0  = 1 - 4*Zca_w/b_w*sin(diedro_w) + 3*(2*Zca_w/b_w)^2*(sin(diedro_w))^2;
Cl_p_w          = betaClp_k*k_Clp/beta*Clp_ClpDiedro0;

% DERIVADA TOTAL

Cl_p            = Cl_p_vert + Cl_p_w;

Stab_Der.Clp_w = Cl_p_w;
Stab_Der.Clp_v = Cl_p_vert;
Stab_Der.Clp = Cl_p;
end