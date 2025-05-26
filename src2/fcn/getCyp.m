function Stab_Der = getCyp(AC_CONFIGURATION,modelo,alpha,Stab_Der,Stab_Der_parts,Trim_ITER, Geo_tier, OUTPUT_read_XLSX)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

% CL          = modelo.ala.CLa; %% CAMBIAR
CL_w1       = Trim_ITER.CL_w1;
CLa_we      = Stab_Der_parts.CLalpha_w1_e_pw;
AR_w        = modelo.ala.AR;
AR_we       = modelo.ala.ARwe;
eOswald_w   = modelo.ala.oswald;
LAMc4_w     = modelo.ala.LAMc4;
TR_w        = modelo.ala.TR;
% Cla_w       = modelo.ala.Cla;
diedro_w    = modelo.ala.diedro;
Zca_w       = -Geo_tier.z_w1_LE; % positive if the wing is below thefuselage centerline.

b_w         = modelo.ala.b;
Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.ala.Xca, 'pchip');

Zca_v       = modelo.vertical.Zca;
Xca_v       = modelo.vertical.Xca;
S_v         = modelo.vertical.S;
b_v         = modelo.vertical.b;

Sref        = modelo.general.Sref;
Minf        = modelo.general.Minf;
Xcg         = modelo.general.Xcg;
Zcg         = modelo.general.Zcg;

eta_v       = modelo.vertical.eta;

Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca, 'pchip');

beta = sqrt(1 - Minf^2);
B_prandtl = sqrt(1 - (Minf^2)*(cos(LAMc4_w))^2);


Cy_beta_vert = Stab_Der_parts.Cy_beta_vert;
if isnan(Cy_beta_vert)
    Cy_beta_vert = 0;
end

%% CY_p
% Ala
% Se han mezclado los métodos del DATCOM y los de Pamadi (pag 406)
aw1         = CLa_we/(pi*AR_we*eOswald_w);
K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
%K           = ((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - (k1*CLa + 2*k2*CL*CLa))/((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - 2*CL*CLa/pi/ARexp_w); 
Cyp_CL_M0   = Cyp_CL_M0_calc(LAMc4_w*180/pi, TR_w);
Cyp_CL_CL0  = (AR_we + 4*cos(LAMc4_w))/(AR_we*B_prandtl + 4*cos(LAMc4_w))*(AR_we*B_prandtl + cos(LAMc4_w))/(AR_we + cos(LAMc4_w))*Cyp_CL_M0;
t_c = Geo_tier.t_c_flap;
%% Perfil revisado para el Ala
airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_w1);
[cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
clalpha_w2_M0 = 1.05*cla_clatheo*clatheo;
k_Clp       = clalpha_w2_M0/2/pi;   
betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
Clp_die0    = betaClp_k*k_Clp/beta;
DCy_p_die   = (3*sin(diedro_w)*(1 - 4*(Zca_w-Zcg)/b_w*sin(diedro_w)))*Clp_die0;
Cy_p_w      = K*Cyp_CL_CL0*CL_w1 + DCy_p_die; %DATCOM 7.1.2.1
% z is the vertical distance between the c.g. and the wing root
% quarter-chord point, positive for the c.g. above the wing root chord.
% DATCOM 7.1.2.1 (pagina 2523)
%DATCOM  7.4.2.1 pag.332/444

% Vertical
Cy_p_vert   = 2/b_w*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)*Cy_beta_vert;

% Twin Vertical Tail configuration
if AC_CONFIGURATION.twin_VTP == 1
    Cy_p_vert = 2*Cy_p_vert;
end

if isnan(Cy_p_vert)
    Cy_p_vert = 0;
end

%% DERIVADA TOTAL
Cy_p = Cy_p_vert + Cy_p_w;
Stab_Der.Cyp_w_diedro = DCy_p_die;
Stab_Der.Cyp_w = Cy_p_w;
Stab_Der.Cyp_v = Cy_p_vert;
Stab_Der.Cyp = Cy_p;

end