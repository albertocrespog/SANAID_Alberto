function [Stab_Der_parts Stab_Der] = get_delta_rv_latdir_deriv_v2(AC_CONFIGURATION,modelo,Stab_Der,Stab_Der_parts,afe,alpha,Aero,Geo_tier,OUTPUT_read_XLSX)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

%% DATOS REQUERIDOS DEL MODELO
y0_b2   = modelo.vertical2.y0_b2;
y1_b2   = modelo.vertical2.y1_b2;
AR      = modelo.vertical2.AR;
TR      = modelo.vertical2.TR;
cf_c    = modelo.vertical2.cm_c;
% t_c     = modelo.vertical2.t_c;
t_c     = Geo_tier.t_c_rudvtr;
Cla     = modelo.vertical2.Cla;
CLa     = modelo.vertical2.CLa;
S_v     = modelo.vertical2.S;

S_w     = Geo_tier.S_w1;
b_w     = modelo.ala.b;

Xca_v    = Geo_tier.x_xbar_vee2;
Zca_v     = Geo_tier.z_zbar_vee2;

Sref    = modelo.general.Sref;
l_fus = modelo.general.l_fus;

prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;
if prop_wash_effect == 1
    rend    = afe.eta_vee2_no_afe;
else
    rend    = 1;
end

D_fus   = modelo.fuselaje.D;
if modelo.vertical2.Xca < modelo.fuselaje.l
    Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical2.Xca);
else
    Dfus_v  = 0;
end

S_v     = modelo.vertical2.S;
b_v     = modelo.vertical2.b;
z_v     = modelo.vertical2.Zca;
CLa_v   = modelo.vertical2.CLa;
lt_v    = modelo.vertical2.Xca;
eta_v   = modelo.vertical2.eta;

S_h     = modelo.horizontal.S;
b_h     = modelo.horizontal.b;
CLa_h   = modelo.horizontal.CLa;
eta_h   = modelo.horizontal.eta;
AR_h    = modelo.horizontal.AR;
Dfus_h  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x,Xca_v, 'pchip');
z_h     = modelo.horizontal.Zca;
TR_h    = modelo.horizontal.TR;
LAMc2_h = modelo.horizontal.LAMc2;
LAMc4_h = modelo.horizontal.LAMc4;
lt_h    = modelo.horizontal.Xca - modelo.horizontal.xca + modelo.horizontal.le_y(end) + modelo.horizontal.ct/2;
diedro_h  = modelo.horizontal.diedro;

S_vee2 = Geo_tier.S_vee2;
CYbeta_vee2 = Aero.CYbeta_vee2;
x_XCG = modelo.general.Xcg;
z_XCG = modelo.general.Zcg;

Mach = modelo.general.Minf;
%% CÁLCULO DE LA DERIVADA DE ESTABILIDAD
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)
k_prima                 = 1;
% The three dimensional ruddervator effectiveness parameter is determined from Figure 8.53 in Airplane Design Part VI 
% and is a function of the V-tail aspect ratio, and ruddervator chord to V-tail chord ratio:
%% Perfil revisado para el VeeTail
airfoil_data =  importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_vee21);
[cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
clalpha_vee2_M0 = 1.05*cla_clatheo*clatheo;
clalpha_vee2 = clalpha_vee2_M0/sqrt(1-Mach^2);
Cla_theo_vee2 = 2*pi + 5.0525*t_c;
Cla_M0 = clalpha_vee2_M0/Cla_theo_vee2;
alphaCL_alphaCl         = alphaCL_alphaCl_calc(cf_c,AR);
% Calculo de Cldelta_Cldeltatheory: Correction Factor for Plain Flap Lift. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.15 (pag 230, 1920 PDF))
Cldelta_Cldeltatheory   = Cldelta_Cldeltatheory_calc(cf_c,Cla_M0);
% Calculo de Cldeltatheory: Lift Effectiveness. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.14 (pag 228, 1918 PDF)
Cldelta_theory          = Cldeltatheory_calc(cf_c,t_c);
% Calculo de Kb: Effect of the taper ratio and flap span. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.4.2, Fig. 8.52 (pag 260, 1950 PDF))
K_b                     = Kb_calc(y0_b2, y1_b2, TR);

% The change in sideslip due to ruddervator directional deflection 
beta_delta_rv = -K_b*Cldelta_Cldeltatheory*Cldelta_theory*(k_prima/clalpha_vee2)*alphaCL_alphaCl;

% The wing-fuselage-"equivalent horizontal tail" interference on the airplane sideforce-coefficient-due-to-sideslip
% derivative of equivalent twin vertical2 tails is obtained from Figure 10.17 in Airplane Design Part VI and is a function of the equivalent 
% twin vertical2 tail span, the fuselage depth at the quarter chord point of the V-tail, the distance between the two equivalent vertical2 tails, 
% and the fuselage length:
Cybetav_Cybetaveff = Cybetav_Cybetaveff_calc(b_h,l_fus,Dfus_v,b_v); 
K_nacareport = get_k_823(Geo_tier.AR_vee2_e,Geo_tier.lambda_vee2);
CL_alpha_wb_vee2_0 = Aero.CL_alpha_vee2_CR;

CYbeta_vee2 = Aero.CYbeta_vee2;
Cy_drv = - Cybetav_Cybetaveff*rend*abs(CYbeta_vee2)*(S_vee2/S_w)*beta_delta_rv*tan(diedro_h);

if isnan(Cy_drv)
    Cy_drv = 0;
end

%% Cálculo preliminar de la derivada de estabilidad C_l_dr.
Cl_drv = Cy_drv*((Zca_v - z_XCG)*cos(alpha)-(Xca_v - x_XCG)*sin(alpha))/b_w;
if isnan(Cl_drv)
    Cl_drv = 0;
end

%% Cálculo preliminar de la derivada de estabilidad C_n_dr.
Cn_drv   = -Cy_drv*((Zca_v - z_XCG)*sin(alpha)+(Xca_v - x_XCG)*cos(alpha))/b_w;

if isnan(Cn_drv)
    Cn_drv = 0;
end

Stab_Der.Cydeltarv2 = Cy_drv;
Stab_Der.Cldeltarv2 = Cl_drv;
Stab_Der.Cndeltarv2 = Cn_drv;

Stab_Der_parts.Cydeltarv2 = Cy_drv;
Stab_Der_parts.Cldeltarv2 = Cl_drv;
Stab_Der_parts.Cndeltarv2 = Cn_drv;


end