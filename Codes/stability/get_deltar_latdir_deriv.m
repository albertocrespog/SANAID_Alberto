function Stab_Der = get_deltar_latdir_deriv(AC_CONFIGURATION,modelo,Stab_Der,Geo_tier,afe,alpha,OUTPUT_read_XLSX, Aero)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

%% DATOS REQUERIDOS DEL MODELO
y0_b2   = modelo.vertical.y0_b2;
y1_b2   = modelo.vertical.y1_b2;
AR      = modelo.vertical.AR;
TR      = modelo.vertical.TR;
cf_c    = modelo.vertical.cm_c;
t_c     = Geo_tier.t_c_rudder;
% Cla     = modelo.vertical.Cla;
% CLa   = modelo.vertical.CLa;
CLa     = Aero.CL_alpha_VTP_CR; %el de aero, sin multiplicar por eta*S/Sref
S_v     = modelo.vertical.S;
Xca_v   = Geo_tier.x_xbar_VTP;
Zca_v   = Geo_tier.z_zbar_VTP;
b       = Geo_tier.b_w1;

prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;

Sref    = modelo.general.Sref;
% %Added 
if prop_wash_effect == 1
    rend    = afe.eta_VTP_no_afe;
else
    rend=1;
end

x_XCG = modelo.general.Xcg;
z_XCG = modelo.general.Zcg;

Mach = modelo.general.Minf;


%% CÁLCULO DE LA DERIVADA DE ESTABILIDAD
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)
%% CÁLCULO DE LA DERIVADA DE ESTABILIDAD
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)
k_prima                 = 1;

%% Perfil revisado para el VTP
airfoil_data = importdata(OUTPUT_read_XLSX.Aerodynamic_Data_flags.airfoil_VTP1);
[cla_clatheo,clatheo] = get_cla_clatheo_airfoil(airfoil_data,t_c);
clalpha_VTP_M0 = 1.05*cla_clatheo*clatheo;
clalpha_VTP = clalpha_VTP_M0/sqrt(1-Mach^2);
Cla_theo_VTP = 2*pi + 5.0525*t_c;
Cla_M0 = clalpha_VTP_M0/Cla_theo_VTP;
% The three dimensional ruddervator effectiveness parameter is determined from Figure 8.53 in Airplane Design Part VI 
% and is a function of the V-tail aspect ratio, and ruddervator chord to V-tail chord ratio:
alphaCL_alphaCl         = alphaCL_alphaCl_calc(cf_c,AR);
% Calculo de Cldelta_Cldeltatheory: Correction Factor for Plain Flap Lift. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.15 (pag 230, 1920 PDF))
Cldelta_Cldeltatheory   = Cldelta_Cldeltatheory_calc(cf_c,Cla_M0);
% Calculo de Cldeltatheory: Lift Effectiveness. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.14 (pag 228, 1918 PDF)
Cldelta_theory          = Cldeltatheory_calc(cf_c,t_c);
% Calculo de Kb: Effect of the taper ratio and flap span. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.4.2, Fig. 8.52 (pag 260, 1950 PDF))
K_b                     = Kb_calc(y0_b2, y1_b2, TR);

% The change in sideslip due to ruddervator directional deflection 
beta_delta_r = -K_b*Cldelta_Cldeltatheory*Cldelta_theory*(k_prima/clalpha_VTP)*alphaCL_alphaCl;


Cy_dr   = -rend*CLa*(S_v/Sref)*beta_delta_r;
if AC_CONFIGURATION.twin_VTP == 1
    Cy_dr = 2*Cy_dr;
end

% Cy_dr   = K_b*CLa*S_v/Sref*rend*k_prima/Cla*alphaCL_alphaCl*Cldelta_Cldeltatheory*Cldelta_theory

if isnan(Cy_dr)
    Cy_dr = 0;
end
Stab_Der.Cydeltar = Cy_dr;

%% Cálculo preliminar de la derivada de estabilidad C_l_dr.
Cl_dr   = ((Zca_v - z_XCG)*cos(alpha)-(Xca_v - x_XCG)*sin(alpha))/b*Cy_dr;
Stab_Der.Cldeltar = Cl_dr;

%% CÁLCULO DE LA DERIVADA DE CONTROL
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 462, 2152 PDF)

Cn_dr   = -((Zca_v - z_XCG)*sin(alpha)+ (Xca_v - x_XCG)*cos(alpha))/b*Cy_dr;
Stab_Der.Cndeltar = Cn_dr;


end