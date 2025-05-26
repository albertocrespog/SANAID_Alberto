% Cálculo preliminar de la derivada de estabilidad C_y_dr.
function Stab_Der = getCldrv(modelo,alpha,Stab_Der)
%% DATOS REQUERIDOS DEL MODELO
y0_b2   = modelo.vertical.y0_b2;
y1_b2   = modelo.vertical.y1_b2;
AR      = modelo.vertical.AR;
TR      = modelo.vertical.TR;
cf_c    = modelo.vertical.cm_c;
t_c     = modelo.vertical.t_c;
Cla     = modelo.vertical.Cla;
CLa     = modelo.vertical.CLa;
S_v     = modelo.vertical.S;

S_w     = modelo.ala.S;
b_w     = modelo.ala.b;

Xca_v    = modelo.vertical.Xca;
Zca_v     = modelo.vertical.Zca;

Sref    = modelo.general.Sref;
l_fus = modelo.general.l_fus;
rend    = 0.95;

D_fus   = modelo.fuselaje.D;
if modelo.vertical.Xca < modelo.fuselaje.l
    Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca);
else
    Dfus_v  = 0;
end

S_v     = modelo.vertical.S;
b_v     = modelo.vertical.b;
z_v     = modelo.vertical.Zca;
CLa_v   = modelo.vertical.CLa;
lt_v    = modelo.vertical.Xca;
eta_v   = modelo.vertical.eta;

S_h     = modelo.horizontal.S;
b_h     = modelo.horizontal.b;
CLa_h   = modelo.horizontal.CLa;
eta_h   = modelo.horizontal.eta;
AR_h    = modelo.horizontal.AR;
Dfus_h  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.horizontal.Xca, 'pchip');
z_h     = modelo.horizontal.Zca;
TR_h    = modelo.horizontal.TR;
LAMc2_h = modelo.horizontal.LAMc2;
LAMc4_h = modelo.horizontal.LAMc4;
lt_h    = modelo.horizontal.Xca - modelo.horizontal.xca + modelo.horizontal.le_y(end) + modelo.horizontal.ct/2;
diedro_h  = modelo.horizontal.diedro;

Cyvee_beta     = modelo.vee.CLa;
S_vee          = modelo.vee.S;
dihedral_vee   = modelo.vee.dihedral_w2;  
Cla_vee        = modelo.vee.Cla;
%         modelo.vee.CLa = CLalpha_wb_w2*(cos(dihedral_w2))^2;
%         modelo.vee.b = b_w2_s;
%         modelo.vee.Xca = x_xbar_w2;
%         modelo.vee.Zca = l_zcg_w2;  


x_XCG = modelo.general.Xcg;
z_XCG = modelo.general.Zcg;

%% CÁLCULO DE LA DERIVADA DE ESTABILIDAD
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)
k_prima                 = 1;
% The three dimensional ruddervator effectiveness parameter is determined from Figure 8.53 in Airplane Design Part VI 
% and is a function of the V-tail aspect ratio, and ruddervator chord to V-tail chord ratio:
alphaCL_alphaCl         = alphaCL_alphaCl_calc(cf_c,AR);
% Calculo de Cldelta_Cldeltatheory: Correction Factor for Plain Flap Lift. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.15 (pag 230, 1920 PDF))
Cldelta_Cldeltatheory   = Cldelta_Cldeltatheory_calc(cf_c,Cla);
% Calculo de Cldeltatheory: Lift Effectiveness. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.2.1, Fig. 8.14 (pag 228, 1918 PDF)
Cldelta_theory          = Cldeltatheory_calc(cf_c,t_c);
% Calculo de Kb: Effect of the taper ratio and flap span. AIRPLANE DESIGN, PART VI; Roskam, Jan; 8.1.4.2, Fig. 8.52 (pag 260, 1950 PDF))
K_b                     = Kb_calc(y0_b2, y1_b2, TR);

% The change in sideslip due to ruddervator directional deflection 
beta_delta_rv = -K_b*Cldelta_Cldeltatheory*Cldelta_theory*(k_prima/Cla_vee)*alphaCL_alphaCl;

% The wing-fuselage-"equivalent horizontal tail" interference on the airplane sideforce-coefficient-due-to-sideslip
% derivative of equivalent twin vertical tails is obtained from Figure 10.17 in Airplane Design Part VI and is a function of the equivalent 
% twin vertical tail span, the fuselage depth at the quarter chord point of the V-tail, the distance between the two equivalent vertical tails, 
% and the fuselage length:
Cybetav_Cybetaveff = Cybetav_Cybetaveff_calc(b_h,l_fus,Dfus_v,b_v);
Cy_drv   = -Cybetav_Cybetaveff*rend*Cyvee_beta*(S_vee/S_w)*beta_delta_rv*tan(dihedral_vee);

% The airplane rolling-moment-coefficient-due-to-ruddervator-deflection derivative is determined from:
Cl_drv   = Cybetav_Cybetaveff*rend*Cyvee_beta*(S_vee/S_w)*beta_delta_rv*tan(dihedral_vee)*((Zca_v - z_XCG)*cos(alpha)-(Xca_v - x_XCG)*sin(alpha))/b_w;
Cy_drv = Stab_Der.Cydeltarv
Cl_drv = Cy_drv*((Zca_v - z_XCG)*cos(alpha)-(Xca_v - x_XCG)*sin(alpha))/b_w
if isnan(Cl_drv)
    Cl_drv = 0;
end
Stab_Der.Cldeltarv = Cl_drv;
end