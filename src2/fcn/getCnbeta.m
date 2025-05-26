% Calculo preliminar de la derivada de estabilidad C_n_beta.
function Stab_Der = getCnbeta(modelo,Stab_Der,Body_Geo,TRIM_RESULTS,Cy_beta_vert,z_zbar_vert,Geo_tier)

%% Cn_beta = Cn_beta_fus + Cn_beta_wing + Cn_beta_vert

length_x_position = Body_Geo.length_x_position;
height_x_position = Body_Geo.height_x_position;
width_x_position = Body_Geo.width_x_position;
alpha = TRIM_RESULTS.trim_alpha;

Sref    = modelo.general.Sref;
W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
qinf    = modelo.general.qinf;
vinf    = modelo.general.Vinf;
rhoinf  = modelo.general.rhoinf;
Xcg     = modelo.general.Xcg;
Zcg     = modelo.general.Zcg;
C_L     = modelo.general.CL;
C_Lw   = modelo.general.CL_w;
% C_Lh   = modelo.general.CL_h;

AR_w    = modelo.ala.AR;
b_w     = modelo.ala.b;
LAMc4_w = modelo.ala.LAMc4;
% z_w     = modelo.ala.Zca;
% Positrive bellow fuselage center line
z_w1_LE1 = -Geo_tier.z_w1_LE;
x_w     = modelo.ala.Xca;
MAC_w   = modelo.ala.MAC;
Dfus_w  = interp1(length_x_position,height_x_position, x_w, 'pchip');

S_v     = modelo.vertical.S;
b_v     = modelo.vertical.b;
CLa_v   = modelo.vertical.CLa;
x_v     = modelo.vertical.Xca;
Dfus_v = interp1(length_x_position,width_x_position,x_v,'pchip');


l_fus   = modelo.fuselaje.l;
S_fusS  = modelo.fuselaje.Sside;
D_fus   = interp1(length_x_position,height_x_position, x_w, 'pchip');
W_fus   = modelo.fuselaje.W;
h1      = interp1(length_x_position,height_x_position, 0.25*l_fus, 'pchip');
h2      = interp1(length_x_position,height_x_position, 0.75*l_fus, 'pchip');

%% APORTE ALA
Cn_beta_w = C_Lw^2 * (1/4/pi/AR_w - (tan(LAMc4_w)/pi/AR_w/(AR_w + 4*cos(LAMc4_w)))*(cos(LAMc4_w)...
- AR_w/2 - AR_w^2/8/cos(LAMc4_w) + (6*(x_w - Xcg)*sin(LAMc4_w))/MAC_w/AR_w));

%% APORTE FUSELAJE
% 2 FORMAS:
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 398, 2088 PDF)
%   PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu; 3.5 (pag 271)

mu_air      = 1.7e-5;  % Viscosidad dinámica del aire
Re          = rhoinf*vinf*l_fus/mu_air;

K_N     = K_N_calc(l_fus, S_fusS, Xcg, D_fus, h1, h2, W_fus);
K_Rl    = K_Rl_calc(Re);

Cn_beta_fus = -180/pi*K_N*K_Rl*S_fusS*l_fus/Sref/b_w;

%   - AIRCRAFT DESIGN: A CONCEPTUAL APPROACH 5th EDITION; Raymer, Daniel P.; 16.4.5 (pag 634)
% Cn_beta_fus = -1.3*modelo.fuselaje.vol/modelo.general.Sref/modelo.ala.b*modelo.fuselaje.D/modelo.fuselaje.W;
% D_fus: Fuselage depth
% W_fus: Fuselage width
% V_fus: Fuselage volume


%% APORTE VERTICAL
% Calculo del side-wash: deflexion de la corriente debida a la presencia
% del ala
sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w1_LE1/D_fus + 0.009*AR_w;
%% REVISE I think the correct is Dfus_w
% sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w/Dfus_w + 0.009*AR_w;
% S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
% LAMc4: flecha del ala en c/4
% z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))


% Calculo de la pendiente de sustentacion del vertical

% X_le        = modelo.vertical.Xac - modelo.vertical.xac + intep1(modelo.vertical.y, modelo.vertical.le_y, modelo.horizontal.Zac,'pchip');
% dx_vle2hac  = horizontal.Xac - X_le;
% dz_fus2hac  = horizontal.Zac;
% 
% ARv_eff = getARv_eff(S_v, modelo.vertical.b, modelo.vertical.TR,...
% modelo.horizontal.S, dx_vle2hac, dz_fus2hac, Dfus_v, modelo.general.vert_conf);
% 
% CLa_v = getCLa_geometryBased(ARv_eff, modelo.vertical.LAMc2, modelo.general.Minf, modelo.vertical.Cla);
 
k       = k_calc(b_v, Dfus_v);
Cn_beta_vert    = - Cy_beta_vert*((z_zbar_vert - Zcg)*sin(alpha)+(x_v - Xcg)*cos(alpha))/b_w;

if isnan(Cn_beta_vert)
    Cn_beta_vert = 0;
end


%% DERIVADA DE ESTABILIDAD COMPLETA

Cn_beta = Cn_beta_w + Cn_beta_fus + Cn_beta_vert;
% 
% Cn_beta_out(1)  = Cn_beta;
% Cn_beta_out(2)  = Cn_beta_vert;
% Cn_beta_out(3)  = Cn_beta_w;
% Cn_beta_out(4)  = Cn_beta_fus;
% Stab_Der.Cnb_diedro = Cnb_diedro;
% Stab_Der.Cnb_flecha = Cnb_flecha;
Stab_Der.Cnb_w = Cn_beta_w;
Stab_Der.Cnb_b = Cn_beta_fus;
Stab_Der.Cnb_wb = Cn_beta_w+Cn_beta_fus;    
Stab_Der.Cnb_v = Cn_beta_vert;    
Stab_Der.Cnb = Cn_beta;
end
