%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   27 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo preliminar de la derivada de estabilidad C_n_beta.

function [Cn_beta_out] = pre_Cnbeta(modelo)

%% Cn_beta = Cn_beta_fus + Cn_beta_wing + Cn_beta_vert

%% APORTE ALA
C_L = modelo.general.W/modelo.general.qinf/modelo.general.Sref;

 Cn_beta_w = C_L^2 * (1/4/pi/modelo.ala.AR - (tan(modelo.ala.LAMc4)/pi/modelo.ala.AR/(modelo.ala.AR + ...
 4*cos(modelo.ala.LAMc4)))*(cos(modelo.ala.LAMc4) - modelo.ala.AR/2 - modelo.ala.AR^2/8/cos(modelo.ala.LAMc4) +...
 (6*(modelo.ala.Xac - modelo.general.Xcg)*sin(modelo.ala.LAMc4))/modelo.ala.MAC/modelo.ala.AR));

Cn_beta_out{2} = Cn_beta_w;

%% APORTE FUSELAJE
% 2 FORMAS:
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 398, 2088 PDF)
%   PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu; 3.5 (pag 271)
% K_N = K_N_calc(l_fus, S_fusS, Xcg, D_fus, h1, h2, W_fus);
% K_Rl = K_Rl_calc(Re);
% Cn_beta_fus = -180/pi*K_N*K_Rl*S_fusS*l_fus/Sref/b_w;

%   - AIRCRAFT DESIGN: A CONCEPTUAL APPROACH 5th EDITION; Raymer, Daniel P.; 16.4.5 (pag 634)
Cn_beta_fus = -1.3*modelo.fuselaje.vol/modelo.general.Sref/modelo.ala.b*modelo.fuselaje.D/modelo.fuselaje.W;
% D_fus: Fuselage depth
% W_fus: Fuselage width
% V_fus: Fuselage volume

Cn_beta_out{3} = Cn_beta_fus;

%% APORTE VERTICAL
% Calculo del side-wash: deflexion de la corriente debida a la presencia
% del ala
sidewash = 0.724 + 3.06*modelo.vertical.S/modelo.general.Sref/(1+cos(modelo.ala.LAMc4)) + 0.4*modelo.ala*Zac/modelo.fuselaje.D + 0.009*modelo.ala.AR;
% S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
% LAMc4: flecha del ala en c/4
% z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))


% Calculo de la pendiente de sustentacion del vertical

Dfus_v = intep1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca,'pchip');

X_le        = modelo.vertical.Xac - modelo.vertical.xac + intep1(modelo.vertical.y, modelo.vertical.le_y, modelo.horizontal.Zac,'pchip');
dx_vle2hac  = horizontal.Xac - X_le;
dz_fus2hac  = horizontal.Zac;

ARv_eff = getARv_eff(modelo.vertical.S, modelo.vertical.b, modelo.vertical.TR,...
modelo.horizontal.S, dx_vle2hac, dz_fus2hac, Dfus_v, modelo.general.vert_conf);

CLa_v = getCLa_geometryBased(ARv_eff, modelo.vertical.LAMc2, modelo.general.Minf, modelo.vertical.Cla);
 

k       = k_calc(modelo.vertical.b, Dfus_v);
l_v     = modelo.vertical.Xca - modelo.general.Xcg;

Cn_beta_vert = k*CLa_v*sidewash*modelo.vertical.S/modelo.general.Sref*l_v/modelo.ala.b;

Cn_beta_out{4} = Cn_beta_vert;

%% DERIVADA DE ESTABILIDAD COMPLETA

Cn_beta = Cn_beta_w + Cn_beta_fus + Cn_beta_vert;

Cn_beta_out{1} =  Cn_beta;

end
