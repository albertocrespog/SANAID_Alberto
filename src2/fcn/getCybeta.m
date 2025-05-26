% Calculo preliminar de la derivada de estabilidad C_l_beta.

% function [Cy_beta_w, Cy_beta_fus, Cy_beta_vert, Cy_beta] = getCybeta(modelo)
function [Stab_Der] = getCybeta(modelo,Stab_Der)
% Datos

qinf    = modelo.general.qinf;
W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
Sref    = modelo.general.Sref;
C_L      = modelo.general.CL;
C_Lw   = modelo.general.CL_w;
% C_Lh   = modelo.general.CL_h;

AR_w    = modelo.ala.AR;

Dfus_w  = interp1(length_x_position,height_x_position,x_xbar_w1, 'pchip');z_w     = modelo.ala.Zca;
LAM_w   = modelo.ala.LAMc2;
LAMc4_w = modelo.ala.LAMc4;
diedro  = modelo.ala.diedro;

S_v     = modelo.vertical.S;
eta_v   = modelo.vertical.eta;
CLa_v   = modelo.vertical.CLa;
b_v     = modelo.vertical.b;
Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca, 'pchip');

vol_fus = modelo.fuselaje.vol;
CLa_fus = modelo.fuselaje.CLa;



%% Cy_beta = Cy_beta_fus + Cy_beta_wing + Cy_beta_vert

%% APORTE ALA
% - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
Cy_beta_w       = -0.00573*abs(diedro*180/pi);

% - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
C_L = W/qinf/Sref;
Cy_beta_w       = C_Lw^2*(6*tan(LAM_w)*sin(LAM_w))/(pi*AR_w*(AR_w + 4*cos(LAM_w)));  

%% APORTE FUSELAJE
Ki  = Ki_calc(z_w, Dfus_w);
% - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
%Cy_beta_fus     = -2*Ki*S0/Sref;

% - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
A_refFus        = vol_fus^(2/3);
Cy_beta_fus     = -Ki*CLa_fus*A_refFus/Sref;

%% APORTE VERTICAL
% AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
% FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

% Calculo del side-wash: deflexion de la corriente debida a la presencia
% del ala
sidewash = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w/Dfus_w + 0.009*AR_w;
% S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
% LAMc4: flecha del ala en c/4
% z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

% The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
% and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
k       = k_calc(b_v, Dfus_v);

% Single Vertical tail
Cy_beta_vert    = -k*CLa_v*sidewash*eta_v*S_v/Sref;
if isnan(Cy_beta_vert)
    Cy_beta_vert = 0;
end

%% DERIVADA TOTAL

Cy_beta         = Cy_beta_w + Cy_beta_fus + Cy_beta_vert;

% Cy_beta_out(1)  = Cy_beta;
% Cy_beta_out(2)  = Cy_beta_vert;
% Cy_beta_out(3)  = Cy_beta_w;
% Cy_beta_out(4)  = Cy_beta_fus;
Stab_Der.Cyb_w = Cy_beta_w;
Stab_Der.Cyb_fus = Cy_beta_fus;
Stab_Der.Cyb_v = Cy_beta_vert;
Stab_Der.Cyb = Cy_beta;

end