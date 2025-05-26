% Calculo preliminar de la derivada de estabilidad C_l_beta.

function Stab_Der = getClbeta(modelo,alpha,Stab_Der)

%% Cl_beta = Cl_beta_wf + Cl_beta_vert
M       = modelo.general.Minf;
% qinf    = modelo.general.qinf;
% W       = modelo.general.mtow*modelo.general.w_w0*9.8065;
Sref    = modelo.general.Sref;
% C_L      = modelo.general.CL;
C_Lw   = modelo.general.CL_w;
% C_Lh   = modelo.general.CL_h;

D_fus   = modelo.fuselaje.D;
if modelo.vertical.Xca < modelo.fuselaje.l
    Dfus_v  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.vertical.Xca);
else
    Dfus_v  = 0;
end
    
AR_w    = modelo.ala.AR;
b_w     = modelo.ala.b;
Dfus_w  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.ala.Xca, 'pchip');
z_w     = modelo.ala.Zca;
TR_w    = modelo.ala.TR;
LAMc2_w = modelo.ala.LAMc2;
LAMc4_w = modelo.ala.LAMc4;
lt_w    = modelo.ala.Xca - modelo.ala.xca + modelo.ala.le_y(end) + modelo.ala.ct/2;
diedro  = modelo.ala.diedro;

S_v     = modelo.vertical.S;
b_v     = modelo.vertical.b;
Zca_v     = modelo.vertical.Zca;
CLa_v   = modelo.vertical.CLa;
Xca_v    = modelo.vertical.Xca;
eta_v   = modelo.vertical.eta;

% S_h     = modelo.horizontal.S;
% b_h     = modelo.horizontal.b;
% CLa_h   = modelo.horizontal.CLa;
% eta_h   = modelo.horizontal.eta;
% AR_h    = modelo.horizontal.AR;
% Dfus_h  = interp1(modelo.fuselaje.x, modelo.fuselaje.D_x, modelo.horizontal.Xca, 'pchip');
% z_h     = modelo.horizontal.Zca;
% TR_h    = modelo.horizontal.TR;
% LAMc2_h = modelo.horizontal.LAMc2;
% LAMc4_h = modelo.horizontal.LAMc4;
% lt_h    = modelo.horizontal.Xca - modelo.horizontal.xca + modelo.horizontal.le_y(end) + modelo.horizontal.ct/2;
% diedro_h  = modelo.horizontal.diedro;

x_XCG = modelo.general.Xcg;
z_XCG = modelo.general.Zcg;

%%  APORTE DEL ALA-FUSELAJE
% Estimacion de Cl_beta_wf extraida de:
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
%   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)
DClbeta_tau         = -180/pi*0.005*sqrt(AR_w)*(Dfus_w/b_w)^2; % 1/rad
DClbeta_zw          = 1.2*sqrt(AR_w)*z_w/b_w*2*Dfus_w/b_w; % 1/rad
Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(TR_w, AR_w, LAMc2_w); % 1/rad
K_MLAM              = K_MLAM_calc(M, AR_w, LAMc2_w);

Kf                  = Kf_calc(lt_w, b_w, AR_w, LAMc2_w); % Corrected with DATCOM
Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_w, TR_w); % 1/rad
Clbeta_die          = 180/pi*Clbeta_die_calc(TR_w, AR_w, LAMc2_w); % 1/rad
K_Mdie              = K_Mdie_calc(M, AR_w, LAMc2_w); % Corrected with DATCOM

% C_L = W/qinf/Sref;
%% APORTE ALA
Cl_beta_wf = C_Lw*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR) + diedro*(Clbeta_die*K_Mdie + DClbeta_tau) + DClbeta_zw;
% die: Angulo de diedro en radianes

%% APORTE DEL VERTICAL
% Calculo del side-wash: deflexion de la corriente debida a la presencia
% del ala
sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w/D_fus + 0.009*AR_w; % From Álvaro
%% REVISE I think the correct is Dfus_w
sidewash        = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w/Dfus_w + 0.009*AR_w;
% S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
% LAMc4: flecha del ala en c/4
% z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

% The empirical factor for estimating airplane sideforce-coefficient-due-to-sideslip derivative is found from Figure 10.12 in Airplane Design Part VI
% and is a function of the vertical tail span and the height of the fuselage at the quarter chord point of the vertical tail section:
k               = k_calc(b_v, Dfus_v);

Cl_beta_vert    = -k*eta_v*CLa_v*sidewash*S_v/Sref*((Zca_v - z_XCG)*cos(alpha)-(Xca_v - x_XCG)*sin(alpha))/b_w;
if isnan(Cl_beta_vert)
    Cl_beta_vert = 0;
end

% S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
% LAMc4: flecha del ala en c/4
% z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

% % Calculo de la relacion de aspecto efectiva del vertical
% AvB_Av      = AvB_Av_calc(b_v, Wfus_v, TR_v);
% KH          = KH_calc(S_h, S_v);
% AvHB_AvB    = AvHB_AvB_calc(xac_h, b_v, zac_h, c_v);
% 
% ARv_eff = AvB_Av*AR_v*(1+KH*(AvHB_AvB-1));
% 
% % b_v:      Envergadura del vertical
% %
% % r1:       Radio medio del fuselaje en donde se encuentra situado el
% %           vertical.
% %
% % lambda_v: Estrechamiento del vertical.
% %
% % ARv:      Relación de aspecto del vertical.
% %
% % St:       Superficie del estabilizador horizontal.
% %
% % Sv:       Superficie del estabilizador vertical.
% %
% % xac_h:    Distancia relativa entre el centro aerodinamico del estabilizador y el borde de ataque del estabilizador vertical.
% %
% % zac_h:    Distancia del plano del estabilizador horizontal a la linea
% %           central del fuselaje.
% %
% % AvB_Av:   Ratio of the vertical tail aspect ratio in the presence of the
% %           fuselage to that of the isolated vertical tail.
% %
% % AvHB_AvB: Ratio of the vertical tail aspect ratio in the presence of the
% %           fuselage and the horizontal tail to that in the presence of the fuselage alone.
% %
% % KH:       Factor that accounts for the relative size of the horizontal and the
% %           vertical tail.
% 
% % Calculo de la pendiente de sustentacion del vertical
% beta    = sqrt(1-M^2);
% Cla_M   = Cla_0/beta; % Cla_0 Pendiente de la curva de sustentacion del perfil a angulo de ataque nulo
% k_cl    = Cla_M/2/pi;
% Cla_AR  = Cla_AR_calc(LAMc2_v, k_cl, beta);
% CLa_v   = ARv_eff*Cla_AR;     %CLa_v = 2*pi*ARv_eff/(2 + sqrt((ARv_eff*beta/k)^2*(1+(tan(LAMvc2)/beta)^2)+4));


%% DERIVADA DE ESTABILIDAD COMPLETA

Cl_beta         = Cl_beta_vert + Cl_beta_wf;

% Cl_beta_out(1)  = Cl_beta;
% Cl_beta_out(2)  = Cl_beta_vert;
% Cl_beta_out(3)  = Cl_beta_wf;
Stab_Der.Clb_wb = Cl_beta_wf;
Stab_Der.Clb_v = Cl_beta_vert;
Stab_Der.Clb = Cl_beta;


end


