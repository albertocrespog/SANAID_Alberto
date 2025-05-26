%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   27 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Calculo preliminar de la derivada de estabilidad C_l_beta.

function Cl_beta = pre_Clbeta()


%% Cl_beta = Cl_beta_wf + Cl_beta_vert

%%  APORTE DEL ALA-FUSELAJE
% Estimacion de Cl_beta_wf extraida de:
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
%   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 10.368 (pag 301)

DClbeta_tau = -180/pi*0.005*sqrt(AR)*(dfus_w/b_w)^2; % 1/rad
DClbeta_zw = 1.2*sqrt(AR)*z_w/b_w*2*dfus_w/b_w; % 1/rad

Clbeta_CL_LAMc2 = 180/pi*Clbeta_CL_LAM_calc(lambda, AR, LAMc2); % 1/rad
K_MLAM = K_MLAM_calc(M, AR, LAMc2);
Kf = Kf_calc(lt, b, AR, LAMc2);
Clbeta_CL_AR = 180/pi*Clbeta_CL_AR_calc(AR, lambda); % 1/rad
Clbeta_die = 180/pi*Clbeta_die_calc(lambda, AR, LAMc2); % 1/rad
K_Mdie = K_Mdie_calc(M, AR, LAMc2);

Cl_beta_wf = C_L*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR) + die*(Clbeta_die*K_Mdie + DClbeta_tau) + DClbeta_zw;
% die: Angulo de diedro en radianes

%% APORTE DEL VERTICAL
% Calculo del side-wash: deflexion de la corriente debida a la presencia
% del ala
sidewash = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*z_w/D_fus + 0.009*AR_w;
% S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
% LAMc4: flecha del ala en c/4
% z_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

% Calculo de la relacion de aspecto efectiva del vertical
AvB_Av      = AvB_Av_calc(b_v,2*r1,lambda_v);
KH          = KH_calc(S_h, S_v);
AvHB_AvB    = AvHB_AvB_calc(xac_h,b_v,zac_h,c_v);

ARv_eff = AvB_Av*AR_v*(1+KH*(AvHB_AvB-1));

% b_v:      Envergadura del vertical
%
% r1:       Radio medio del fuselaje en donde se encuentra situado el
%           vertical.
%
% lambda_v: Estrechamiento del vertical.
%
% ARv:      Relación de aspecto del vertical.
%
% St:       Superficie del estabilizador horizontal.
%
% Sv:       Superficie del estabilizador vertical.
%
% xac_h:    Distancia relativa entre el centro aerodinamico del estabilizador y el borde de ataque del estabilizador vertical.
%
% zac_h:    Distancia del plano del estabilizador horizontal a la linea
%           central del fuselaje.
%
% AvB_Av:   Ratio of the vertical tail aspect ratio in the presence of the
%           fuselage to that of the isolated vertical tail.
%
% AvHB_AvB: Ratio of the vertical tail aspect ratio in the presence of the
%           fuselage and the horizontal tail to that in the presence of the fuselage alone.
%
% KH:       Factor that accounts for the relative size of the horizontal and the
%           vertical tail.

% Calculo de la pendiente de sustentacion del vertical
beta    = sqrt(1-M^2);
Cla_M   = Cla_0/beta; % Cla_0 Pendiente de la curva de sustentacion del perfil a angulo de ataque nulo
k_cl    = Cla_M/2/pi;
Cla_AR  = Cla_AR_calc(LAMc2_v, k_cl, beta);
CLa_v   = ARv_eff*Cla_AR;     %CLa_v = 2*pi*ARv_eff/(2 + sqrt((ARv_eff*beta/k)^2*(1+(tan(LAMvc2)/beta)^2)+4));
k       = k_calc(b_v, 2*r1);

Cl_beta_vert = -k*CLa_v*sidewash*S_v/Sref*z_v/b_w;

%% DERIVADA DE ESTABILIDAD COMPLETA

Cl_beta = Cl_beta_vert + Cl_beta_wf;

end


