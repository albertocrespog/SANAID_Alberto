function modelo = calc_derivadasEstabilidad(modelo,alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha       = alpha*pi/180;

g           = 9.8065;
Minf        = modelo.general.Minf;
Vinf        = modelo.general.Vinf;
rhoinf      = modelo.general.rhoinf;
qinf        = modelo.general.qinf;
Sref        = modelo.general.Sref;
Xcg         = modelo.general.Xcg;
MTOW        = modelo.general.mtow*g;
W_W0        = modelo.general.w_w0;
W           = MTOW*W_W0;
CD0         = modelo.general.polar(1);
k1          = modelo.general.polar(2);
k2          = modelo.general.polar(3);

% CL          = W/(0.5*qinf*Sref)
% CD          = CD0 + k1*CL + k2*CL^2
% D           = 0.5*rhoinf*Vinf^2*Sref*CD
% P_req       = Vinf*D

% PROPULSION

Z_Prop      = modelo.propulsion.Z;
X_Prop      = modelo.propulsion.X;
rend_prop   = modelo.propulsion.rendProp;
A_power     = modelo.propulsion.Acoef;
B_power     = modelo.propulsion.Bcoef;
C_power     = modelo.propulsion.Ccoef;
n_eng       = modelo.propulsion.n;
D_prop      = modelo.propulsion.Dprop;
n_blades    = modelo.propulsion.nBlades;
beta_blade  = modelo.propulsion.betaBlade;

if isempty(modelo.propulsion.w_R_30) || isempty(modelo.propulsion.w_R_60) || isempty(modelo.propulsion.w_R_90)
    w_R_30  = 0.0525*2;
    w_R_60  = 0.073*2;
    w_R_90  = 0.045*2;
else
    w_R_30  = modelo.propulsion.w_R_30;
    w_R_60  = modelo.propulsion.w_R_60;
    w_R_90  = modelo.propulsion.w_R_90;
end


% FUSELAJE
Wmax_fus    = modelo.fuselaje.W;
Dmax_fus    = modelo.fuselaje.D;
Sfront_fus  = modelo.fuselaje.Sfront;
Sside_fus   = modelo.fuselaje.Sside;
l_fus       = modelo.fuselaje.l;
Vol_fus     = modelo.fuselaje.vol;
x_fus       = modelo.fuselaje.x;
Wx_fus      = modelo.fuselaje.W_x;
Dx_fus      = modelo.fuselaje.D_x;
Sdistr_fus  = modelo.fuselaje.S_x;
h1_fus      = interp1(x_fus,Dx_fus,0.25*l_fus,'pchip');
h2_fus      = interp1(x_fus,Dx_fus,0.75*l_fus,'pchip');
dSdX_fus    = modelo.fuselaje.dSdX;%diff(Sdistr_fus)./diff(x_fus);
CLa_nose    = modelo.fuselaje.CLa;
alphax_fus  = modelo.fuselaje.alpha_x;

% Check if there is a minimum number of sections (40) to ensure smoothness 
% the smoothness of the different parameters (x_fus, dSdX_fus, ...)
if length(x_fus)<40
    xInt_fus        = linspace(0,l_fus,50);
    dSdX_fus        = interp1(x_fus,dSdX_fus,xInt_fus,'pchp');
    Sdistr_fus      = interp1(x_fus,Sdistr_fus,xInt_fus,'pchp');
    [~,iSmax]       = min(abs(dSdX_fus));
    x_Smax          = xInt_fus(iSmax);
    Sfront_fus      = Sdistr_fus(iSmax);
    Wx_fus          = interp1(x_fus,Wx_fus,xInt_fus,'pchp');
    Dx_fus          = interp1(x_fus,Dx_fus,xInt_fus,'pchp');
    alphax_fus      = interp1(x_fus,alphax_fus,xInt_fus,'pchp');
    x_fus           = xInt_fus;
else
    [~,iSmax]       = min(abs(dSdX_fus));
    x_Smax          = x_fus(iSmax);
end


% ALA
b_w         = modelo.ala.b;
S_w         = modelo.ala.S;
S_we        = modelo.ala.Se;
AR_w        = modelo.ala.AR;
TR_w        = modelo.ala.TR;
cMAC_w      = modelo.ala.MAC;
Xca_w       = modelo.ala.Xca;
Zca_w       = modelo.ala.Zca;
xca_w       = modelo.ala.xca; % Posición del centro aerodinámico respecto al punto mas adelantado del ala
Wala_fus    = interp1(x_fus, Wx_fus, Xca_w, 'pchip'); % Ancho del fuselaje en el encastre del ala
Dala_fus    = interp1(x_fus, Dx_fus, Xca_w, 'pchip'); % Altura del fuselaje en el encastre del ala
XLE_w       = Xca_w - xca_w + interp1(modelo.ala.y, modelo.ala.le_y, Wala_fus/2, 'linear');
b_we        = b_w - Wala_fus;
cr_we       = interp1(modelo.ala.y, modelo.ala.c_y, Wala_fus/2, 'linear');
AR_we       = b_we^2/S_we;
LAM_w       = modelo.ala.LAM;
LAMc4_w     = modelo.ala.LAMc4;
LAMc2_w     = modelo.ala.LAMc2;
TR_we       = modelo.ala.ct/cr_we;
Xt_w        = Xca_w - xca_w + modelo.ala.le_y(end) + modelo.ala.ct/2; % Posición de la punta del ala

CLa_w       = modelo.ala.CLa;
Cla_w       = modelo.ala.Cla;
CLa_we      = CLa_w;
CL0_w       = modelo.ala.CL0;
CM0_w       = modelo.ala.CM0;
eta_w       = modelo.ala.eta;

i_w         = modelo.ala.i;
diedro_w    = modelo.ala.diedro;

alpha0_w    = -CL0_w/CLa_w;

y1_b2_w     = modelo.ala.y1_b2;
y0_b2_w     = modelo.ala.y0_b2;
cm_c_w      = modelo.ala.cm_c;
t_c_w       = modelo.ala.t_c;

% Coeficiente de eficiencia de Oswald del ala (Pamadi, pag. 392)
lambda1_w = AR_w*TR_w/cos(LAM_w);
a1 = 0.0004; a2=-0.008; a3=0.0501; a4=0.8642;
R_w = a1*lambda1_w^3 + a2*lambda1_w^2 + a3*lambda1_w + a4;
eOswald_w = 1.1*CLa_w/(R_w*CLa_w + (1-R_w)*pi*AR_w);

% HORIZONTAL
S_h         = modelo.horizontal.S;
b_h         = modelo.horizontal.b;
AR_h        = modelo.horizontal.AR;
cMAC_h      = modelo.horizontal.MAC;
LAM_h       = modelo.horizontal.LAM;
Xca_h       = modelo.horizontal.Xca;
Zca_h       = modelo.horizontal.Zca;
TR_h        = modelo.horizontal.TR;
eta_h       = modelo.horizontal.eta;
i_h         = modelo.horizontal.i;
y0_b2_h     = modelo.horizontal.y0_b2;
y1_b2_h     = modelo.horizontal.y1_b2;
cm_c_h      = modelo.horizontal.cm_c;
t_c_h       = modelo.horizontal.t_c;

CLa_h       = modelo.horizontal.CLa;
Cla_h       = modelo.horizontal.Cla;
CL0_h       = modelo.horizontal.CL0;
%CM0_h       = modelo.horizontal.CM0;

downw       = downwash_calc(Xca_w, Xca_h, AR_w, TR_w, LAMc4_w, b_w, Zca_h - Zca_w);

% % Coeficiente de eficiencia de Oswald del horizontal (Pamadi, pag. 392)
% lambda1_h = AR_h*TR_h/cos(LAM_h);
% a1 = 0.0004; a2=-0.008; a3=0.0501; a4=0.8642;
% R_h = a1*lambda1_h^3 + a2*lambda1_h^2 + a3*lambda1_h + a4;
% eOswald_h = 1.1*CLa_h/(R_h*CLa_h + (1-R_h)*pi*AR_h)

% CANARD
S_c         = modelo.canard.S;
cMAC_c      = modelo.canard.MAC;
Xca_c       = modelo.canard.Xca;
TR_c        = modelo.canard.TR;
AR_c        = modelo.canard.AR;
LAM_c       = modelo.canard.LAM;
eta_c       = modelo.canard.eta;
i_c         = modelo.canard.i;
y0_b2_c     = modelo.canard.y0_b2;
y1_b2_c     = modelo.canard.y1_b2;
cm_c_c      = modelo.canard.cm_c;
t_c_c       = modelo.canard.t_c;

CLa_c       = modelo.canard.CLa;
Cla_c       = modelo.canard.Cla;
CL0_c       = modelo.canard.CL0;
%CM0_c       = modelo.canard.CM0;

upw         = upwash_calc(AR_w, Xca_w, Xca_c, cMAC_w);

% lambda1_c = AR_c*TR_c/cos(LAM_c);
% a1 = 0.0004; a2=-0.008; a3=0.0501; a4=0.8642;
% R_c = a1*lambda1_c^3 + a2*lambda1_c^2 + a3*lambda1_c + a4;
% eOswald_c = 1.1*CLa_c/(R_c*CLa_c + (1-R_c)*pi*AR_c)


% VERTICAL

S_v         = modelo.vertical.S;
b_v         = modelo.vertical.b;
cMAC_v      = S_v/b_v;
AR_v        = modelo.vertical.AR;
TR_v        = modelo.vertical.TR;
Xca_v       = modelo.vertical.Xca;
Zca_v       = modelo.vertical.Zca;
if Xca_v > l_fus
    Dvert_fus   = Dx_fus(end);
    Wvert_fus   = Wx_fus(end);
else
    Dvert_fus   = interp1(x_fus, Dx_fus, Xca_v, 'pchip'); % Altura del fuselaje donde se encuentra el vertical
    Wvert_fus   = interp1(x_fus, Wx_fus, Xca_v, 'pchip'); % Anchura del fuselaje donde se encuentra el vertical
end
CLa_v       = modelo.vertical.CLa;
Cla_v       = modelo.vertical.Cla;
eta_v       = modelo.vertical.eta;
LAMc2_v     = modelo.vertical.LAM;

cm_c_v      = modelo.vertical.cm_c;
t_c_v       = modelo.vertical.t_c;
y0_b2_v     = modelo.vertical.y0_b2;
y1_b2_v     = modelo.vertical.y1_b2;

AvB_Av      = AvB_Av_calc(b_v,Dvert_fus,TR_v);
AvHB_AvB    = AvHB_AvB_calc((Xca_h-Xca_v)+0.25*cMAC_v, b_v, Zca_h, cMAC_v);
KH          = KH_calc(S_h, S_v);
AReff_v     = AvB_Av*AR_v*(1 + KH*(AvHB_AvB - 1));
CLa_v       = getCLa_geometryBased(AReff_v,LAMc2_v,Minf,Cla_v);


%% PRELIMINARY CALCULATIONS

if isempty(k2) || isnan(k2) || k2==0
    %k2 = (S_c/Sref)*(1/(pi*eOswald_c*AR_c)) + (S_w/Sref)*(1/(pi*eOswald_w*AR_w)) + (S_h/Sref)*(1/(pi*eOswald_h*AR_h));
    k2 = (S_w/Sref)*(1/(pi*eOswald_w*AR_w));
end

if isempty(k1)
    k1 = 0;
end

Pmax        = (A_power*Vinf^2 + B_power*Vinf + C_power)*rend_prop;

CL          = W/(qinf*Sref);
CD          = CD0 + k1*CL + k2*CL^2;
D           = qinf*Sref*CD;
P_req       = Vinf*D;
delta_T     = P_req/Pmax;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% ESTABILIDAD LONGITUDINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUSELAGE LIFT SLOPE

% f_k2_k1     = k2_k1_calc(l_fus, Dmax_fus);      % Fuselage apparent mass coefficient. Pamadi, Figure 3.6
% CLa_nose  = 2*f_k2_k1*(Sfront_fus/S_w);         % Pendiente de sustentación del morro aislado

KN  = (CLa_nose/CLa_we)*S_w/S_we;
KWB = 0.1714*(Wmax_fus/b_w)^2 + 0.8326*(Wmax_fus/b_w) + 0.9974;
KBW = 0.7810*(Wmax_fus/b_w)^2 + 1.1976*(Wmax_fus/b_w) + 0.0088;

CLa_wb = CLa_we*KWB*(S_we/S_w); % Ecuacion 3.33
CLa_bw = CLa_we*KBW*(S_we/S_w); % Ecuacion 3.34
CLa_WB = (KN + KWB + KBW)*CLa_we*S_we/S_w; % Pendiente de sustentación del conjunto ala-fuselaje Ecuacion 3.24


%% AERODINAMIC CENTER

% Fuselaje
% i=1;
% while dSdX_fus(i) > 0
%     i=i+1;
% end
% x_Smax  = x_fus(i);
l_N     = XLE_w; 
int_Sb      = @(x) interp1(x_fus, dSdX_fus,x).*(l_N - x);
Xca_cre_N   = -1/(cr_we*Sfront_fus)*quad(@(x) int_Sb(x), 0, x_Smax);
Xca_cre_WB  = (Xca_w-l_N)/cr_we;

beta    = sqrt(1-Minf^2);
chi1    = chi_calc(b_w, Wmax_fus);   

if beta*AR_we >= 4
    Xca_cre_BW = 1/4 + (b_w - Wmax_fus)/(2*cr_we)*chi1*tan(LAMc4_w); % Eq 3.39 Pamadi
else
    % Metodo extraido de Pamadi, pagina 189, último parrafo
    % Interpolación lineal con los valores de ARexp_w*beta = 0,4
    Xca_cre_BW4     = 1/4 + (b_w - Wmax_fus)/(2*cr_we)*chi*tan(LAMc4_w);
    Xca_cre_BW0     = xac_cre_BW0_calc(AR_we, TR_we, LAM_w);
    Xca_cre_BW      = interp1([0 4], [Xca_cre_BW0 Xca_cre_BW4], beta*AR_we, 'linear');
end

% Centro aerodinámico del conjunto ala-fuselaje
Xca_WB_cre  = (Xca_cre_N*CLa_nose + Xca_cre_WB*CLa_wb + Xca_cre_BW*CLa_bw)/CLa_WB; % Eq 3.32 Pamadi
Xca_WB      = Xca_WB_cre*cr_we + XLE_w;

% Pendiente de sustentación de la aeronave completa
if isempty(S_h) || strcmp(modelo.conf,'canard')
    CLatotal_h = 0;
else
    CLatotal_h = CLa_h*(1-downw)*eta_h*S_h/Sref;
end

if isempty(S_c) || strcmp(modelo.conf,'convencional')
    CLatotal_c = 0;
else
    CLatotal_c = CLa_c*(1+upw)*eta_c*S_c/Sref;
end

CLa_total = CLa_WB*S_w/Sref + CLatotal_h + CLatotal_c;


% Cálculo del punto neutro de la aeronave completa con todas las
% superficies estabilizadoras horizontales
XNP_WB = CLa_WB*S_w/Sref*Xca_WB;

if isempty(S_c) || strcmp(modelo.conf,'convencional')
    XNP_c = 0;
else
    XNP_c = CLatotal_c*Xca_c;
end

if isempty(S_h) || strcmp(modelo.conf,'canard')
    XNP_h = 0;
else
    XNP_h = CLatotal_h*Xca_h;
end

X_NP    = (XNP_WB + XNP_h + XNP_c)/CLa_total;
SM      = (X_NP - Xcg)/cMAC_w;

% Cálculo del CL0 del avión completo
epsilon0_h  = 2*CL0_w/(pi*AR_w);
epsilon0_c  = epsilon0_h;

if isempty(S_h) || strcmp(modelo.conf,'canard')
    CL0total_h = 0;
else
    CL0total_h = CL0_h*eta_h*S_h/Sref + CLa_h*(i_h - epsilon0_h)*eta_h*S_h/Sref;
end

if isempty(S_c) || strcmp(modelo.conf,'convencional')
    CL0total_c = 0;
else
    CL0total_c = CL0_c*eta_c*S_c/Sref + CLa_c*(i_c + epsilon0_c)*eta_c*S_c/Sref;
end


CL0_total = CL0_w*eta_w*S_w/Sref + CLa_we*eta_w*i_w*S_we/Sref + ...
            CL0total_h + CL0total_c;
        
% Cálculo de las derivadas de control longitudinal
% Horizontal
if isempty(S_h) || strcmp(modelo.conf,'canard')
    modelo.derivadas.CL_de  = 0;
    modelo.derivadas.CD_de  = 0;    
    modelo.derivadas.CM_de  = 0;
else
    Kb_h                = Kb_calc(y0_b2_h, y1_b2_h, TR_h);              % Fig 3.36 Pamadi
    Clde_Cldetheo_h     = Cldelta_Cldeltatheory_calc(cm_c_h, Cla_h);    % Fig 3.37b Pamadi
    Cldetheo_h          = Cldeltatheory_calc(cm_c_h, t_c_h);            % Fig 3.37a Pamadi
    alphaCL_alphaCl_h   = alphaCL_alphaCl_calc(cm_c_h ,AR_h);           % Fig 3.35 Pamadi
    alpha_de_h          = Kb_h*Clde_Cldetheo_h*Cldetheo_h/Cla_h*alphaCL_alphaCl_h;
    CLde_h              = alpha_de_h*CLa_h*eta_h*S_h/Sref;
    CDde_h              = (2*k2*CL + k1)*CLde_h;
    CMde_h              = - CLde_h*(Xca_h - Xcg)/cMAC_w;

    modelo.derivadas.CL_de  = CLde_h;
    modelo.derivadas.CD_de  = CDde_h;  
    modelo.derivadas.CM_de  = CMde_h;
end

% Canard
if isempty(S_c) || strcmp(modelo.conf,'convencional')
    modelo.derivadas.CL_dc      = 0;
    modelo.derivadas.CD_dc      = 0;
    modelo.derivadas.CM_dc      = 0;
else
    Kb_c                = Kb_calc(y0_b2_c, y1_b2_c, TR_c);              % Fig 3.36 Pamadi
    Clde_Cldetheo_c     = Cldelta_Cldeltatheory_calc(cm_c_c, Cla_c);    % Fig 3.37b Pamadi
    Cldetheo_c          = Cldeltatheory_calc(cm_c_c, t_c_c);            % Fig 3.37a Pamadi
    alphaCL_alphaCl_c   = alphaCL_alphaCl_calc(cm_c_c ,AR_c);
    alpha_de_c          = Kb_c*Clde_Cldetheo_c*Cldetheo_c/Cla_c*alphaCL_alphaCl_c;
    CLde_c              = alpha_de_c*CLa_c*eta_c*S_c/Sref;
    CDde_c              = (2*k2*CL + k1)*CLde_c;
    CMde_c              = - CLde_c*(Xca_c - Xcg)/cMAC_w;
    
    modelo.derivadas.CL_dc  = CLde_c;
    modelo.derivadas.CD_dc  = CDde_c;
    modelo.derivadas.CM_dc  = CMde_c;
end


% Calculo de CMa_fus
% Pamadi, página 170, ecuacion 3.7
f_k2_k1         = k2_k1_calc(l_fus, Dmax_fus);
intUp_lim       = l_fus;
intLow_lim      = 0;
Winterp_fus     = @(x) interp1(x_fus, Wx_fus, x, 'pchip');
deps_dalpha     = @(x) depsilon_dalpha_calc(x, cr_we, XLE_w, l_fus, downw,CLa_WB);
int_CMa2         = @(x) (Winterp_fus(x)).^2.*deps_dalpha(x);
int_CMa        = @(x) (Winterp_fus(x)).^2;
CMa_B2         = f_k2_k1*pi/(2*S_w*cMAC_w)*quad(@(x) int_CMa(x), intLow_lim, intUp_lim); % Peor resultado que el empleado finalmente
CMa_B3        = pi/(2*S_w*cMAC_w)*quad(@(x) int_CMa2(x), intLow_lim, intUp_lim);

% Calculo de CMa_B pagina 399 Pamadi, DATCOM 4.2.2.1 3 (PDF 887)
VB_1        = Vol_fus/(Sfront_fus*l_fus);       % Eq 4.513 Pamadi
xm_1        = Xcg/l_fus;                        % Eq 4.513 Pamadi

int_Sb1     = @(x) interp1(x_fus, Sdistr_fus, x).*(x);
Xc_fus      = (1/Vol_fus)*quad(@(x) int_Sb1(x),0,l_fus);
xc_1        = Xc_fus/l_fus;                % Eq 4.513 Pamadi


int_CmaB    = @(x) interp1(x_fus, dSdX_fus, x).*(Xcg - x);     
CMa_B       = 2*f_k2_k1/Vol_fus*quad(@(x) int_CmaB(x),0,x_Smax);

% Calculo de CM0_fus
% Pamadi, pagina 171, ecuacion 3.8

alphax_fus_interp       = @(x) interp1(x_fus, alphax_fus, x, 'spline');
int_CM0                 = @(x) (Winterp_fus(x)).^2.*(alpha0_w + alphax_fus_interp(x));
CM0_fus                 = f_k2_k1/(cMAC_w)*quad(@(x) int_CM0(x), intLow_lim, intUp_lim);


% Cálculo de CM_T1. Momento generado por los motores
% T_prop = rend_prop*Pmax*delta_T/Vinf;
T_prop  = D;
CMT1 = -T_prop*Z_Prop/(qinf*Sref*cMAC_w);


% Cálculo de CM0 del avión completo
if isempty(S_h) || strcmp(modelo.conf,'canard')
    CM0_h = 0;
else
    CM0_h = ((Xcg - Xca_h)/cMAC_w)*(eta_h*S_h/Sref)*(CL0_h + CLa_h*(i_h - epsilon0_h));
end

if isempty(S_c) || strcmp(modelo.conf,'convencional')
    CM0_c = 0;
else
    CM0_c = ((Xcg - Xca_c)/cMAC_w)*(eta_c*S_c/Sref)*(CL0_c + CLa_c*(i_c + epsilon0_c));
end

CM0_total = CM0_fus*S_w/Sref + CM0_w*S_w/Sref + CM0_h + ...
CM0_c + ((Xcg - Xca_w)/cMAC_w)*(CL0_w*S_w/Sref + CLa_we*i_w*S_we/Sref) + ...
CM0_h + CM0_c + CMT1;


% Coeficiente de momentos de la traccion de los motores
CMa_T = 0;

% Coefiente de momentos de la aeronave completa

%CMa_total = -CLa_total*SM + CMa_T + CMa_B*S_w/Sref;
CMa_total = -CLa_total*SM + CMa_T;

modelo.derivadas.CM0    = CM0_total;
modelo.derivadas.CL0    = CL0_total;

modelo.derivadas.CL     = CL;
modelo.derivadas.CD     = CD;

%% ALPHA STABILITY DERIVATIVES
% CD = CD0 + k1*CL + k2*CL^2

CLa         = CLa_total;
CDa         = (k1 + 2*k2*CL)*CLa;
CMa         = CMa_total;


modelo.derivadas.CL_a    = CLa;
modelo.derivadas.CD_a    = CDa;
modelo.derivadas.CM_a    = CMa;


%% SPEED STABILITY DERIVATIVES

if Minf < 0.3
    % Para valores de Mach por debajo de 0.3 ni el coeficiente de 
    % sustentación ni el de resistencia experimenta variaciones con la velocidad
    CDu = 0;
    CLu = 0;
    CMu = 0;
else
    CLu = Minf^2/(1-Minf^2)*CL; % Airplane Flight Dynamics and Automatic Flight Controls; Roskam, Jan; Ecuación 3.119
    CDu = (k1 + 2*k2*CL)*CLu;   % Airplane Flight Dynamics and Automatic Flight Controls; Roskam, Jan; Ecuación 3.107
    CMu = 0;                    % No se considera la variación del centro aerodinámico de las superficies sustentadoras con la velocidad
end

modelo.derivadas.CL_u    = CLu;
modelo.derivadas.CD_u    = CDu;
modelo.derivadas.CM_u    = CMu;


%% PITCH RATE STABILITY DERIVATIVES
% CDq
CDq     = 0; % Airplane Design; Roskam, Jan; Part VI, pag 424 (2114 PDF)

% CLq
CLq_w   = (0.5 + 2*(Xca_w - Xcg)/cMAC_w)*CLa_we;
CLq_B   = 2*CLa_nose*Vol_fus^(2/3)/Sfront_fus*(1-Xcg/l_fus);    
CLq_WB  = (KWB + KBW)*S_we/Sref*CLq_w + CLq_B*Sfront_fus*l_fus/(Sref*cMAC_w); % Eq 4.488 pamadi

if isempty(S_h) || strcmp(modelo.conf,'canard')
    CLq_h   = 0;
else
    CLq_h   = 2*CLa_h*eta_h*(Xca_h - Xcg)/cMAC_w*S_h/Sref;      % Airplane Design; Roskam, Jan; Part VI, pag 426 (2115 PDF)
end
if isempty(S_c) || strcmp(modelo.conf,'convencional')
    CLq_c   = 0; 
else
    CLq_c   = -2*CLa_c*eta_c*(Xcg - Xca_c)/cMAC_w*S_c/Sref;     % Airplane Design; Roskam, Jan; Part VI, pag 426 (2115 PDF)
end

CLq     = CLq_WB + CLq_h + CLq_c;

% CMq
B_prandtl   = sqrt(1 - (Minf*cos(LAMc4_w))^2);
% Se definen una serie de coeficientes dados en la pagina 398 del Pamadi
c1          = AR_we^3*(tan(LAMc4_w))^2; % Eq 4.505 Pamadi
c2          = 3/B_prandtl; % Eq 4.506 Pamadi
c3          = AR_we*B_prandtl + 6*cos(LAMc4_w); % Eq 4.507 Pamadi
c4          = AR_we + 6*cos(LAMc4_w); % Eq 4.508 Pamadi
c5          = AR_we + 2*cos(LAMc4_w); % Eq 4.509 Pamadi

xi         = (Xca_w - Xcg)/cMAC_w;
CMq_e_M02   = -0.7*Cla_w*cos(LAMc4_w)*((AR_we*(0.5*xi + 2*xi^2))/c5 + (c1/(24*c4)) + 1/8); % Eq 4.504 Pamadi
CMq_e       = CMq_e_M02*(c1/c3 + c2)/(c1/c4 + 3);   % Eq 4.503 Pamadi %% REVISADA CON DATCOM, sección 7.1.1.2 (2489 pdf)
    


CMq_B       = 2*CMa_B*VB_1*((1 - xm_1)^2 - Vol_fus*(xc_1-xm_1))/(1-xm_1-VB_1);          % Eq 4.512 Pamadi
CMq_B       = 0; %% Se elimina el efecto del fuselaje
CMq_WB      = (KWB + KBW)*(S_we/Sref)*CMq_e + CMq_B*Sfront_fus/Sref*(l_fus/cMAC_w)^2;  % Eq 4.502 Pamadi

if isempty(S_h) || strcmp(modelo.conf,'canard')
    CMq_h       = 0;
else
    CMq_h       = -2*CLa_h*eta_h*S_h/Sref*(Xca_h - Xcg)^2/cMAC_w^2;  % Airplane Design; Roskam, Jan; Part VI, pag 426 (2116 PDF)
end

if isempty(S_c) || strcmp(modelo.conf,'convencional')
    CMq_c       = 0;
else
    CMq_c       = -2*CLa_c*eta_c*S_c/Sref*(Xca_c - Xcg)^2/cMAC_w^2;  % Airplane Design; Roskam, Jan; Part VI, pag 426 (2116 PDF)
end

% LOS APORTES DE CANARD Y HORIZONTAL SALEN CON EL MISMO SIGNO, HABRIA QUE CONFIRMAR QUE ESTO SEA CORRECTO

CMq         = CMq_WB + CMq_c + CMq_h;

modelo.derivadas.CL_q   = CLq;
modelo.derivadas.CD_q   = CDq;
modelo.derivadas.CM_q   = CMq;

%% ACCELERATION DERIVATIVE ALPHA DOT
% El método del DATCOM es aplicable para alas triangulares, las cuales
% llevan asociadas relaciones de aspecto pequeñas. Debe ser por esto que se
% dispara el valor del coeficiente C_L(g), el cual depende de beta*AR

% Dado este inconveniente se recurre a los métodos empleados en Airplane Design
% CDaDot
CDaDot     = 0; % Airplane Design; Roskam, Jan; Part VI, pag 381 (2071 PDF)

% CLaDot
CLaDot_fus  = 2*CLa_nose*Vol_fus^(2/3+1)/Sfront_fus^2/l_fus;    % Pamadi, 4.530
CLaDot_fus  = 0; %% Se elimna el efecto del fuselaje

CLaDot_e = 1.5*xca_w/cr_we*CLa_we;

CLaDot_WB   = (KWB + KBW)*(S_we/Sref)*(cMAC_w/cr_we)*CLaDot_e + CLaDot_fus*Sfront_fus*l_fus/Sref/cMAC_w;

if isempty(S_h) || strcmp(modelo.conf,'canard')
    CLaDot_h = 0;
else
    CLaDot_h    = 2*eta_h*CLa_h*(Xca_h - Xcg)/cMAC_w*S_h/Sref*(downw); % Airplane Design; Roskam, Jan; Part VI, pag 381 (2071 PDF)
end
                                                                
CLaDot      = CLaDot_WB + CLaDot_h;                             % Pamadi, 4.527

% CMaDot    % REVISAR ESTA DERIVADA, FALTA UN TERMINO
CMaDot_e        = -81/32*(xca_w/cr_we)^2*CLa_we + (XLE_w- Xcg)/cr_we*CLaDot_e;
CMaDot_B1        = 2*CMa_B*VB_1^2*(xc_1-xm_1)/(1-xm_1-VB_1);
CMaDot_B        = 0; %% Se elimina el efecto del fuselaje
CMaDot_WB       = (KWB + KBW)*(S_we/Sref)*(cMAC_w/cr_we)^2*CMaDot_e + CMaDot_B*Sfront_fus/Sref*(l_fus/cMAC_w)^2;

CMaDot_h     = -2*eta_h*CLa_h*(Xca_h - Xcg)^2/cMAC_w^2*S_h/Sref*(downw); % Airplane Design; Roskam, Jan; Part VI, pag 381 (2071 PDF)

CMaDot      = CMaDot_WB + CMaDot_h;

modelo.derivadas.CL_alphaDot    = CLaDot;
modelo.derivadas.CD_alphaDot    = CDaDot;
modelo.derivadas.CM_alphaDot    = CMaDot;

%% PROPULSIVE DERIVATIVES

% Coeffcicients of second order Power model
% Coeficientes del Roskam, ver diapositivas longitudinal - Tema 14.2
% P = A_power*V^2 + B_power*V + C_power



%% THRUST VS AOA DERIVATIVE
% Airplane Design; Roskam, Jan; Part VI, pag 381 (2071 PDF)
% Airplane Design; Roskam, Jan; Part VI, pag 340 (2030 PDF)

rho_met2imp = 0.00194032033; % slugs/ft3
m22ft2      = 10.7639104;
m2ft        = 3.2808399;
N2lbf       = 0.224808943;
W2hp        = 0.00134102209;


Pprop_SET   = (C_power + B_power*Vinf + A_power*Vinf^2)/n_eng;
Pprop_SHP   = Pprop_SET*W2hp;

% for jj = 1:n_eng
%     K_T(jj)     = 550*Pprop_SHP(jj)*sqrt(rhoinf*rho_met2imp)/(sqrt((2*W*N2lbf/(Sref*m22ft2))^3)*(D_prop(jj)*m2ft)^2); 
%     dTc_dCL(jj) = (3/2)*rend_prop(jj)*K_T(jj)*sqrt(CL);
% end

for jj = 1:n_eng
    K_T(jj)     = 550*Pprop_SHP*sqrt(rhoinf*rho_met2imp)/(sqrt((2*W*N2lbf/(Sref*m22ft2))^3)*(D_prop*m2ft)^2); 
    dTc_dCL(jj) = (3/2)*rend_prop*K_T(jj)*sqrt(CL);
end

dCM_dCL_TL = 2/((Sref*m22ft2*cMAC_w*m2ft))*sum(((D_prop*m2ft)^2)*(-Z_Prop)*m2ft*dTc_dCL);

% for jj = 1:n_eng
%     depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop(jj), cr_we, XLE_w, l_fus, downw); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
%     CNa_807         = CNa_807_calc(n_blades(jj), beta_blade(jj));
%     KN              = 262*(2*w_R_30(jj)/D_prop(jj)) + 262*(2*w_R_60(jj)/D_prop(jj)) + 135*(2*w_R_90(jj)/D_prop(jj));
%     dCN_dalpha(jj)  = CNa_807*(1 + 0.8*(KN/80.7 - 1)); 
% end
% dCM_dCL_N   = (pi/4)*(Sref*m22toft2*cMAC_w*m2ft*CLa_w)*sum(dCN_dalpha.*(1 + depsu_dalpha).*(Xcg - X_Prop)*m2ft.*(D_prop*m2ft)^2);

for jj = 1:n_eng
    depsu_dalpha(jj)= depsilon_dalpha_calc(X_Prop, cr_we, XLE_w, l_fus, downw, CLa_WB*S_w/Sref); %% HAY QUE ESTUDIAR EL SIGNO DEL DOWNWASH
    CNa_807         = CNa_807_calc(n_blades, beta_blade);
    KN              = 262*(2*w_R_30/D_prop) + 262*(2*w_R_60/D_prop) + 135*(2*w_R_90/D_prop);
    dCN_dalpha(jj)  = CNa_807*(1 + 0.8*(KN/80.7 - 1)); 
end

dCM_dCL_N   = (pi/4)*(Sref*m22ft2*cMAC_w*m2ft*CLa_w)*sum((dCN_dalpha.*(1 + depsu_dalpha))*(Xcg - X_Prop)*m2ft*(D_prop*m2ft)^2);

Delta_CM_CL_T = dCM_dCL_TL + dCM_dCL_N;

CMTa = Delta_CM_CL_T*CLa;

modelo.derivadas.CM_Ta  = CMTa;

%% THRUST VS SPEED DERIVATIVES

% Airplane Design; Roskam, Jan; Part VI, pag 377 (2067 PDF)
% Diapositivas TEMA 14_2 Derivadas Estabilidad Longitudinal

pres_met2imp    = N2lbf/m22ft2;
CTx1    = CD;
if modelo.propulsion.F_OEI == 1.25
    CTxu    = (2*A_power*Vinf + B_power)*W2hp/m2ft/(qinf*pres_met2imp*Sref*m22ft2) - 2*CTx1;
else
    CTxu    = - 3*CTx1;
end
CMTu    = -(Z_Prop/cMAC_w)*CTxu;
CMT1    = -(Z_Prop/cMAC_w)*CTx1;

CTxa    = 0;

modelo.derivadas.CT_x1  = CTx1;
modelo.derivadas.CT_xu  = CTxu;
modelo.derivadas.CT_xa  = CTxa;
modelo.derivadas.CM_Tu  = CMTu;
modelo.derivadas.CM_T1  = CMT1;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% ESTABILIDAD LATERAL-DIRECCIONAL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIDESLIP ANGLE STABILITY DERIVATIVES
% Cy_beta = Cy_beta_fus + Cy_beta_wing + Cy_beta_vert

% Ala
% - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
Cy_beta_w       = -0.00573*abs(diedro_w*180/pi);

% - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
Cy_beta_w       = CL^2*(6*tan(LAM_w)*sin(LAM_w))/(pi*AR_we*(AR_we + 4*cos(LAM_w)));  


% Fuselaje
Ki  = Ki_calc(Zca_w, Dala_fus);
% - METODO 1:  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 383, 2073 PDF)
%Cy_beta_fus     = -2*Ki*S0/Sref;

% - METODO 2:  FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 321)
A_refFus        = Vol_fus^(2/3);
Cy_beta_fus     = -Ki*CLa_nose*A_refFus/Sref;

% Vertical
% AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 386, 2076 PDF)
% FLIGHT VEHICLE PERFORMANCE AND AERODYNAMIC CONTROL; Smetana, Frederik (pag 322)

% Calculo del side-wash: deflexion de la corriente debida a la presencia
% del ala
sidewash = 0.724 + 3.06*S_v/Sref/(1+cos(LAMc4_w)) + 0.4*Zca_w/Dala_fus + 0.009*AR_we;
% S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
% LAMc4: flecha del ala en c/4
% Zca_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

k               = k_calc(b_v, Dvert_fus);

Cy_beta_vert    = -k*CLa_v*sidewash*S_v/Sref;

if isnan(Cy_beta_vert)
    Cy_beta_vert = 0;
end

% DERIVADA TOTAL
Cy_beta         = Cy_beta_w + Cy_beta_fus + Cy_beta_vert;

modelo.derivadas.Cy_beta    = Cy_beta;

%% CyT_beta Airplane Design (2088 PDF)
CyT_beta    = (pi/4)/(Sref*m22ft2)*sum(((D_prop*m2ft).^2).*dCN_dalpha);

modelo.derivadas.Cy_Tbeta   = CyT_beta;

%% Cl_beta = Cl_beta_wf + Cl_beta_vert

% Ala-fuselaje
% Estimacion de Cl_beta_wf extraida de:
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Eq. 10.34 (pag 392, 2082 PDF)
%   - PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu N.; 3.3 Eq. 3.368 (pag 301)

DClbeta_tau         = -180/pi*0.005*sqrt(AR_we)*(Dala_fus/b_w)^2; % 1/rad
DClbeta_zw          = 1.2*sqrt(AR_we)*Zca_w/b_w*2*Dala_fus/b_w; % 1/rad

Clbeta_CL_LAMc2     = 180/pi*Clbeta_CL_LAM_calc(TR_we, AR_we, LAMc2_w); % 1/rad
K_MLAM              = K_MLAM_calc(Minf, AR_we, LAMc2_w);
Kf                  = Kf_calc(Xt_w, b_w, AR_we, LAMc2_w);
Clbeta_CL_AR        = 180/pi*Clbeta_CL_AR_calc(AR_we, TR_we); % 1/rad
Clbeta_die          = 180/pi*Clbeta_die_calc(TR_we, AR_we, LAMc2_w); % 1/rad
K_Mdie              = K_Mdie_calc(Minf, AR_we, LAMc2_w);

%C_L = W/qinf/Sref;

Cl_beta_wf = CL*(Clbeta_CL_LAMc2*K_MLAM*Kf + Clbeta_CL_AR) + diedro_w*(Clbeta_die*K_Mdie + DClbeta_tau) + DClbeta_zw;
% diedro: Angulo de diedro en radianes

% Vertical

Cl_beta_vert    = -k*CLa_v*sidewash*S_v/Sref*Zca_v/b_w;

if isnan(Cl_beta_vert)
    Cl_beta_vert = 0;
end

% DERIVADA TOTAL
Cl_beta         = Cl_beta_vert + Cl_beta_wf;

modelo.derivadas.Cl_beta    = Cl_beta;

%% Cn_beta = Cn_beta_fus + Cn_beta_wing + Cn_beta_vert

% Ala
Cn_beta_w = CL^2*(1/4/pi/AR_we - (tan(LAMc4_w)/pi/AR_w/(AR_we + 4*cos(LAMc4_w)))*(cos(LAMc4_w)...
- AR_we/2 - AR_we^2/8/cos(LAMc4_w) + (6*(Xca_w - Xcg)*sin(LAMc4_w))/cMAC_w/AR_we));


% Fuselaje
% 2 FORMAS:
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 398, 2088 PDF)
%   PERFORMANCE, STABILITY, DYNAMICS AND CONTROL; Pamadi, Bandu; 3.5 (pag 271)

mu_air      = 1.7e-5;  % Viscosidad dinámica del aire
Re          = rhoinf*Vinf*l_fus/mu_air;

K_N     = K_N_calc(l_fus, Sside_fus, Xcg, Dmax_fus, h1_fus, h2_fus, Wmax_fus);
K_Rl    = K_Rl_calc(Re);

Cn_beta_fus = -180/pi*K_N*K_Rl*Sside_fus*l_fus/Sref/b_w;

%   - AIRCRAFT DESIGN: A CONCEPTUAL APPROACH 5th EDITION; Raymer, Daniel P.; 16.4.5 (pag 634)
% Cn_beta_fus = -1.3*modelo.fuselaje.vol/modelo.general.Sref/modelo.ala.b*modelo.fuselaje.D/modelo.fuselaje.W;
% D_fus: Fuselage depth
% W_fus: Fuselage width
% V_fus: Fuselage volume


% Vertical

% S_v: effective vertical surface (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 387, 2077 PDF))
% LAMc4: flecha del ala en c/4
% Zca_w: distancia del ala a linea central de fuselaje, >0 cuando es ala baja (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1, Fig. 10.9 (pag 384, 2074 PDF))

%Cn_beta_vert2 = k*CLa_v*sidewash*S_v/Sref*(Xca_v - Xcg)/b_w;     
Cn_beta_vert  = - Cy_beta_vert*(Xca_v - Xcg)/b_w;         %% (AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.2.4.1 (pag 398, 2088 PDF))

if isnan(Cn_beta_vert)
    Cn_beta_vert = 0;
end

% DERIVADA DE ESTABILIDAD COMPLETA
Cn_beta     = Cn_beta_w + Cn_beta_fus + Cn_beta_vert;

modelo.derivadas.Cn_beta    = Cn_beta;

%% CnT_beta Airplane Design (2088 PDF)

CnT_beta    = (pi/4)/(Sref*m22ft2)/(b_w*m2ft)*sum(((D_prop*m2ft).^2).*(X_Prop - Xcg).*dCN_dalpha);

modelo.derivadas.Cn_Tbeta   = CnT_beta;
%% ALEIRON DEFLECTION STABILITY DERIVATIVES

%% Cy_da = 0  AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 442, 2132 PDF)

Cy_da   = 0;

modelo.derivadas.Cy_da      = Cy_da;
%% Cl_da

%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 442, 2132 PDF)    
% Metodo 1
betaCldelta_kf  = betaCldelta_k_calc(y1_b2_w, TR_w, AR_w, LAMc4_w, Minf);
betaCldelta_k0  = betaCldelta_k_calc(y0_b2_w, TR_w, AR_w, LAMc4_w, Minf);
Cldelta_prima   = Cla_w/(2*pi*beta)*(betaCldelta_kf - betaCldelta_k0);

Cldelta         = Cldelta_Cldeltatheory_calc(cm_c_w, Cla_w)*Cldeltatheory_calc(cm_c_w, t_c_w);
alpha_delta     = Cldelta/Cla_w;

Cl_da           = Cldelta_prima*alpha_delta*S_w/Sref;

% Metodo 2 Diapositivas Calculo de Aeronaves 14.4

K1              = Kb_calc(y0_b2_w, y1_b2_w, TR_w);
K2              = alphaCL_alphaCl_calc(cm_c_w,AR_w);
Cl_da2           = K1*K2*Cldelta;
modelo.derivadas.Cl_da      = Cl_da;
%% Cn_da
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.5 (pag 448, 2138 PDF)

Ka      = Ka_calc(y0_b2_w,AR_w,TR_w);

Cn_da   = 2*Ka*CL*Cl_da;

modelo.derivadas.Cn_da      = Cn_da;

%% AILERON DEFLECTION STABILITY DERIVATIVES

%% Cy_dr
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)

% 1951
k_prima                 = 1; % Figura 8.13 Airplane Design (1918 PDF)
alphaCL_alphaCl         = alphaCL_alphaCl_calc(cm_c_v, AReff_v);
Cldelta_Cldeltatheory   = Cldelta_Cldeltatheory_calc(cm_c_v, Cla_v);
Cldelta_theory          = Cldeltatheory_calc(cm_c_v, t_c_v);
K_b                     = Kb_calc(y0_b2_v, y1_b2_v, TR_v);

Cy_dr   = K_b*CLa_v*S_v/Sref*eta_v*k_prima/Cla_v*Cldelta_Cldeltatheory*Cldelta_theory*alphaCL_alphaCl; % Corrección de la expresion dada en Airplane Design, mirar informe de erratas

if isnan(Cy_dr)
    Cy_dr = 0;
end

modelo.derivadas.Cy_dr      = Cy_dr;

%% Cl_dr
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 461, 2151 PDF)

Cl_dr   = (Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w*Cy_dr;

if isnan(Cl_dr)
    Cl_dr = 0;
end

modelo.derivadas.Cl_dr      = Cl_dr;

%% Cn_dr
%   - AIRPLANE DESIGN, PART VI; Roskam, Jan; 10.3.8 (pag 462, 2152 PDF)

Cn_dr   = -(Zca_v*sin(alpha) + (Xca_v - Xcg)*cos(alpha))/b_w*Cy_dr;

if isnan(Cn_dr)
    Cn_dr = 0;
end

modelo.derivadas.Cn_dr      = Cn_dr;

%% ROLL RATE STABILITY DERIVATIVES
%% Cy_p
% DATCOM 7.1.2.1 (pagina 2523)

% Vertical
Cy_p_vert   = 2/b_w*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)*Cy_beta_vert;

if isnan(Cy_p_vert)
    Cy_p_vert = 0;
end

% Ala
% Se han mezclado los métodos del DATCOM y los de Pamadi (pag 406)
aw1         = CLa_we/(pi*AR_we*eOswald_w);
K           = (1-aw1)/(1-eOswald_w*aw1);    % Pamadi Ecuacion 4.549
%K           = ((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - (k1*CLa + 2*k2*CL*CLa))/((CLa*tan(alpha) + CL*(1 + (tan(alpha))^2)) - 2*CL*CLa/pi/ARexp_w); 
Cyp_CL_M0   = Cyp_CL_M0_calc(LAMc4_w*180/pi, TR_w);
Cyp_CL_CL0  = (AR_we + 4*cos(LAMc4_w))/(AR_we*B_prandtl + 4*cos(LAMc4_w))*(AR_we*B_prandtl + cos(LAMc4_w))/(AR_we + cos(LAMc4_w))*Cyp_CL_M0;

k_Clp       = Cla_w/2/pi;   
betaClp_k   = C_lp_calc(TR_w, AR_w, LAMc4_w, Minf);
Clp_die0    = betaClp_k*k_Clp/beta;
DCy_p_die   = (3*sin(diedro_w)*(1 - 4*Zca_w/b_w*sin(diedro_w)))*Clp_die0;

Cy_p_w      = K*Cyp_CL_CL0*CL + DCy_p_die;

% DERIVADA TOTAL

Cy_p = Cy_p_vert + Cy_p_w;

modelo.derivadas.Cy_p       = Cy_p;

%% Cl_p
% Vertical
% Pamadi Ecuacion 4.579, Pagina 412

Cl_p_vert           = abs(2*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w)*(((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w))*Cy_beta_vert;

if isnan(Cy_p_vert)
    Cy_p_vert = 0;
end
% Ala
% Pamadi Ecuacion 4.576, Pagina 412
% Effect of drag on rolling moment is ignored
Clp_ClpDiedro0  = 1 - 4*Zca_w/b_w*sin(diedro_w) + 3*(2*Zca_w/b_w)^2*(sin(diedro_w))^2;
Cl_p_w          = betaClp_k*k_Clp/beta*Clp_ClpDiedro0;

% DERIVADA TOTAL

Cl_p            = Cl_p_vert + Cl_p_w;

modelo.derivadas.Cl_p       = Cl_p;

%% Cn_p
% Ala
 
Cnp_CL_CL0_M0   = -(AR_we + 6*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)*(xi/AR_we + tan(LAMc4_w)/12))/(6*(AR_we + 4*cos(LAMc4_w)));
Cnp_CL_CL0      = ((AR_we + 4*cos(LAMc4_w))/((AR_we*beta + 4*cos(LAMc4_w))))*...
                ((AR_we*beta + 0.5*((AR_we*beta + cos(LAMc4_w)))*(tan(LAMc4_w)^2))/...
                (AR_we + 0.5*(AR_we + cos(LAMc4_w))*tan(LAMc4_w)^2))*Cnp_CL_CL0_M0;

Cn_p_w          = Cl_p_w*tan(alpha)*(K-1) + K*Cnp_CL_CL0*CL; %Pamadi 4.594
% Vertical

Cn_p_vert       = -2/b_w*(Zca_v*sin(alpha) + (Xca_v - Xcg)*cos(alpha))*((Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha)) - Zca_v)/b_w*Cy_beta_vert;

if isnan(Cn_p_vert)
    Cy_p_vert = 0;
end

% DERIVADA TOTAL
Cn_p            = Cn_p_vert + Cn_p_w;

modelo.derivadas.Cn_p       = Cn_p;

%% YAW RATE STABILITY DERIVATIVES

%% Cy_r

Cy_r_w      = 0;
Cy_r_v      = -(2/b_w)*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))*Cy_beta_vert;

if isnan(Cy_r_v)
    Cy_r_v = 0;
end

Cy_r        = Cy_r_w + Cy_r_v;

modelo.derivadas.Cy_r       = Cy_r;

%% Cl_r

% Ala
Num_clr     = 1 + (AR_we*(1 - beta^2))/(2*beta*(AR_we*beta + 2*cos(LAMc4_w))) + ...
            ((AR_we*beta + 2*cos(LAMc4_w))/(AR_we*beta + 4*cos(LAMc4_w)))*...
            ((tan(LAMc4_w))^2)/8;
Den_clr     = 1 + ((AR_we + 2*cos(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)))*...
            ((tan(LAMc4_w))^2)/8;
Clr_CL_M0   = Clr_CL_M0_calc(AR_we, TR_we, LAMc4_w); % Pamadi Fig. 4.28
Clr_CL_CL0  = Num_clr/Den_clr*Clr_CL_M0; % Pamadi Ec. 4.614

DClr_diedro = 1/12*(pi*AR_we*sin(LAMc4_w))/(AR_we + 4*cos(LAMc4_w)); % Pamadi Ec. 4.617

Cl_r_w      = CL*Clr_CL_CL0 + DClr_diedro*diedro_w; % Pamadi Ec. 4.613


% Vertical
% Pamadi 4.618
Cl_r_vert   = -2/b_w^2*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))*(Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))*Cy_beta_vert;

if isnan(Cl_r_vert)
    Cl_r_vert = 0;
end
% DERIVADA TOTAL
Cl_r        = Cl_r_w + Cl_r_vert;

modelo.derivadas.Cl_r       = Cl_r;

%% Cn_r

% Ala
% Este metodo es una mierda, hay que cambiarlo por el indicado en Pamadi
% pagina 419
Cn_r_w      = -(1/3)*(CD0 + k1*CL + k2*CL^2);

% Vertical
Cn_r_vert   = (2/b_w^2)*(((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))^2)*Cy_beta_vert;

if isnan(Cn_r_vert)
    Cn_r_vert = 0;
end

% DERIVADA TOTAL
Cn_r = Cn_r_vert + Cn_r_w;

modelo.derivadas.Cn_r       = Cn_r;

%% SIDESLIP RATE STABILITY DERIVATIVES

% Pamadi pagina 419

%% Cy_betaDot
% DATCOM 7.4.4.4 (2825 PDF)
sigma_beta_alpha    = sigma_beta_alpha_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, TR_w);
sigma_beta_diedro   = sigma_beta_diedro_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, Minf);
sigma_beta_WB       = sign(Zca_w)*sigma_beta_wb_calc(2*Zca_v/b_w, 180/pi*LAMc4_w, TR_w, Wmax_fus/b_w);

sigma_beta          = sigma_beta_alpha*alpha*180/pi + sigma_beta_diedro*diedro_w + sigma_beta_WB;
Cy_betaDot          = 2*CLa_v*sigma_beta*S_v/Sref/b_w*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha)); 

if isnan(Cy_betaDot)
    Cy_betaDot = 0;
end

modelo.derivadas.Cy_betaDot     = Cy_betaDot;


%% Cl_betaDot

Cl_betaDot          = Cy_betaDot*(Zca_v*cos(alpha) - (Xca_v - Xcg)*sin(alpha))/b_w;

if isnan(Cl_betaDot)
    Cl_betaDot = 0;
end

modelo.derivadas.Cl_betaDot     = Cl_betaDot;

%% Cn_betaDot

Cn_betaDot          = -Cy_betaDot*((Xca_v - Xcg)*cos(alpha) + Zca_v*sin(alpha))/b_w;

if isnan(Cn_betaDot)
    Cn_betaDot = 0;
end

modelo.derivadas.Cn_betaDot     = Cn_betaDot;

end







