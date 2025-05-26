function Stab_Dyn_Long = longitudinal_analysis(V,rho,S_w,b_w,cmac_w,m_TO_1,Stab_Der,...
    q_inf,theta1,conv_UNITS,Stability_Pamadi,Weight_tier)


g = conv_UNITS.g;
in2m = conv_UNITS.in2m;

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%COEFICIENTES Stab_DerLIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ixx = Weight_tier.Ixx;
Iyy = Weight_tier.Iyy;
Izz = Weight_tier.Izz;
Ixz = Weight_tier.Ixz;

c_w = cmac_w;
cmed = c_w/(2*V);                  %adimensionalización de la cuerda
m1 = 2*m_TO_1/(rho*S_w*V);               %adimensionalización de la masa

CD_alpha = Stab_Der.CD_alpha;
CD = Stab_Der.CD;

CXalfa = Stab_Der.CXalfa;
CZalfa = Stab_Der.CZalfa;
CMalfa = Stab_Der.CMalfa;

CDu = Stab_Der.CDu;
CXu = Stab_Der.CXu;
CLu = Stab_Der.CLu;
CZu = Stab_Der.CZu;
CMu = Stab_Der.CMu;

% CLq_w2 = Stab_Der.CLq_h;
% CLq_wb = Stab_Der.CLq_wb;
% CMq_h = Stab_Der.CMq_h;
% CMq_wb = Stab_Der.CMq_wb;
CXq = Stab_Der.CXq;
CLq = Stab_Der.CLq;
CZq = Stab_Der.CZq;
CMq = Stab_Der.CMq;

CDq = -CXq;
% Stab_Der.CL_alphapunto_t = Stab_Der.CL_alphapunto_t;
% CL_alphapunto_WB = Stab_Der.CL_alphapunto_WB;
% CM_alphapunto_t = Stab_Der.CM_alphapunto_t;
% CM_alphapunto_WB = Stab_Der.CM_alphapunto_WB;
CLalphapunto = Stab_Der.CLalphapunto;
CXalfapunto = Stab_Der.CXalfapunto;
CZalfapunto = Stab_Der.CZalfapunto;
CMalfapunto = Stab_Der.CMalphapunto;

CXteta = Stab_Der.CXteta;
CZteta = Stab_Der.CZteta;
CMteta = Stab_Der.CMteta;

CXdeltae = Stab_Der.CXdeltae;
CDdeltae = Stab_Der.CDdeltae;
CZdeltae = Stab_Der.CZdeltae;
CMdeltae = Stab_Der.CMdeltae;
CL_delta_e = Stab_Der.CL_delta_e;

CTx1 = Stab_Der.CTx1;
CTxu = Stab_Der.CTxu;
CTxalpha = Stab_Der.CTxalpha;
CMt1 = Stab_Der.CMt1;
CMtu = Stab_Der.CMtu;
CMtalpha = Stab_Der.CMtalpha;

% CTx1 = 0;
% CTxu = 0;
% CTxalpha = 0;
% CMt1 = 0;
% CMtu = 0;
% CMtalpha = 0;

CL = Stab_Der.CL;
CD = Stab_Der.CD;
CM = Stab_Der.CM;

CL_alpha_w1w2b = Stab_Der.CL_alpha_w1w2b;
CM_alpha_w1w2b = Stab_Der.CM_alpha_w1w2b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CONSTRUCCIÓN DE LA MATRIZ A y B%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inercias de catia
%inercias adimensionales
Ix1=Ixx/(0.5*rho*V^2*S_w*b_w); 
Iz1=Izz/(0.5*rho*V^2*S_w*b_w); 
Ixz1=Ixz/(0.5*rho*V^2*S_w*b_w); 
Iy1=Iyy/(0.5*rho*V^2*S_w*c_w);       %adimensionalización de la inercia

%inercias de trabajo
Ix1prima=Ix1/(Ix1*Iz1-Ixz1^2);
Iz1prima=Iz1/(Ix1*Iz1-Ixz1^2);
Ixz1prima=Ixz1/(Ix1*Iz1-Ixz1^2);


% Fran Analysis
if Stability_Pamadi == 0
    
   
    
    mu = 2*m_TO_1/(rho*c_w*S_w); % Masa adimensional
    I_y_gorro = Ixz/(rho*S_w*(c_w/2)^3); % Inercia adimensional

     
    % SISTEMA LONGITUDINAL COMPLETO
    
    
        
    matriz_masa = [2*mu 0 0 0;
        0 2*mu-CZalfapunto 0 0;
        0 -CMalfapunto I_y_gorro 0;
        0 0 0 1];
    
    C_xs = 0; % Para vuelo de crucero. Corregir si se analiza ascenso/descenso
    C_zs = -(2*g*m_TO_1) / (rho*S_w*V^2); % Coeficiente de sustentación de vuelo
    matriz_A = [2*C_xs+CXu CXalfa 0 C_zs;
        CZu CZalfa 2*mu+CZq -C_xs;
        CMu CMalfa CMq 0;
        0 0 1 0];
    
    matriz_B = [0;
        CZdeltae;
        CMdeltae;
        0];     
    
    A_long = matriz_masa\matriz_A;
    B_long = matriz_masa\matriz_B;
    
    
    poles = eig(A_long)*2*V/c_w; % Hay que dimensionalizar
%         poles = eig(A_long); % Hay que dimensionalizar
    k=length(poles);
    
elseif Stability_Pamadi == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%CONSTRUCCIÓN DE LA MATRIZ A y B%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %se procede a colocar cada coeficiente por su sitio dentro de la matriz
    %además de mu, Iy1 y cmed se necesitan otros coeficientes adim
    
    chi1=CXalfapunto*cmed/(m1-CZalfapunto*cmed);
    chi2=CMalfapunto*cmed/(m1-CZalfapunto*cmed);
    
    %los coeficientes son por orden:
    % a11=((Cxu+Ctxu) + chi1*Czu)/(m1);
    a11=((CXu) + chi1*CZu)/(m1);
    a12=(CXalfa+chi1*CZalfa)/(m1);
    a13=(CXq*cmed+chi1*(m1+CZq*cmed))/(m1);
    a14=(CXteta+chi1*CZteta)/(m1);
    a21=(CZu)/(m1-CZalfapunto*cmed);
    a22=(CZalfa)/(m1-CZalfapunto*cmed);
    a23=(m1+CZq*cmed)/(m1-CZalfapunto*cmed);
    a24=(CZteta)/(m1-CZalfapunto*cmed);
    % a31=((Cmu+Cmtu)+chi2*Czu)/(Iy1);
    a31=((CMu)+chi2*CZu)/(Iy1);
    % a32=((Cmalfa+Cmtalpha)+chi2*Czalfa)/(Iy1);
    a32=((CMalfa)+chi2*CZalfa)/(Iy1);
    a33=(CMq*cmed+chi2*(m1+CZq*cmed))/(Iy1);
    a34=(chi2*CZteta)/(Iy1);
    a41=0;
    a42=0;
    a43=1;
    a44=0;
    
%     A_long_pamadi=[a11 a12 a13 a14;a21 a22 a23 a24;a31 a32 a33 a34;a41 a42 a43 a44];
    A_long=[a11 a12 a13 a14;a21 a22 a23 a24;a31 a32 a33 a34;a41 a42 a43 a44];
    
    %para la matriz B los coeficientes son:
    b1=(CXdeltae+chi1*CZdeltae)/(m1);
    b2=(CZdeltae)/(m1-cmed*CZalfapunto);
    b3=(CMdeltae+chi2*CZdeltae)/(Iy1);
    b4=0;
    
    B_long=[b1;b2;b3;b4];
    
elseif Stability_Pamadi == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Longitudinal Dimensional Stability Derivatives %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    q1barw = q_inf;
    Sw = S_w;
    m = m_TO_1;
    u1 = V;
    cbar = cmac_w;
    cm1 = CM;
    
    Xu=(-q1barw*Sw*(CDu+(2*CD)))/(m*u1)
    %     Xu=(q1barw*Sw*CXu)/(m*u1)
    Xtu=(q1barw*Sw*(CTxu+(2*CTx1)))/(m*u1)
    Xalpha=(-q1barw*Sw*(CD_alpha-CL))/m
    %     Xalpha=(q1barw*Sw*CXalfa)/m
    Xq=(-q1barw*Sw*CDq)/(2*m*u1)
    
    Xdeltae=(-q1barw*Sw*CDdeltae)/m
    
    Zu=(-q1barw*Sw*(CLu+(2*CL)))/(m*u1)
    %     Zu=(q1barw*Sw*CZu)/(m*u1)
    Zalpha=(-q1barw*Sw*(CL_alpha_w1w2b + CD))/m
    %     Zalpha=(q1barw*Sw*CZalfa)/m
    Zalphadot=(-q1barw*Sw*CLalphapunto*cbar)/(2*m*u1)
    %     Zalphadot=(q1barw*Sw*CZalfapunto*cbar)/(2*m*u1)
    Zq=(-q1barw*Sw*CLq*cbar)/(2*m*u1)
    
   
    Zdeltae=(-q1barw*Sw*CL_delta_e)/m
    %     Zdeltae=(q1barw*Sw*CZdeltae)/m
    Mu=(q1barw*Sw*cbar*(CMu+(2*CM)))/(Iyy*u1)
    Mtu=(q1barw*Sw*cbar*(CMtu+(2*CMt1)))/(Iyy*u1)
    Malpha=(q1barw*Sw*cbar*CM_alpha_w1w2b)/Iyy
    Malphadot=(q1barw*Sw*(cbar^2)*CMalfapunto)/(2*Iyy*u1)
    Mtalpha=(q1barw*Sw*cbar*CMtalpha)/Iyy
    Mq=(q1barw*Sw*(cbar^2)*CMq)/(2*Iyy*u1)
    Mdeltae=(q1barw*Sw*cbar*CMdeltae)/Iyy
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % terms of the Matrix that finds the poles of the system response %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    aa = Xu + Xtu;
    bb = Xalpha;
    cc = Xq;
    dd = -g*cos(theta1*D2R);
    
    ee = Zu/(u1 - Zalphadot);
    ff = Zalpha/(u1 - Zalphadot);
    gg = (Zq + u1)/(u1 - Zalphadot);
    hh = (-g*sin(theta1*D2R))/(u1 - Zalphadot);
    
    ii = Mu + Mtu + ((Zu*Malphadot)/(u1 - Zalphadot));
    jj = Malpha + Mtalpha +((Zalpha*Malphadot)/(u1 - Zalphadot));
    kk = Mq + (((Zq+u1)*Malphadot)/(u1 - Zalphadot));
    ll = -(g*sin(theta1*D2R)*Malphadot)/(u1 - Zalphadot);

    ee_r = Zu;
    ff_r = Zalpha;
    gg_r = Zq + u1;
    hh_r = -g*sin(theta1*D2R);
    
    ii_r = Mu + Mtu;
    jj_r = Malpha + Mtalpha;
    kk_r = Mq;
    ll_r = 0;

    mm = 0;
    nn = 0;
    oo = 1;
    pp = 0;
   
    M_long = [1 0 0 0; 0 (u1 - Zalphadot) 0 0; 0 -Malphadot 1 0; 0 0 0 1];
    Ra_long = [aa bb cc dd; ee_r ff_r gg_r hh_r; ii_r jj_r kk_r ll_r; mm nn oo pp];
    
    A_long = inv(M_long)*Ra_long;
%     A_long = [aa bb cc dd; ee ff gg hh; ii jj kk ll; mm nn oo pp];

    qq = Xdeltae;
    rr = (Zdeltae/(u1 - Zalphadot));
    ss = Mdeltae + ((Zdeltae*Malphadot)/(u1 - Zalphadot));
    tt = 0;

    qq_r = Xdeltae;
    rr_r = Zdeltae;
    ss_r = Mdeltae;
    tt_r = 0;

    Rb_long=[qq_r; rr_r; ss_r; tt_r];
    B_long = inv(M_long)*Rb_long;
%     B_long=[qq; rr; ss; tt];
    
end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %solve for the eigenvalues of the system response for the longitudinal %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [vector,poleslon] = eig(Matrix_lon);

        poles = eig(A_long); % Hay que dimensionalizar
    k=length(poles);

%checks which set of poles are the shorp period and which one the Phugoid
pole1=poles(1);
pole1a=poles(2);
pole2=poles(3);
pole2a=poles(4);

Stab_Dyn_Long.pole1 = pole1;
Stab_Dyn_Long.pole2 = pole1a;
Stab_Dyn_Long.pole3 = pole2;
Stab_Dyn_Long.pole4 = pole2a;


if abs(pole1)>abs(pole2)
    SP_poles=poles(1);
    phugoid_poles=poles(3);
else
    SP_poles=poles(3);
    phugoid_poles=poles(1);
end

% Calculation of the Natural Frequency Damping ratio and 
% time constants for the Short Period Mode

% real and imaginary parts
re_sp=real(SP_poles);
im_sp=imag(SP_poles);

Stab_Dyn_Long.SP_poles = SP_poles;
Stab_Dyn_Long.phugoid_poles = phugoid_poles;


% Stab_Dyn_Long.re_sp = re_sp;
% Stab_Dyn_Long.im_sp = im_sp;

% damping ratio
ZETA_SP=abs(re_sp/sqrt(re_sp^2 + im_sp^2));
% natural frequency
WN_SP=sqrt(re_sp^2 + im_sp^2);

% WN_SP_approx=sqrt(Zalpha*Mq/u1 - Malpha)
% ZETA_SP_approx = -(Mq + (Zalpha/u1) + Malphadot)/(2*WN_SP_approx)

Stab_Dyn_Long.ZETA_SP = ZETA_SP;
% Stab_Dyn_Long.ZETA_SP_approx = ZETA_SP_approx;
Stab_Dyn_Long.WN_SP = WN_SP;
% Stab_Dyn_Long.WN_SP_approx = WN_SP_approx;

% period
T_SP=(2*pi)/(WN_SP*sqrt(1-ZETA_SP^2));
%time to half or double
ta_SP=0.693/abs(re_sp);

Stab_Dyn_Long.T_SP = T_SP;
Stab_Dyn_Long.ta_SP = ta_SP;

% Calculation of the Natural Frequency Damping ratio and 
% time constants for the Phugoid Mode
re_ph=real(phugoid_poles);
im_ph=imag(phugoid_poles);

Stab_Dyn_Long.re_ph = re_ph;
Stab_Dyn_Long.im_ph = im_ph;

% natural frequency
WN_PH=sqrt(re_ph^2 + im_ph^2);
% WN_PH_approx=-sqrt(-g*Zu/u1)
% damping ratio
ZETA_PH=abs(re_ph/sqrt(re_ph^2 + im_ph^2));
% ZETA_PH_approx= -Xu/(2*WN_PH_approx)

Stab_Dyn_Long.ZETA_PH = ZETA_PH;
Stab_Dyn_Long.WN_PH = WN_PH;

% period
T_PH=(2*pi)/(WN_PH*sqrt(1-ZETA_PH^2));
%time to half or double
ta_PH=0.693/abs(re_ph);

Stab_Dyn_Long.T_PH = T_PH;
Stab_Dyn_Long.ta_PH = ta_PH;

% Forwards Speed Stability
crit1 = CTxu-CDu;
if crit1 < 0
    Forward_Speed_Stability =1;
else
    Forward_Speed_Stability =0;
end

% Vertical Speed Stability
crit3 = - CL_alpha_w1w2b;
if crit3 < 0
    Vertical_Speed_Stability =1;
else
    Vertical_Speed_Stability =0;
end

% Angle of Attack Stability 
crit4 = CM_alpha_w1w2b + CMtalpha;
if crit4 < 0
    AoA_Stability =1;
else
    AoA_Stability =0;
end

% Pitch Rate Stablity
crit7 = CMq;
if crit7 < 0
    Q_Stability =1;
else
    Q_Stability =0;
end

% Forward Speed in Pitching Moment
crit9 = CMu + CMtu;
if crit9 > 0
    u_Q_Stability =1;
else
    u_Q_Stability =0;
end

Stab_Dyn_Long.Forward_Speed_Stability = Forward_Speed_Stability;
Stab_Dyn_Long.Vertical_Speed_Stability = Vertical_Speed_Stability;
Stab_Dyn_Long.AoA_Stability = AoA_Stability;
Stab_Dyn_Long.Q_Stability = Q_Stability;
Stab_Dyn_Long.u_Q_Stability = u_Q_Stability;
