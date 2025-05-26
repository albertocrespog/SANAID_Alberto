function model = dynamic_analisis(model)
g       = 9.8065;

u1      = model.general.Vinf;
q1      = model.general.qinf;
Sref    = model.general.Sref;
CD0     = model.general.polar(1);
k1      = model.general.polar(2);
k2      = model.general.polar(3);
W       = model.general.W;
m       = W/g;
theta1  = model.general.alpha;

Ixx     = model.pesos.Ixx;
Iyy     = model.pesos.Iyy;
Izz     = model.pesos.Izz;
Ixz     = model.pesos.Ixz;

c_w     = model.ala.MAC;
b_w     = model.ala.b;

CL      = model.derivadas.CL;
CD      = model.derivadas.CD;
CM      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% ESTABILIDAD LONGITUDINAL %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
CLa     = model.derivadas.CL_a;
CDa     = model.derivadas.CD_a;
CMa     = model.derivadas.CM_a;

CLq     = model.derivadas.CL_q;
CDq     = model.derivadas.CD_q;
CMq     = model.derivadas.CM_q;

CLu     = model.derivadas.CL_u;
CDu     = model.derivadas.CD_u; 
CMu     = model.derivadas.CM_u;

CLaDot  = model.derivadas.CL_alphaDot;
CDaDot  = model.derivadas.CD_alphaDot;
CMaDot  = model.derivadas.CM_alphaDot;

CTx1    = model.derivadas.CT_x1;
CTxu    = model.derivadas.CT_xu;
CTxa    = model.derivadas.CT_xa;
CMt1    = model.derivadas.CM_T1;
CMtu    = model.derivadas.CM_Tu;
CMta    = model.derivadas.CM_Ta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CONSTRUCCIÓN DE LA MATRIZ A y B%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Longitudinal Dimensional Stability Derivatives %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xu      = (-q1*Sref*(CDu+(2*CD)))/(m*u1);
Xtu     = (q1*Sref*(CTxu+(2*CTx1)))/(m*u1);
Xa      = (-q1*Sref*(CDa-CL))/m;
Xq      = (-q1*Sref*CDq)/(2*m*u1); %% Revisar esta derivada
Xq      = 0;

Zu      = (-q1*Sref*(CLu+(2*CL)))/(m*u1);
Za      = (-q1*Sref*(CLa + CD))/m;
ZaDot   = (-q1*Sref*CLaDot*c_w)/(2*m*u1);
Zq      = (-q1*Sref*CLq*c_w)/(2*m*u1);

Mu      = (q1*Sref*c_w*(CMu+(2*CM)))/(Iyy*u1);
Mtu     = (q1*Sref*c_w*(CMtu+(2*CMt1)))/(Iyy*u1);
Ma      = (q1*Sref*c_w*CMa)/Iyy;
MaDot   = (q1*Sref*(c_w^2)*CMaDot)/(2*Iyy*u1);
Mta     = (q1*Sref*c_w*CMta)/Iyy;
Mq      = (q1*Sref*(c_w^2)*CMq)/(2*Iyy*u1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms of the Matrix that finds the poles of the system response %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ra_long     =   [Xu+Xtu,    Xa,         Xq,     -g*cos(theta1);
                Zu,         Za,         Zq+u1,  -g*sin(theta1);
                Mu+Mtu,     Ma+Mta      Mq,     0;
                0,          0,          1,      0];

Ma_long     =   [1,         0,          0,      0;
                1,          u1-ZaDot,   0,      0;
                0,          -MaDot,     1,      0;
                0,          0,          0,      1];

A_long      =   inv(Ma_long)*Ra_long;

long_poles = eig(A_long); % Hay que dimensionalizar


%checks which set of poles are the shorp period and which one the Phugoid
pole1=long_poles(1);
pole1a=long_poles(2);
pole2=long_poles(3);
pole2a=long_poles(4);

Stab_Dyn_Long.pole1 = pole1;
Stab_Dyn_Long.pole2 = pole1a;
Stab_Dyn_Long.pole3 = pole2;
Stab_Dyn_Long.pole4 = pole2a;


if abs(pole1)>abs(pole2)
    SP_poles=long_poles(1);
    phugoid_poles=long_poles(3);
else
    SP_poles=long_poles(3);
    phugoid_poles=long_poles(1);
end

% Calculation of the Natural Frequency Damping ratio and
% time constants for the Short Period Mode

% real and imaginary parts
re_sp=real(SP_poles);
im_sp=imag(SP_poles);

Stab_Dyn_Long.SP_poles = SP_poles;
Stab_Dyn_Long.phugoid_poles = phugoid_poles;

% damping ratio
ZETA_SP=abs(re_sp/sqrt(re_sp^2 + im_sp^2));
% natural frequency
WN_SP=sqrt(re_sp^2 + im_sp^2);


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
crit3 = - CL_alpha_wbh;
if crit3 < 0
    Vertical_Speed_Stability =1;
else
    Vertical_Speed_Stability =0;
end

% Angle of Attack Stability
crit4 = CM_alpha_wbh + CMtalpha;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ESTABILIDAD LATERAL-DIRECCIONAL %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cyb     = model.derivadas.Cy_beta; 
Clb     = model.derivadas.Cl_beta;
Cnb     = model.derivadas.Cn_beta;

Cyp     = model.derivadas.Cy_p;
Clp     = model.derivadas.Cl_p;
Cnp     = model.derivadas.Cn_p;

Cyr     = model.derivadas.Cy_r;
Clr     = model.derivadas.Cl_r;
Cnr     = model.derivadas.Cn_r;

CybDot  = model.derivadas.Cy_betaDot;
ClbDot  = model.derivadas.Cl_betaDot;
CnbDot  = model.derivadas.Cn_betaDot;

CyTb    = model.derivadas.Cy_Tbeta;
CnTb    = model.derivadas.Cn_Tbeta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Lateral-Directional Dimensional Stability Derivatives %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ybeta   = (q1*S_w*Cyb)/m;
Yp      = (q1*S_w*b_w*Cyp)/(2*m*u1);
Yr      = (q1*S_w*b_w*Cyr)/(2*m*u1);
   
Lbeta   = (q1*S_w*b_w*Clb)/Ixx;
Lp      = (q1*S_w*b_w^2*Clp)/(2*Ixx*u1);
Lr      = (q1*S_w*b_w^2*Clr)/(2*Ixx*u1);

Nbeta   = (q1*S_w*b_w*Cnb)/Izz;
Ntbeta  = (q1*S_w*b_w*CnTb)/Izz;
Np      = (q1*S_w*b_w^2*Cnp)/(2*Izz*u1);
Nr      = (q1*S_w*b_w^2*Cnr)/(2*Izz*u1);


Ra_lat      =   [Ybeta/u1,          Yp,         Yr-u1,          -g*cos(theta1),     0;
                Lbeta/u1,           Lp,         Lr,             0,                  0;
                (Nbeta+Ntbeta)/u1,  Np,         Nr,             0,                  0;
                0                   1,          tan(theta1),    0,                  0;   
                0,                  0,          sec(theta1),    0,                  0];

Ma_lat      =   [1,         0,          0,          0,      0;
                0,          1,          -Ixz/Ixx,   0,      0;
                0,          -Ixz/Izz,   1,          0,      0;
                0,          0,          0,          1,      0;
                0,          0,          0,          0,      1];

A_lat       =   inv(Ma_lat)*Ra_lat;

lat_poles   = eig(A_lat);


