function results = dynamic_analysis(model)
g       = 9.8065;

u1      = model.general.Vinf;
q1      = model.general.qinf;
Sref    = model.general.Sref;
m       = model.general.mtow*model.general.w_w0;
%W       = m*g;
% theta1  = model.general.alpha;
theta1  = 0;

Ixx     = model.general.Ixx;
Iyy     = model.general.Iyy;
Izz     = model.general.Izz;
Ixz     = model.general.Ixz;

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
CMt1    = 0;
CMtu    = 0;
% CMta    = model.derivadas.CM_Ta;
CMta    = 0;


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
try
Ra_long     =   [Xu+Xtu,    Xa,         Xq,     -g*cos(theta1);
                Zu,         Za,         Zq+u1,  -g*sin(theta1);
                Mu+Mtu,     Ma+Mta      Mq,     0;
                0,          0,          1,      0];

Ma_long     =   [1,         0,          0,      0;
                1,          u1-ZaDot,   0,      0;
                0,          -MaDot,     1,      0;
                0,          0,          0,      1];

A_long      =   Ma_long\Ra_long;

long_poles = eig(A_long); % Hay que dimensionalizar


%checks which set of poles are the shorp period and which one the Phugoid
pole1=long_poles(1);
pole1a=long_poles(2);
pole2=long_poles(3);
pole2a=long_poles(4);

results.long.pole1 = pole1;
results.long.pole2 = pole1a;
results.long.pole3 = pole2;
results.long.pole4 = pole2a;


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

results.long.SP_poles = SP_poles;
results.long.PH_poles = phugoid_poles;

% damping ratio
results.long.SP_damp    = abs(re_sp/sqrt(re_sp^2 + im_sp^2));
% natural frequency
results.long.SP_wn      = abs(SP_poles);
DAMP_SP                 = results.long.SP_damp;
WN_SP                   = results.long.SP_wn;

 
% period
T_SP=(2*pi)/(WN_SP*sqrt(1-DAMP_SP^2));
%time to half or double
ta_SP=0.693/abs(re_sp);

results.long.SP_T    = T_SP;
results.long.SP_T2   = ta_SP;

% Calculation of the Natural Frequency Damping ratio and
% time constants for the Phugoid Mode
re_ph=real(phugoid_poles);
im_ph=imag(phugoid_poles);

% natural frequency

WN_PH=abs(phugoid_poles);
% WN_PH_approx=-sqrt(-g*Zu/u1)
% damping ratio
ZETA_PH=abs(re_ph/sqrt(re_ph^2 + im_ph^2));
ZETA_PH=abs(re_ph/WN_PH);
%ZETA_PH_approx= -Xu/(2*WN_PH_approx)

results.long.PH_damp    = ZETA_PH;
results.long.PH_wn      = WN_PH;

% period
T_PH=(2*pi)/(WN_PH*sqrt(1-ZETA_PH^2));
%time to half or double
ta_PH=0.693/abs(re_ph);

results.long.PH_T   = T_PH;
results.long.PH_T2  = ta_PH;

% % Forwards Speed Stability
% crit1 = CTxu-CDu;
% if crit1 < 0
%     Forward_Speed_Stability =1;
% else
%     Forward_Speed_Stability =0;
% end
% 
% % Vertical Speed Stability
% crit3 = - CL_alpha_wbh;
% if crit3 < 0
%     Vertical_Speed_Stability =1;
% else
%     Vertical_Speed_Stability =0;
% end
% 
% % Angle of Attack Stability
% crit4 = CM_alpha_wbh + CMtalpha;
% if crit4 < 0
%     AoA_Stability =1;
% else
%     AoA_Stability =0;
% end
% 
% % Pitch Rate Stablity
% crit7 = CMq;
% if crit7 < 0
%     Q_Stability =1;
% else
%     Q_Stability =0;
% end
% 
% % Forward Speed in Pitching Moment
% crit9 = CMu + CMtu;
% if crit9 > 0
%     u_Q_Stability =1;
% else
%     u_Q_Stability =0;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%% ESTABILIDAD LATERAL-DIRECCIONAL %%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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
% 
% CyTb    = model.derivadas.Cy_Tbeta;
% CnTb    = model.derivadas.Cn_Tbeta;

CyTb    = 0;
CnTb    = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Lateral-Directional Dimensional Stability Derivatives %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ybeta   = (q1*Sref*Cyb)/m;
Yp      = (q1*Sref*b_w*Cyp)/(2*m*u1);
Yr      = (q1*Sref*b_w*Cyr)/(2*m*u1);
   
Lbeta   = (q1*Sref*b_w*Clb)/Ixx;
Lp      = (q1*Sref*b_w^2*Clp)/(2*Ixx*u1);
Lr      = (q1*Sref*b_w^2*Clr)/(2*Ixx*u1);

Nbeta   = (q1*Sref*b_w*Cnb)/Izz;
Ntbeta  = (q1*Sref*b_w*CnTb)/Izz;
Np      = (q1*Sref*b_w^2*Cnp)/(2*Izz*u1);
Nr      = (q1*Sref*b_w^2*Cnr)/(2*Izz*u1);


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

A_lat       =   Ma_lat\Ra_lat;

lat_poles   = eig(A_lat);

%checks which set of poles are the shorp period and which one the Phugoid
pole1=lat_poles(1);
pole2=lat_poles(2);
pole3=lat_poles(3);
pole4=lat_poles(4);
pole5=lat_poles(5);

results.lat.pole1 = pole1;
results.lat.pole2 = pole2;
results.lat.pole3 = pole3;
results.lat.pole4 = pole4;
results.lat.pole5 = pole5;

real_poles = real(lat_poles);

pole_check = 0;

for ii=1:length(real_poles)
    for jj=1:length(real_poles)
        if real_poles(ii) == real_poles(jj)
            check_ij(ii,jj)=1;
        else
            check_ij(ii,jj)=0;
        end
    end
end

% for ii=1:length(real_poles)
%     same_poles(ii) = sum(check_ij(ii,:));
% end

same_poles = sum(check_ij');

ii=0;
jj=0;
k=0;
for ii=1:length(lat_poles)
    if same_poles(ii) == 2
        jj=jj+1;
        DR_poles(jj) = lat_poles(ii);
    else 
        k=k+1;
        Rest_poles(k)= lat_poles(ii);
    end
end

for ii=1:length(Rest_poles)
    Rest_poles_vec(ii) = real(Rest_poles(ii));
end

[B,I] = sort(abs(Rest_poles_vec));

Yaw_pole = Rest_poles(I(1));
Spiral_pole = Rest_poles(I(2));
Roll_pole = Rest_poles(I(3));

%Stab_Dyn_LatDir.DR_poles = DR_poles(1);
results.lat.Spiral_pole = Spiral_pole;
results.lat.Roll_pole = Roll_pole;
results.lat.Yaw_pole = Yaw_pole;

% real and imaginary parts
re_dr=real(DR_poles(1));
im_dr=imag(DR_poles(1));

% damping ratio
ZETA_DR=abs(re_dr/sqrt(re_dr^2 + im_dr^2));
% natural frequency
WN_DR=sqrt(re_dr^2 + im_dr^2);

results.lat.DR_damp     = ZETA_DR;
results.lat.DR_wn       = WN_DR;

% period
T_DR=(2*pi)/(WN_DR*sqrt(1-ZETA_DR^2));
%time to half or double
ta_DR=0.693/abs(re_dr);

results.lat.DR_T    = T_DR;
results.lat.DR_T2   = ta_DR;

results.lat.ROL_T2  = 0.693/abs(Roll_pole);
results.lat.ESP_T2  = 0.693/abs(Spiral_pole);

catch me
    if (strcmp(me.identifier,'MATLAB:catenate:dimensionMismatch'))
        msg = 'Stability Derivatives have to be calculated and saved before perform any dynamic stability analysis';
        causeException = MException('MATLAB:myCode:dimensions',msg);
        me = addCause(me,causeException);
    end
    rethrow(me)
    %disp(me);
end

% if Roll_pole < 0
%     ta_half_roll = 0.693/abs(Roll_pole);
%     Stab_Dyn_LatDir.ta_half_roll = ta_half_roll;
% else
%     ta_double_roll = 0.693/abs(Roll_pole);
%     Stab_Dyn_LatDir.ta_double_roll = ta_double_roll;
% end
% 
% if Spiral_pole < 0
%     ta_half_spiral = 0.693/abs(Spiral_pole);
%     Stab_Dyn_LatDir.ta_half_spiral = ta_half_spiral;
% else
%     ta_double_spiral = 0.693/abs(Spiral_pole);
%     Stab_Dyn_LatDir.ta_double_spiral = ta_double_spiral;
% end

% if Stability_Pamadi == 0
%     % Aproximations
%     WN_DR_approx = sqrt(Nbeta + (1/V)*(Ybeta*Nr - Nbeta*Yr))
%     ZETA_DR_approx = - (Nr + Ybeta/V)/(2*WN_DR)
%     
%     Spiral_approx = (Lbeta*Nr-Nbeta*Lr)/(Lbeta + Nbeta*A1)
% end

% % Side Speed Stability
% crit2 = Cyb - CyTb;
% if crit2 < 0
%     Sideslip_Speed_Stability =1;
% else
%     Sideslip_Speed_Stability =0;
% end
% 
% % Angle of Sideslip Stability 
% crit5 = Cnb + CNTb;
% if crit5 > 0
%     Beta_Stability =1;
% else
%     Beta_Stability =0;
% end
% 
% % Roll Rate Stablity
% crit6 = Clp;
% if crit6 < 0
%     P_Stability =1;
% else
%     P_Stability =0;
% end
% 
% % Yaw Rate Stablity
% crit8 = Cnr;
% if crit8 < 0
%     R_Stability =1;
% else
%     R_Stability =0;
% end
% 
% % Sideslip on Rolling Moment Stablity
% crit10 = Clb;
% if crit10 < 0
%     beta_L_Stability =1;
% else
%     beta_L_Stability =0;
% end
% 
% % Sideslip on Yawing Momwnt Stablity
% crit11 = Cnb;
% if crit11 > 0
%     beta_N_Stability =1;
% else
%     beta_N_Stability =0;
% end
% 
% % Spiral Root Stablity
% crit11 = (Lbeta*Nr-Nbeta*Lr);
% if crit11 > 0
%     spiral_Stability =1;
% else
%     spiral_Stability =0;
% 
% 
