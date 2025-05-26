clear all
close all
model=load('model.mat');
g       = 9.8065;

Delta = 1;
u1      = model.general.Vinf*Delta;
q1      = model.general.qinf*Delta^2;
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

CL_de = model.derivadas.CL_de;
CD_de = model.derivadas.CD_de;
CM_de = model.derivadas.CM_de;
% model.derivadas.CL_dc
% model.derivadas.CD_dc
% model.derivadas.CM_dc

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
% Xq      = 0;

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

Xdeltae = (-q1*Sref*CD_de)/m;
Zdeltae = (-q1*Sref*CL_de)/m;
Mdeltae = (q1*Sref*c_w*CM_de)/Iyy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms of the Matrix that finds the poles of the system response %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
% Xu+Xtu,    Xa,         Xq,     -g*cos(theta1)
% Zu,         Za,         Zq+u1,  -g*sin(theta1)
% Mu+Mtu,     Ma+Mta      Mq,     0
Ra_long     =   [Xu+Xtu,    Xa,         Xq,     -g*cos(theta1);
                Zu,         Za,         Zq+u1,  -g*sin(theta1);
                Mu+Mtu,     Ma+Mta      Mq,     0;
                0,          0,          1,      0];

Ma_long     =   [1,         0,          0,      0;
                1,          u1-ZaDot,   0,      0;
                0,          -MaDot,     1,      0;
                0,          0,          0,      1];

Matrix_lon      =   Ma_long\Ra_long;

qq = Xdeltae;
rr = (Zdeltae/(u1 - ZaDot));
ss = Mdeltae + ((Zdeltae*MaDot)/(u1 - ZaDot));
tt = 0;

Deflection_lon=[qq; rr; ss; tt];

long_poles = eig(Matrix_lon); % Hay que dimensionalizar

save Longitudinal.mat  Matrix_lon Deflection_lon theta1 u1

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

%% Short Period
WN_SP = abs(sqrt(re_sp^2 + im_sp^2));
WN_SP_approx=sqrt((Za*Mq/u1) - Ma);
DAMP_SP = abs(re_sp/WN_SP);
DAMP_SP_approx= - (Mq + (Za/u1) + MaDot)/(2*WN_SP);
re_sp_approx = -WN_SP_approx*DAMP_SP_approx;
im_sp_approx = WN_SP_approx*sqrt(1 - DAMP_SP_approx^2);

% damping ratio
results.long.SP_wn      = WN_SP;
results.long.SP_damp    = DAMP_SP;
 
% period
T_SP=(2*pi)/(WN_SP*sqrt(1-DAMP_SP^2));
%time to half or double
ta_SP=0.693/abs(re_sp);
% cycles to double or half
cycles_SP = 0.110*(sqrt(1-DAMP_SP^2))/abs(DAMP_SP);  
% Logarithmic decrement
decrement_SP = 0.693/cycles_SP;

results.long.SP_T    = T_SP;
results.long.SP_T2   = ta_SP;
results.long.cycles_SP    = cycles_SP;
results.long.decrement_SP   = decrement_SP;

% Eigenvalues
SP_msg = strcat('The Short Period (SP) Analysis');
SP_msg2 = strcat('------------------------------');
poles_sp_msg = strcat('The SP Poles are: real part: ',num2str(re_sp),' & imaginary part: ',num2str(im_sp),'i');
poles_sp_approx_msg = strcat('The SP Poles aproximation are: real part: ',num2str(re_sp_approx),' & imaginary part: ',num2str(im_sp_approx),'i');
WN_sp_msg = strcat('The SP natural frequency is: ',num2str(WN_SP),' and the approximation: ',num2str(WN_SP_approx));
DAMP_sp_msg = strcat('The SP damping is: ',num2str(DAMP_SP),' and the approximation: ',num2str(DAMP_SP_approx));
Period_sp_msg = strcat('The SP period is: ',num2str(T_SP),' s');
ta_sp_msg = strcat('The SP time to half/double is: ',num2str(ta_SP),' s');
cycles_sp_msg = strcat('The SP Cycles to double or half is: ',num2str(cycles_SP));
decrement_sp_msg = strcat('The SP Logarithmic decrement is: ',num2str(decrement_SP));

disp(SP_msg)
disp(SP_msg2)
disp(poles_sp_msg)
disp(poles_sp_approx_msg)
disp(WN_sp_msg)
disp(DAMP_sp_msg)
disp(Period_sp_msg)
disp(ta_sp_msg)
disp(cycles_sp_msg)
disp(decrement_sp_msg)
disp(SP_msg2)

Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
disp(Warning)
pause


% Calculation of the Natural Frequency Damping ratio and
% time constants for the Phugoid Mode
re_ph=real(phugoid_poles);
im_ph=imag(phugoid_poles);

WN_PH=abs(sqrt(re_ph^2 + im_ph^2));
WN_PH_approx=sqrt(-g*Zu/u1);
DAMP_PH = abs(re_ph/WN_PH);
DAMP_PH_approx= -(Xu+Xtu)/(2*WN_PH_approx);
re_ph_approx = -WN_PH_approx*DAMP_PH_approx;
im_ph_approx = WN_PH_approx*sqrt(1 - DAMP_PH_approx^2);

results.long.PH_wn      = WN_PH;
results.long.PH_damp    = DAMP_PH;

% period
T_PH=(2*pi)/(WN_PH*sqrt(1-DAMP_PH^2));
%time to half or double
ta_PH=0.693/abs(re_ph);
% cycles to double or half
cycles_PH = 0.110*(sqrt(1-DAMP_PH^2))/abs(DAMP_PH);  
% Logarithmic decrement
decrement_PH = 0.693/cycles_PH;

results.long.PH_T   = T_PH;
results.long.PH_T2  = ta_PH;
results.long.cycles_PH    = cycles_PH;
results.long.decrement_PH   = decrement_PH;

% Eigenvalues
PH_msg = strcat('The Phugoid (PH) Analysis');
PH_msg2 = strcat('------------------------------');
poles_ph_msg = strcat('The PH Poles are: real part: ',num2str(re_ph),' & imaginary part: ',num2str(im_ph),'i');
poles_ph_approx_msg = strcat('The PH Poles aproximation are: real part: ',num2str(re_ph_approx),' & imaginary part: ',num2str(im_ph_approx),'i');
WN_ph_msg = strcat('The PH natural frequency is: ',num2str(WN_PH),' and the approximation: ',num2str(WN_PH_approx));
DAMP_ph_msg = strcat('The PH damping is: ',num2str(DAMP_PH),' and the approximation: ',num2str(DAMP_PH_approx));
Period_ph_msg = strcat('The PH period is: ',num2str(T_PH),' s');
ta_ph_msg = strcat('The PH time to half/double is: ',num2str(ta_PH),' s');
cycles_ph_msg = strcat('The PH Cycles to double or half is: ',num2str(cycles_PH));
decrement_ph_msg = strcat('The PH Logarithmic decrement is: ',num2str(decrement_PH));

disp(PH_msg)
disp(PH_msg2)
disp(poles_ph_msg)
disp(poles_ph_approx_msg)
disp(WN_ph_msg)
disp(DAMP_ph_msg)
disp(Period_ph_msg)
disp(ta_ph_msg)
disp(cycles_ph_msg)
disp(decrement_ph_msg)
disp(PH_msg2)

Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
disp(Warning)
pause

% Forwards Speed Stability
crit1 = CTxu-CDu;
if crit1 < 0
    Forward_Speed_Stability =1;
    Forward_Speed_Stability_msg = 'Satisfies with Forward Speed Stability';
else
    Forward_Speed_Stability =0;
    Forward_Speed_Stability_msg = 'Does NOT Satisfy with Forward Speed Stability';
end
% Vertical Speed Stability
crit3 = - CLa;
if crit3 < 0
    Vertical_Speed_Stability =1;
    Vertical_Speed_Stability_msg = 'Satisfies with Vertical Speed Stability';
else
    Vertical_Speed_Stability =0;
    Vertical_Speed_Stability_msg = 'Does NOT Satisfy with Vertical Speed Stability';
end

% Angle of Attack Stability
crit4 = CMa + CMta;
if crit4 < 0
    AoA_Stability =1;
    AoA_Stability_msg = 'Satisfies with Angle of Attack Stability';
else
    AoA_Stability =0;
    AoA_Stability_msg = 'Does NOT Satisfy with Angle of Attack Stability';
end
    
% Pitch Rate Stablity
crit7 = CMq;
if crit7 < 0
    Q_Stability =1;
    Q_Stability_msg = 'Satisfies with Pitch Rate Stability';
else
    Q_Stability =0;
    Q_Stability_msg = 'Does NOT Satisfy with Pitch Rate Stability';
end

% Forward Speed on Pitching Moment
crit7 = CMu;
if crit7 > 0
    u_Q_Stability =1;
    u_Q_Stability_msg = 'Satisfies with Forward Speed on Pitching Moment Stability';
else
    u_Q_Stability =0;
    u_Q_Stability_msg = 'Does NOT Satisfy with Forward Speed on Pitching Moment Stability';
end

% Stability Analysis
Long_criteria_msg = strcat('Longitudinal Criteria for Longitudinal Stability');
long_msg2         = strcat('------------------------------------------------');

disp(Long_criteria_msg)
disp(long_msg2)
disp(Forward_Speed_Stability_msg)
disp(Vertical_Speed_Stability_msg)
disp(AoA_Stability_msg)
disp(Q_Stability_msg)
disp(u_Q_Stability_msg)
disp(long_msg2)

Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
disp(Warning)
pause
    
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

Cy_dr = model.derivadas.Cy_dr;
Cl_dr = model.derivadas.Cl_dr;
Cn_dr = model.derivadas.Cn_dr;
Cy_da = model.derivadas.Cy_da;
Cl_da = model.derivadas.Cl_da;
Cn_da = model.derivadas.Cn_da;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Lateral-Directional Dimensional Stability Derivatives %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1=Ixz/Ixx;
B1=Ixz/Izz;

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

% % Formulations are equivalent
% Ybeta_hat = (Ybeta)/u1;
% Yp_hat = Yp;
% Yr_hat = Yr-u1;
% Yphi_hat = g*cos(theta1);
% Ypsi_hat= 0;
% 
% Lbeta_hat = ((Lbeta + A1*(Nbeta+Ntbeta))/(1-A1*B1))/u1;
% Lp_hat = (Lp + A1*Np)/(1-A1*B1);
% Lr_hat = (Lr + A1*Nr)/(1-A1*B1);
% Lphi_hat = 0;
% Lpsi_hat = 0;
% 
% Nbeta_hat = ((Nbeta + Ntbeta + B1*Lbeta)/(1-A1*B1))/u1;
% Np_hat = (Np + B1*Lp)/(1-A1*B1);
% Nr_hat = (Nr + B1*Lr)/(1-A1*B1);
% Nphi_hat = 0;
% Npsi_hat = 0;
% 
% Phibeta_hat = 0;
% Phip_hat = 1;
% Phir_hat = tan(theta1);
% Phiphi_hat = 0;
% Phipsi_hat= 0;
% 
% Psibeta_hat = 0;
% Psip_hat = 0;
% Psir_hat = 1/cos(theta1);
% Psiphi_hat = 0;
% Psipsi_hat = 0;

Ra_lat      =   [Ybeta/u1,          Yp,         Yr-u1,          g*cos(theta1),     0;
                Lbeta/u1,           Lp,         Lr,             0,                  0;
                (Nbeta+Ntbeta)/u1,  Np,         Nr,             0,                  0;
                0                   1,          tan(theta1),    0,                  0;   
                0,                  0,          sec(theta1),    0,                  0];

Ma_lat      =   [1,         0,          0,          0,      0;
                0,          1,          -Ixz/Ixx,   0,      0;
                0,          -Ixz/Izz,   1,          0,      0;
                0,          0,          0,          1,      0;
                0,          0,          0,          0,      1];

Matrix_lat       =   Ma_lat\Ra_lat;
lat_poles   = eig(Matrix_lat);

% % Formulations are equivalent
% Matrix_lat  =  [Ybeta_hat,          Yp_hat,         Yr_hat,          Yphi_hat,     Ypsi_hat;
%                 Lbeta_hat,          Lp_hat,         Lr_hat,          Lphi_hat,     Lpsi_hat;
%                 Nbeta_hat,          Np_hat,         Nr_hat,          Nphi_hat,     Npsi_hat;
%                 Phibeta_hat,        Phip_hat,       Phir_hat,        Phiphi_hat,   Phipsi_hat;
%                 Psibeta_hat,        Psip_hat,       Psir_hat,        Psiphi_hat,   Psipsi_hat;];
% lat_poles   = eig(Matrix_lat)

Yda = (q1*Sref*Cy_da)/(m);
Ydr = (q1*Sref*Cy_dr)/(m);
Lda = (q1*Sref*b_w*Cl_da)/(Ixx);
Ldr = (q1*Sref*b_w*Cl_dr)/(Ixx);
Nda = (q1*Sref*b_w*Cn_da)/(Izz);
Ndr = (q1*Sref*b_w*Cn_dr)/(Izz);

qql = Yda;
rrl = Ydr;
ssl = (Lda+A1*Nda)/(1-A1*B1);
ttl = (Ldr+A1*Ndr)/(1-A1*B1);
uul = (B1*Lda+Nda)/(1-A1*B1);
vvl = (B1*Ldr+Ndr)/(1-A1*B1);
wwl = 0;
xxl = 0;
yyl = 0;
zzl = 0;

Deflection_lat=[qql rrl; ssl ttl; uul vvl; wwl xxl; yyl zzl];
save Lateral.mat Matrix_lat Deflection_lat

%checks which set of poles are the dutch roll, spiral or roll
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

same_poles = sum(check_ij');

ii=0;
jj=0;
k=0;
DR_poles = [];
Rest_poles = [];
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

if ~isempty(DR_poles) && ~isempty(Rest_poles)

    [B,I] = sort(abs(Rest_poles_vec));

    Spiral_approx = (Lbeta*Nr-Nbeta*Lr)/(Lbeta + Nbeta*A1);
    Rolling_approx = Lp;

    Yaw_pole = Rest_poles(I(1));

    if abs(Spiral_approx) < abs(Rolling_approx)
        Spiral_pole = Rest_poles(I(2));
        Roll_pole = Rest_poles(I(3));
    else
        Spiral_pole = Rest_poles(I(3));
        Roll_pole = Rest_poles(I(2));
    end

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

    % cycles to double or half
    cycles_DR = 0.110*(sqrt(1-ZETA_DR^2))/abs(ZETA_DR);
    % Logarithmic decrement
    decrement_DR = 0.693/cycles_DR;

    results.lat.DR_T    = T_DR;
    results.lat.DR_T2   = ta_DR;

    results.lat.cycles_DR    = cycles_DR;
    results.lat.decrement_DR   = decrement_DR;

    results.lat.ROL_T2  = 0.693/abs(Roll_pole);
    results.lat.ESP_T2  = 0.693/abs(Spiral_pole);

else
    
    results.lat.DR_damp     = [];
    results.lat.DR_wn       = [];
    
    results.lat.DR_T        = [];
    results.lat.DR_T2       = [];

    results.lat.ROL_T2      = [];
    results.lat.ESP_T2      = [];
    
    results.lat.Spiral_pole = [];
    results.lat.Roll_pole   = [];
    results.lat.Yaw_pole    = [];
    
end

    
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

% Aproximations
WN_DR;
WN_DR_approx = sqrt(Nbeta + (1/u1)*(Ybeta*Nr - Nbeta*Yr));
ZETA_DR;
ZETA_DR_approx = - (Nr + Ybeta/u1)/(2*WN_DR_approx);
re_dr_approx = -WN_DR_approx*ZETA_DR_approx;
im_dr_approx = WN_DR_approx*sqrt(1 - ZETA_DR_approx^2);
Spiral_approx = (Lbeta*Nr-Nbeta*Lr)/(Lbeta + Nbeta*A1);
T_spiral_approx = 1/abs(Spiral_approx);
Rolling_approx = Lp;
T_rolling_approx = abs(1/Lp);

% Eigenvalues
DR_msg = strcat('The Dutch Roll (DR) Analysis');
DR_msg2 = strcat('------------------------------');
poles_dr_msg = strcat('The DR Poles are: real part: ',num2str(re_dr),' & imaginary part: ',num2str(im_dr),'i');
poles_dr_approx_msg = strcat('The DR Poles aproximation are: real part: ',num2str(re_dr_approx),' & imaginary part: ',num2str(im_dr_approx),'i');
WN_dr_msg = strcat('The DR natural frequency is: ',num2str(WN_DR),' and the approximation: ',num2str(WN_DR_approx));
DAMP_dr_msg = strcat('The DR damping is: ',num2str(ZETA_DR),' and the approximation: ',num2str(ZETA_DR_approx));
Period_dr_msg = strcat('The DR period is: ',num2str(T_DR),' s');
ta_dr_msg = strcat('The DR time to half/double is: ',num2str(ta_DR),' s');
cycles_dr_msg = strcat('The DR Cycles to double or half is: ',num2str(cycles_DR));
decrement_dr_msg = strcat('The DR Logarithmic decrement is: ',num2str(decrement_DR));
poles_rolling_msg = strcat('The Roll pole is: ',num2str(Roll_pole), 'and the Roll pole approximation is: ',num2str(Rolling_approx));
Period_rolling_msg = strcat('The Rolling period is: ',num2str(T_rolling_approx),' s');
poles_spiral_msg = strcat('The Spiral pole is: ',num2str(Spiral_pole), 'and the Spiral pole approximation is: ',num2str(Spiral_approx));
Period_spiral_msg = strcat('The Spiral period is: ',num2str(T_spiral_approx),' s');

disp(DR_msg)
disp(DR_msg2)
disp(poles_dr_msg)
disp(poles_dr_approx_msg)
disp(WN_dr_msg)
disp(DAMP_dr_msg)
disp(Period_dr_msg)
disp(ta_dr_msg)
disp(cycles_dr_msg)
disp(decrement_dr_msg)
disp(poles_rolling_msg)
disp(Period_rolling_msg)
disp(poles_spiral_msg)
disp(Period_spiral_msg)
disp(DR_msg2)

Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
disp(Warning)
pause

% Side Speed Stability
crit2 = Cyb - CyTb;
if crit2 < 0
    Sideslip_Speed_Stability =1;
    Sideslip_Speed_Stability_msg = 'Satisfies with Side Speed Stability';
else
    Sideslip_Speed_Stability =0;
    Sideslip_Speed_Stability_msg = 'Does NOT Satisfy with Side Speed Stability';
end
    
% Angle of Sideslip Stability 
crit5 = Cnb + CnTb;
if crit5 > 0
    Beta_Stability =1;
    Beta_Stability_msg = 'Satisfies with Angle of Sideslip Stability';
else
    Beta_Stability =0;
    Beta_Stability_msg = 'Does NOT Satisfy with Angle of Sideslip Stability';
end

% Roll Rate Stablity
crit6 = Clp;
if crit6 < 0
    P_Stability =1;
    P_Stability_msg = 'Satisfies with Roll Rate Stability';
else
    P_Stability =0;
    P_Stability_msg = 'Does NOT Satisfy with Roll Rate Stability';
end

% Yaw Rate Stablity
crit8 = Cnr;
if crit8 < 0
    R_Stability =1;
    R_Stability_msg = 'Satisfies with Yaw Rate Stability';
else
    R_Stability =0;
    R_Stability_msg = 'Does NOT Satisfy with Yaw Rate Stability';
end

% Sideslip on Rolling Moment Stablity
crit10 = Clb;
if crit10 < 0
    beta_L_Stability =1;
    beta_L_Stability_msg = 'Satisfies with Sideslip on Rolling Momen Stability';
else
    beta_L_Stability =0;
    beta_L_Stability_msg = 'Does NOT Satisfy with Sideslip on Rolling Momen Stability';
end

% Sideslip on Yawing Momwnt Stablity
crit11 = Cnb;
if crit11 > 0
    beta_N_Stability =1;
    beta_N_Stability_msg = 'Satisfies with Sideslip on Yawing Moment Stability';
else
    beta_N_Stability =0;
    beta_N_Stability_msg = 'Does NOT Satisfy with Sideslip on Yawing Moment Stability';
end

% Spiral Root Stablity
crit11 = (Lbeta*Nr-Nbeta*Lr);
if crit11 > 0
    spiral_Stability =1;
    spiral_Stability_msg = 'Satisfies with Spiral Root Stablity';
else
    spiral_Stability =0;
    spiral_Stability_msg = 'Does NOT Satisfy with Spiral Root Stablity';
end

% Stability Analysis
Lat_criteria_msg = strcat('Lateral-Directiona Criteria for Longitudinal Stability');
lat_msg2         = strcat('-----------------------------------------------------');

disp(Lat_criteria_msg)
disp(lat_msg2)
disp(Sideslip_Speed_Stability_msg)
disp(Beta_Stability_msg)
disp(P_Stability_msg)
disp(R_Stability_msg)
disp(beta_L_Stability_msg)
disp(beta_N_Stability_msg)
disp(spiral_Stability_msg)