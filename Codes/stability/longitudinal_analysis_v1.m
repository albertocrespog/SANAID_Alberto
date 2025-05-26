function Stab_Dyn_Long = logitudinal_analysis_v1(Performance,Stab_Der,conv_UNITS,StabilityModel,Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions)

theta1 = TRIM_RESULTS.trim_alpha;
in2m = conv_UNITS.in2m;
g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

if Variable_Study == 1
    V = conditions.V;
    m_TOW = conditions.m_TOW;
else
    V = Performance.V;
    m_TOW = Weight_tier.m_TOW;
end

% Performance
u1 = V;
rho = Performance.rho;
% Geometric data
S_ref = Geo_tier.S_ref;
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
cmac_w1 = Geo_tier.cmac_w1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%COEFICIENTES Stab_DerLIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ixx = Weight_tier.Ixx;
Iyy = Weight_tier.Iyy;
Izz = Weight_tier.Izz;
Ixz = Weight_tier.Ixz;

cmed = cmac_w1/(2*V);                  %adimensionalización de la cuerda
m1 = 2*m_TOW/(rho*S_w1*V);               %adimensionalización de la masa

CD_alpha = Stab_Der.CD_alpha;
CD = Stab_Der.CD;

CX_alpha = Stab_Der.CXalfa;
CZ_alpha = Stab_Der.CZalfa;
CM_alpha = Stab_Der.CMalfa;

CD_u = Stab_Der.CDu;
CX_u = Stab_Der.CXu;
CL_u = Stab_Der.CLu;
CZ_u = Stab_Der.CZu;
CM_u = Stab_Der.CMu;

CX_q = Stab_Der.CXq;
CL_q = Stab_Der.CLq;
CZ_q = Stab_Der.CZq;
CM_q = Stab_Der.CMq;
CD_q = -CX_q;

CL_alphaDot = Stab_Der.CLalphapunto;
CX_alphaDot = Stab_Der.CXalfapunto;
CZ_alphaDot = Stab_Der.CZalfapunto;
CM_alphaDot = Stab_Der.CMalphapunto;

CX_theta = Stab_Der.CXteta;
CZ_theta = Stab_Der.CZteta;
CM_theta = Stab_Der.CMteta;

CX_delta_e = Stab_Der.CXdeltae;
CD_delta_e = Stab_Der.CDdeltae;
CZ_delta_e = Stab_Der.CZdeltae;
CM_delta_e = Stab_Der.CMdeltae;
CL_delta_e = Stab_Der.CL_delta_e;

CTx1 = Stab_Der.CTx1;
CTx_u = Stab_Der.CTxu;
CTx_alpha = Stab_Der.CTxalpha;
CMT1 = Stab_Der.CMt1;
CMT_u = Stab_Der.CMtu;
CMT_alpha = Stab_Der.CMtalpha;

CL = Stab_Der.CL;
CD = Stab_Der.CD;
CM = Stab_Der.CM;

CL_alpha_ac = Stab_Der.CL_alpha_ac;
CM_alpha_ac = Stab_Der.CM_alpha_ac;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CONSTRUCCIÓN DE LA MATRIZ A y B%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inercias de catia
%inercias adimensionales
Ix1=Ixx/(0.5*rho*V^2*S_w1*b_w1);
Iz1=Izz/(0.5*rho*V^2*S_w1*b_w1);
Ixz1=Ixz/(0.5*rho*V^2*S_w1*b_w1);
Iy1=Iyy/(0.5*rho*V^2*S_w1*cmac_w1);       %adimensionalización de la inercia

%inercias de trabajo
Ix1prima=Ix1/(Ix1*Iz1-Ixz1^2);
Iz1prima=Iz1/(Ix1*Iz1-Ixz1^2);
Ixz1prima=Ixz1/(Ix1*Iz1-Ixz1^2);

q1barw = 0.5*rho*V^2;

Xu      = (-q1barw*S_w1*(CD_u+(2*CD)))/(m_TOW*V);
Xtu     = (q1barw*S_w1*(CTx_u+(2*CTx1)))/(m_TOW*V);
Xa      = (-q1barw*S_w1*(CD_alpha-CL))/m_TOW;
Xq      = (-q1barw*S_w1*CD_q)/(2*m_TOW*V); %% Revisar esta derivada

Zu      = (-q1barw*S_w1*(CL_u+(2*CL)))/(m_TOW*V);
Za      = (-q1barw*S_w1*(CL_alpha_ac + CD))/m_TOW;
ZaDot   = (-q1barw*S_w1*CL_alphaDot*cmac_w1)/(2*m_TOW*V);
Zq      = (-q1barw*S_w1*CL_q*cmac_w1)/(2*m_TOW*V);

Mu      = (q1barw*S_w1*cmac_w1*(CM_u+(2*CM)))/(Iyy*V);
Mtu     = (q1barw*S_w1*cmac_w1*(CMT_u+(2*CMT1)))/(Iyy*V);
Ma      = (q1barw*S_w1*cmac_w1*CM_alpha_ac)/Iyy;
MaDot   = (q1barw*S_w1*(cmac_w1^2)*CM_alphaDot)/(2*Iyy*V);
Mta     = (q1barw*S_w1*cmac_w1*CMT_alpha)/Iyy;
Mq      = (q1barw*S_w1*(cmac_w1^2)*CM_q)/(2*Iyy*V);

Xdeltae = (-q1barw*S_w1*CD_delta_e)/m_TOW;
Zdeltae = (-q1barw*S_w1*CL_delta_e)/m_TOW;
Mdeltae = (q1barw*S_w1*cmac_w1*CM_delta_e)/Iyy;

switch StabilityModel
    case 1
        V = V;
        chi1=CX_alphaDot*cmed/(m1-CZ_alphaDot*cmed);
        chi2=CM_alphaDot*cmed/(m1-CZ_alphaDot*cmed);
        
        %los coeficientes son por orden:
        % a11=((Cxu+Ctxu) + chi1*Czu)/(m1);
        a11=((CX_u) + chi1*CZ_u)/(m1);
        a12=(CX_alpha+chi1*CZ_alpha)/(m1);
        a13=(CX_q*cmed+chi1*(m1+CZ_q*cmed))/(m1);
        a14=(CX_theta+chi1*CZ_theta)/(m1);
        a21=(CZ_u)/(m1-CZ_alphaDot*cmed);
        a22=(CZ_alpha)/(m1-CZ_alphaDot*cmed);
        a23=(m1+CZ_q*cmed)/(m1-CZ_alphaDot*cmed);
        a24=(CZ_theta)/(m1-CZ_alphaDot*cmed);
        %         a31=((CM_u + CMT_u)+chi2*CZ_u)/(Iy1);
        a31=((CM_u)+chi2*CZ_u)/(Iy1);
        %         a32=((CM_alpha + CMT_alpha)+chi2*CZ_alpha)/(Iy1);
        a32=((CM_alpha)+chi2*CZ_alpha)/(Iy1);
        a33=(CM_q*cmed+chi2*(m1+CZ_q*cmed))/(Iy1);
        a34=(chi2*CZ_theta)/(Iy1);
        a41=0;
        a42=0;
        a43=1;
        a44=0;
        
        Matrix_lon=[a11 a12 a13 a14;a21 a22 a23 a24;a31 a32 a33 a34;a41 a42 a43 a44];
        
        %para la matriz B los coeficientes son:
        b1=(CX_delta_e+chi1*CZ_delta_e)/(m1);
        b2=(CZ_delta_e)/(m1-cmed*CZ_alphaDot);
        b3=(CM_delta_e+chi2*CZ_delta_e)/(Iy1);
        b4=0;
        
        Deflection_lon=[b1;b2;b3;b4];
        
        long_poles = eig(Matrix_lon); % Hay que dimensionalizar
        
        
        if Variable_Study == 1
            % no saving for variable study
        else
            save Longitudinal.mat  Matrix_lon Deflection_lon theta1 u1
        end
        
    case 2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Longitudinal Dimensional Stability Derivatives %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        q1barw = 0.5*rho*V^2;
        
        S_w1 = S_w1;
        m_TOW = m_TOW;
        V = V;
        cmac_w1 = cmac_w1;
        CM = CM;
        
        Xu      = (-q1barw*S_w1*(CD_u+(2*CD)))/(m_TOW*V);
        Xtu     = (q1barw*S_w1*(CTx_u+(2*CTx1)))/(m_TOW*V);
        Xa      = (-q1barw*S_w1*(CD_alpha-CL))/m_TOW;
        Xq      = (-q1barw*S_w1*CD_q)/(2*m_TOW*V); %% Revisar esta derivada
        
        Zu      = (-q1barw*S_w1*(CL_u+(2*CL)))/(m_TOW*V);
        Za      = (-q1barw*S_w1*(CL_alpha_ac + CD))/m_TOW;
        ZaDot   = (-q1barw*S_w1*CL_alphaDot*cmac_w1)/(2*m_TOW*V);
        Zq      = (-q1barw*S_w1*CL_q*cmac_w1)/(2*m_TOW*V);
        
        Mu      = (q1barw*S_w1*cmac_w1*(CM_u+(2*CM)))/(Iyy*V);
        Mtu     = (q1barw*S_w1*cmac_w1*(CMT_u+(2*CMT1)))/(Iyy*V);
        Ma      = (q1barw*S_w1*cmac_w1*CM_alpha_ac)/Iyy;
        MaDot   = (q1barw*S_w1*(cmac_w1^2)*CM_alphaDot)/(2*Iyy*V);
        Mta     = (q1barw*S_w1*cmac_w1*CMT_alpha)/Iyy;
        Mq      = (q1barw*S_w1*(cmac_w1^2)*CM_q)/(2*Iyy*V);
        
        Xdeltae = (-q1barw*S_w1*CD_delta_e)/m_TOW;
        Zdeltae = (-q1barw*S_w1*CL_delta_e)/m_TOW;
        Mdeltae = (q1barw*S_w1*cmac_w1*CM_delta_e)/Iyy;
        
        Xu      = (-q1barw*S_w1*(CD_u+(2*CD)))/(m_TOW*V);
        Xtu     = (q1barw*S_w1*(CTx_u+(2*CTx1)))/(m_TOW*V);
        Xa      = (-q1barw*S_w1*(CD_alpha-CL))/m_TOW;
        Xq      = (-q1barw*S_w1*CD_q)/(2*m_TOW*V); %% Revisar esta derivada
        
        Zu      = (-q1barw*S_w1*(CL_u+(2*CL)))/(m_TOW*V);
        Za      = (-q1barw*S_w1*(CL_alpha_ac + CD))/m_TOW;
        ZaDot   = (-q1barw*S_w1*CL_alphaDot*cmac_w1)/(2*m_TOW*V);
        Zq      = (-q1barw*S_w1*CL_q*cmac_w1)/(2*m_TOW*V);
        
        Mu      = (q1barw*S_w1*cmac_w1*(CM_u+(2*CM)))/(Iyy*V);
        Mtu     = (q1barw*S_w1*cmac_w1*(CMT_u+(2*CMT1)))/(Iyy*V);
        Ma      = (q1barw*S_w1*cmac_w1*CM_alpha_ac)/Iyy;
        MaDot   = (q1barw*S_w1*(cmac_w1^2)*CM_alphaDot)/(2*Iyy*V);
        Mta     = (q1barw*S_w1*cmac_w1*CMT_alpha)/Iyy;
        Mq      = (q1barw*S_w1*(cmac_w1^2)*CM_q)/(2*Iyy*V);
        
        Xdeltae = (-q1barw*S_w1*CD_delta_e)/m_TOW;
        Zdeltae = (-q1barw*S_w1*CL_delta_e)/m_TOW;
        Mdeltae = (q1barw*S_w1*cmac_w1*CM_delta_e)/Iyy;
        
        Ra_long     =   [Xu+Xtu,    Xa,         Xq,     -g*cos(theta1);
            Zu,         Za,         Zq+V,  -g*sin(theta1);
            Mu+Mtu,     Ma+Mta      Mq,     0;
            0,          0,          1,      0];
        
        Ma_long     =   [1,         0,          0,      0;
            1,          V-ZaDot,   0,      0;
            0,          -MaDot,     1,      0;
            0,          0,          0,      1];
        
        Matrix_lon      =   Ma_long\Ra_long;
        
        qq = Xdeltae;
        rr = (Zdeltae/(V - ZaDot));
        ss = Mdeltae + ((Zdeltae*MaDot)/(V - ZaDot));
        tt = 0;
        
        Deflection_lon=[qq; rr; ss; tt];
        
        long_poles = eig(Matrix_lon); % Hay que dimensionalizar
        
        if Variable_Study == 1
            % no saving for variable study
        else
            save Longitudinal.mat  Matrix_lon Deflection_lon theta1 u1
        end
        
end

%checks which set of poles are the shorp period and which one the Phugoid
pole1=long_poles(1);
pole1a=long_poles(2);
pole2=long_poles(3);
pole2a=long_poles(4);

Stab_Dyn_Long.long.pole1 = pole1;
Stab_Dyn_Long.long.pole2 = pole1a;
Stab_Dyn_Long.long.pole3 = pole2;
Stab_Dyn_Long.long.pole4 = pole2a;

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

Stab_Dyn_Long.long.SP_poles = SP_poles;
Stab_Dyn_Long.long.PH_poles = phugoid_poles;

%% Short Period
WN_SP = abs(sqrt(re_sp^2 + im_sp^2));
WN_SP_approx=sqrt((Za*Mq/V) - Ma);
DAMP_SP = abs(re_sp/WN_SP);
DAMP_SP_approx= - (Mq + (Za/V) + MaDot)/(2*WN_SP);
re_sp_approx = -WN_SP_approx*DAMP_SP_approx;
im_sp_approx = WN_SP_approx*sqrt(1 - DAMP_SP_approx^2);

% damping ratio
Stab_Dyn_Long.long.SP_wn      = WN_SP;
Stab_Dyn_Long.long.SP_damp    = DAMP_SP;
 
% period
T_SP=(2*pi)/(WN_SP*sqrt(1-DAMP_SP^2));
%time to half or double
ta_SP=0.693/abs(re_sp);
% cycles to double or half
cycles_SP = 0.110*(sqrt(1-DAMP_SP^2))/abs(DAMP_SP);  
% Logarithmic decrement
decrement_SP = 0.693/cycles_SP;

Stab_Dyn_Long.long.SP_T    = T_SP;
Stab_Dyn_Long.long.SP_T2   = ta_SP;
Stab_Dyn_Long.long.cycles_SP    = cycles_SP;
Stab_Dyn_Long.long.decrement_SP   = decrement_SP;

Poles_msg = strcat('The poles of the longitudinal open-loop system are');
Poles_msg2 = strcat('-------------------------------------');

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

% Calculation of the Natural Frequency Damping ratio and
% time constants for the Phugoid Mode
re_ph=real(phugoid_poles);
im_ph=imag(phugoid_poles);

WN_PH=abs(sqrt(re_ph^2 + im_ph^2));
WN_PH_approx=sqrt(-g*Zu/V);
DAMP_PH = abs(re_ph/WN_PH);
DAMP_PH_approx= -(Xu+Xtu)/(2*WN_PH_approx);
re_ph_approx = -WN_PH_approx*DAMP_PH_approx;
im_ph_approx = WN_PH_approx*sqrt(1 - DAMP_PH_approx^2);

Stab_Dyn_Long.long.PH_wn      = WN_PH;
Stab_Dyn_Long.long.PH_damp    = DAMP_PH;

% period
T_PH=(2*pi)/(WN_PH*sqrt(1-DAMP_PH^2));
%time to half or double
ta_PH=0.693/abs(re_ph);
% cycles to double or half
cycles_PH = 0.110*(sqrt(1-DAMP_PH^2))/abs(DAMP_PH);  
% Logarithmic decrement
decrement_PH = 0.693/cycles_PH;

Stab_Dyn_Long.long.PH_T   = T_PH;
Stab_Dyn_Long.long.PH_T2  = ta_PH;
Stab_Dyn_Long.long.cycles_PH    = cycles_PH;
Stab_Dyn_Long.long.decrement_PH   = decrement_PH;

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

% Forwards Speed Stability
crit1 = CTx_u-CD_u;
if crit1 < 0
    Forward_Speed_Stability =1;
    Forward_Speed_Stability_msg = 'Satisfies with Forward Speed Stability';
else
    Forward_Speed_Stability =0;
    Forward_Speed_Stability_msg = 'Does NOT Satisfy with Forward Speed Stability';
end
% Vertical Speed Stability
crit3 = - CL_alpha_ac;
if crit3 < 0
    Vertical_Speed_Stability =1;
    Vertical_Speed_Stability_msg = 'Satisfies with Vertical Speed Stability';
else
    Vertical_Speed_Stability =0;
    Vertical_Speed_Stability_msg = 'Does NOT Satisfy with Vertical Speed Stability';
end

% Angle of Attack Stability
crit4 = CM_alpha_ac + CMT_alpha;
if crit4 < 0
    AoA_Stability =1;
    AoA_Stability_msg = 'Satisfies with Angle of Attack Stability';
else
    AoA_Stability =0;
    AoA_Stability_msg = 'Does NOT Satisfy with Angle of Attack Stability';
end
    
% Pitch Rate Stablity
crit7 = CM_q;
if crit7 < 0
    Q_Stability =1;
    Q_Stability_msg = 'Satisfies with Pitch Rate Stability';
else
    Q_Stability =0;
    Q_Stability_msg = 'Does NOT Satisfy with Pitch Rate Stability';
end

% Forward Speed on Pitching Moment
crit7 = CM_u;
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

Stab_Dyn_Long.Forward_Speed_Stability =Forward_Speed_Stability;
Stab_Dyn_Long.Vertical_Speed_Stability =Vertical_Speed_Stability;
Stab_Dyn_Long.AoA_Stability =AoA_Stability;
Stab_Dyn_Long.Q_Stability =Q_Stability;
Stab_Dyn_Long.u_Q_Stability =u_Q_Stability;

% Does not display message for variable study
if Variable_Study == 1
    % No message displayed
else
    % Shos the poles
    disp(Poles_msg)
    disp(long_poles)

    % Short Period Messages
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
    
    %     Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
    %     disp(Warning)
    
    % Phugoid Messages
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
    
    %     Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
    %     disp(Warning)
    
    % Stability Analysis
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
end

