function [Propulsion] = get_EngineProperties_v2(V,h,Fdes,n_eng)

% Flight Conditons
[Temp,rho,p,a]=atmos_inter(h);
Mach = V/a;

% Parameter advanced ratio
% Engine = data_engine(k);
load('Data.mat')
load('Data_Prop.mat')
load('Geo_tier.mat')

% Number of Prop used
k = 4;

D_prop=0.6043;        % Propeller Diameter [m].
D= D_prop;
D_propin=D_prop*100/2.54; % Propeller Diameter [in].
A_prop=pi*(D_prop/2)^2;          % Swept area by propeller.

% Propeller data
% Datos genéricos Hélice para 22x12W - They are used to scale the Prop
b_p = (22*2.54/100);
c_p = 3/100;
b_p_c_p = b_p/c_p;
c_prop = D_prop/b_p_c_p;
S_prop = D_prop*c_prop;
AR_prop = (D_prop^2)/S_prop;
RPM_max = 145000/D_propin; % Max RPM by engine builder.

RPM_max = Geo_tier.RPM_max;
% Calculates the Revolutions per second for each engine
Fdes_engine = Fdes/n_eng;
n = Selection_J_CT_F_des(V,CT,k,Geo_tier,Fdes_engine,rho);
check_RPMMAX = n*60; % changes to RPM
if gt(check_RPMMAX,RPM_max)
    msg = 'RPM required greater than RPM max';
    error(msg)
    pause
end

Propulsion.n = n;
Propulsion.RPM_max = RPM_max;

eta_gear=0.96;                   % Gear box efficiency.
eta_m=0.88;                      % Engine efficiency (output/input).
eta_esc=0.98;                    % Speed controller efficiency.
eta_dist=0.96;                   % Shaft efficiency.

% Calculates the Advanced PArameter Ratio
J = V/(n*D);

% Polyfit Coefficients for APC 22x10
% CT_Polyfit{k} = CT{k}(3) + CT{k}(2).*XData_CT{k} + CT{k}(1).*XData_CT{k}.^2;
CT2 = CT{k}(1);
CT1 = CT{k}(2);
CT0 = CT{k}(3);
CT = CT{k}(3) + CT{k}(2)*J + CT{k}(1)*J^2;
Propulsion.CT0 = CT0;
Propulsion.CT1 = CT1;
Propulsion.CT2 = CT2;

% CP_Polyfit{k} = CP{k}(4) + CP{k}(3).*XData_CP{k} + CP{k}(2).*XData_CP{k}.^2 + CP{k}(1).*XData_CP{k}.^3;         
CP3 = CP{k}(1);
CP2 = CP{k}(2);
CP1 = CP{k}(3);
CP0 = CP{k}(4);
CP = CP{k}(4) + CP{k}(3)*J + CP{k}(2)*J^2 + CP{k}(1)*J^3;         
Propulsion.CP0 = CP0;
Propulsion.CP1 = CP1;
Propulsion.CP2 = CP2;
Propulsion.CP3 = CP3;

% CQ_Polyfit{k} = CQ{k}(4) + CQ{k}(3).*XData_CQ{k} + CQ{k}(2).*XData_CQ{k}.^2 + CQ{k}(1).*XData_CQ{k}.^3;
CQ3 = CQ{k}(1);
CQ2 = CQ{k}(2);
CQ1 = CQ{k}(3);
CQ0 = CQ{k}(4);
CQ = CQ{k}(4) + CQ{k}(3)*J + CQ{k}(2)*J^2 + CQ{k}(1)*J^3;
Propulsion.CQ0 = CQ0;
Propulsion.CQ1 = CQ1;
Propulsion.CQ2 = CQ2;
Propulsion.CQ3 = CQ3;

% ethamp_Polyfit{k} = ethamp{k}(4) + ethamp{k}(3).*XData_ethamp{k} + ethamp{k}(2).*XData_ethamp{k}.^2 + ethamp{k}(1).*XData_ethamp{k}.^3;
etha_mp3 = ethamp{k}(1);
etha_mp2 = ethamp{k}(2);
etha_mp1 = ethamp{k}(3);
etha_mp0 = ethamp{k}(4);
ethamp = ethamp{k}(4) + ethamp{k}(3)*J + ethamp{k}(2)*J^2 + ethamp{k}(1)*J^3;
Propulsion.etha_mp0 = etha_mp0;
Propulsion.etha_mp1 = etha_mp1;
Propulsion.etha_mp2 = etha_mp2;
Propulsion.etha_mp3 = etha_mp3;

Ti = CT*rho*(n^2)*(D^4);  
Pi = CP*rho*(n^3)*(D^5); 
Qi = CQ*rho*(n^2)*(D^5);
etha_emp = ethamp*eta_m*eta_gear*eta_esc*eta_dist;
Pe = Pi/etha_emp;

Propulsion.Ti = Ti;
Propulsion.RPM = n*60;
Propulsion.Pi = Pi;
Propulsion.Qi = Qi;
Propulsion.Pe = Pe;
Propulsion.J = J;
Propulsion.CT = CT;
Propulsion.CP = CP;
Propulsion.CQ = CQ;
Propulsion.ethamp = ethamp;
Propulsion.etha_emp = etha_emp;
