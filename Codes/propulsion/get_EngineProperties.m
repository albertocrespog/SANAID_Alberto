function [Propulsion fig] = get_EngineProperties(h,V,rho,alpha,Geo_tier,Fdes,Prop_data,n_eng)

% rho = Data_Trim.rho;
% Parameter advanced ratio
% Engine = data_engine(k);
% load('Data.mat')
% load('Data_Prop.mat')

D = Prop_data.D_prop;
RPM_max = Prop_data.RPM_max;
% Calculates the Revolutions per second
T_eng = Fdes/n_eng;
n = Selection_J_CT_F_des(V,Prop_data,T_eng,rho);
check_RPMMAX = n*60; % changes to RPM
if gt(n,check_RPMMAX)
    msg = 'RPM required greater than RPM max';
    error(msg)
    pause
end

Propulsion.n = n;
Propulsion.RPM_max = RPM_max;

eta_gear = Prop_data.eta_gear;
eta_m = Prop_data.eta_m;
eta_esc = Prop_data.eta_esc;
eta_dist = Prop_data.eta_dist;

% Thrust coefficients from wind tunnel
% Polyfit Coefficients for APC 22x10
% CT_Polyfit{k} = CT{k}(3) + CT{k}(2).*XData_CT{k} + CT{k}(1).*XData_CT{k}.^2;
CT0 = Prop_data.CT0;
CT1 = Prop_data.CT1;
CT2 = Prop_data.CT2;
% CP_Polyfit{k} = CP{k}(4) + CP{k}(3).*XData_CP{k} + CP{k}(2).*XData_CP{k}.^2 + CP{k}(1).*XData_CP{k}.^3;         
CP0 = Prop_data.CP0;
CP1 = Prop_data.CP1;
CP2 = Prop_data.CP2;
CP3 = Prop_data.CP3;
% CQ_Polyfit{k} = CQ{k}(4) + CQ{k}(3).*XData_CQ{k} + CQ{k}(2).*XData_CQ{k}.^2 + CQ{k}(1).*XData_CQ{k}.^3;
CQ0 = Prop_data.CQ0;
CQ1 = Prop_data.CQ1;
CQ2 = Prop_data.CQ2;
CQ3 = Prop_data.CQ3;
% ethamp_Polyfit{k} = ethamp{k}(4) + ethamp{k}(3).*XData_ethamp{k} + ethamp{k}(2).*XData_ethamp{k}.^2 + ethamp{k}(1).*XData_ethamp{k}.^3;
etha_mp0 = Prop_data.etha_mp0;
etha_mp1 = Prop_data.etha_mp1;
etha_mp2 = Prop_data.etha_mp2;
etha_mp3 = Prop_data.etha_mp3;

% Calculates the Advanced PArameter Ratio
J = V/(n*D);

CT = CT0 + CT1*J + CT2*J^2;
CP = CP0 + CP1*J + CP2*J^2 + CP3*J^3;
CQ = CQ0 + CQ1*J + CQ2*J^2 + CQ3*J^3;
ethamp = etha_mp0 + etha_mp1*J + etha_mp2*J^2 + etha_mp3*J^3;

Ti = CT*rho*(n^2)*(D^4);  
Pi = CP*rho*(n^3)*(D^5); 
Qi = CQ*rho*(n^2)*(D^5);
etha_emp = ethamp*eta_m*eta_gear*eta_esc*eta_dist;
Pe = Pi/etha_emp;

% For fixed pitch propellers
d_etap_d_J = etha_mp1 + 2*etha_mp2*J + 3*etha_mp3*J^2;% 

Propulsion.Ti = Ti;
Propulsion.Pi = Pi;
Propulsion.Qi = Qi;
Propulsion.Pe = Pe;
Propulsion.J = J;
Propulsion.CT = CT;
Propulsion.CP = CP;
Propulsion.CQ = CQ;
Propulsion.ethamp = ethamp;

Propulsion.CT0 = CT0;
Propulsion.CT1 = CT1;
Propulsion.CT2 = CT2;

Propulsion.CP0 = CP0;
Propulsion.CP1 = CP1;
Propulsion.CP2 = CP2;
Propulsion.CP3 = CP3;

Propulsion.CQ0 = CQ0;
Propulsion.CQ1 = CQ1;
Propulsion.CQ2 = CQ2;
Propulsion.CQ3 = CQ3;

Propulsion.etha_mp0 = etha_mp0;
Propulsion.etha_mp1 = etha_mp1;
Propulsion.etha_mp2 = etha_mp2;
Propulsion.etha_mp3 = etha_mp3;

Propulsion.d_etap_d_J = d_etap_d_J;