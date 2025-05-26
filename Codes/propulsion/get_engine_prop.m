function Propulsion = get_engine_prop(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX) % 1: TIPO DE MOTOR -->  3_PISTON

h_SL = 0;
[Temp_SL,rho_SL,p_SL,a_SL]=atmos_inter_mio(h_SL);
[Temp,rho,p,a]=atmos_inter_mio(h);

% Mach angle at Sea Level and at altitude
M_SL = V/a_SL;
M = V/a;

S_ref = Geo_tier.S_ref;
cmac = Geo_tier.cmac_w1;

qbar = 0.5*rho*V^2;

D_prop = Prop_data.D_prop;
A_prop = Prop_data.A_prop;
R_heli = D_prop/2;
D = Prop_data.D_prop;

propul = AC_CONFIGURATION.propulsion;

% propul(1) = 3; % Type of engine
% propul(2) = n_eng; % Number of engines
% propul(3) = 5; % Thrust (lbf) or Power (shp) per engine
% propul(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% propul(5) = 0; % By-pass
% propul(6) = 0.82; % Prop efficiency
% propul(7) = 1; % Normativa
% propul(8) = 0; % Capacidad calorifica del combustible
% propul(9) = 28*2.54/100; % Diámetro de la hélice

% PROPULSION
tipo_motor = propul(1);
n_eng = propul(2);
P_SL = propul(3);
c_SL = propul(4);
civil_mil = propul(5);
eta_p = propul(6);
derivacion = propul(7);

RPM_min = OUTPUT_read_XLSX.Propulsive_flags.RPM_min; %
RPM_max = OUTPUT_read_XLSX.Propulsive_flags.RPM_max; %
beta_pitch = OUTPUT_read_XLSX.Propulsive_flags.beta_pitch; %
beta_variable_pitch = OUTPUT_read_XLSX.Propulsive_flags.beta_variable_pitch; %

% Calculates the Revolutions per second
T_eng = Fdes/n_eng;
if beta_variable_pitch
    n = get_propeller_models(beta_pitch);
    % n = Selection_J_CT_F_des_v1(V,Prop_data,T_eng,rho);
end


RPM = n*60;

RPM_max = 150000/(D*100/2.54); % Max RPM by engine builder - https://www.apcprop.com/technical-information/rpm-limits/
% Solves for delta_T as a function of Max RPM
RPM_min = 0;
Delta_RPM = RPM_max - RPM_min;
delta_T = n*60/Delta_RPM;

Omega = RPM*(2*pi/60);
% n = Selection_J_CT_F_des(V,Prop_data,T_eng,rho)
check_RPMMAX = n*60; % changes to RPM
% if gt(n*60,RPM_max)
%     msg = 'RPM required greater than RPM max';
%     delta_T
%     disp(msg)
% end

eta_gear = Prop_data.eta_gear;
eta_m = Prop_data.eta_m;
eta_esc = Prop_data.eta_esc;
eta_dist = Prop_data.eta_dist;

N_order_CT = Prop_data.N_order_CT;
N_order_CP = Prop_data.N_order_CP;
N_order_CQ = Prop_data.N_order_CQ;
N_order_etamp = Prop_data.N_order_etamp;

CT_Polyfit = Prop_data.CT_Polyfit;
CP_Polyfit = Prop_data.CP_Polyfit;
CQ_Polyfit = Prop_data.CQ_Polyfit;
etamp_Polyfit = Prop_data.etamp_Polyfit;

% Calculates the Advanced PArameter Ratio
J = V/(n*D);
% Determines CT as a funcion of approximation polynomial
Ct_total =  CT_Polyfit(N_order_CT+1);
for j=1:N_order_CT
    Ct_intermediate = CT_Polyfit(j)*J.^(N_order_CT+1-j);
    Ct_total = Ct_total + Ct_intermediate;
end
% CT = CT0 + CT1*J + CT2*J^2;
CT = Ct_total;

% Determines CP as a funcion of approximation polynomial
Cp_total =  CP_Polyfit(N_order_CP+1);
for j=1:N_order_CP
    Cp_intermediate = CP_Polyfit(j)*J.^(N_order_CP+1-j);
    Cp_total = Cp_total + Cp_intermediate;
end
% CP = CP0 + CP1*J + CP2*J^2 + CP3*J^3;
CP = Cp_total;

% Determines CP as a funcion of approximation polynomial
Cq_total =  CQ_Polyfit(N_order_CQ+1);
for j=1:N_order_CQ
    Cq_intermediate = CQ_Polyfit(j)*J.^(N_order_CQ+1-j);
    Cq_total = Cq_total + Cq_intermediate;
end
% CQ = CQ0 + CQ1*J + CQ2*J^2 + CQ3*J^3;
CQ = Cq_total;

% Determines CP as a funcion of approximation polynomial
etamp_total =  etamp_Polyfit(N_order_etamp+1);
for j=1:N_order_etamp
    etamp_intermediate = etamp_Polyfit(j)*J.^(N_order_etamp+1-j);
    etamp_total = etamp_total + etamp_intermediate;
end
% ethamp = etha_mp0 + etha_mp1*J + etha_mp2*J^2 + etha_mp3*J^3
ethamp = etamp_total;

Ti_eng = CT*rho*(n^2)*(D^4);
Ti = Ti_eng*n_eng;
Pi_eng = CP*rho*(n^3)*(D^5);
Pi = Pi_eng*n_eng;
Qi_eng = CQ*rho*(n^2)*(D^5);
etha_emp = ethamp*eta_m*eta_gear*eta_esc*eta_dist;
Pe = Pi/etha_emp;
Pe_eng = Pe/n_eng;

% Determines d CT/d J
CT_dJ_total =  CT_Polyfit(N_order_CT);
for j=1:(N_order_CT-1)
    CT_dJ_intermediate = ((N_order_CT+1-j))*CT_Polyfit(j)*J.^(N_order_CT-j);
    CT_dJ_total = CT_dJ_total + CT_dJ_intermediate;
end
d_CT_d_J = CT_dJ_total;
d_T_d_J = n_eng*d_CT_d_J*rho*(n^2)*(D^4);
d_Teng_d_J = d_CT_d_J*rho*(n^2)*(D^4);

% Determines d CT/d V
CT_dV_total =  CT_Polyfit(N_order_CT)/(n*D);
for j=1:(N_order_CT-1)
    CT_dV_intermediate = ((N_order_CT+1-j))*CT_Polyfit(j)*J.^(N_order_CT-j);
    CT_dV_total = CT_dV_total + CT_dV_intermediate;
end
d_CT_d_V = CT_dV_total;
d_T_d_V = n_eng*d_CT_d_V*rho*(n^2)*(D^4);
d_Teng_d_V = d_CT_d_V*rho*(n^2)*(D^4);

% Determines d CP/d J
CP_dJ_total =  CP_Polyfit(N_order_CP);
for j=1:(N_order_CP-1)
    CP_dJ_intermediate = ((N_order_CP+1-j))*CP_Polyfit(j)*J.^(N_order_CP-j);
    CP_dJ_total = CP_dJ_total + CP_dJ_intermediate;
end
d_CP_d_J = CP_dJ_total;
d_P_d_J = n_eng*d_CP_d_J*rho*(n^3)*(D^5);
d_Peng_d_J = d_CP_d_J*rho*(n^3)*(D^5);

% Determines d CP/d V
CP_dV_total =  CP_Polyfit(N_order_CP)/(n*D);
for j=1:(N_order_CP-1)
    CP_dV_intermediate = ((N_order_CP+1-j))*CP_Polyfit(j)*J.^(N_order_CP-j);
    CP_dV_total = CP_dV_total + CP_dV_intermediate;
end
d_CP_d_V = CP_dV_total;
d_P_d_V = n_eng*d_CP_d_V*rho*(n^3)*(D^5);
d_Peng_d_V = d_CP_d_V*rho*(n^3)*(D^5);

% Determines deta/dJ
etamp_dJ_total =  etamp_Polyfit(N_order_etamp);
for j=1:(N_order_etamp-1)
    etamp_dJ_intermediate = ((N_order_etamp+1-j))*etamp_Polyfit(j)*J.^(N_order_etamp-j);
    etamp_dJ_total = etamp_dJ_total + etamp_dJ_intermediate;
end
d_etap_d_J = etamp_dJ_total;
% d_etap_d_J = etha_mp1 + 2*etha_mp2*J + 3*etha_mp3*J^2

% Determines deta/dV
etamp_dV_total =  etamp_Polyfit(N_order_etamp)/(n*D);
for j=1:(N_order_etamp-1)
    etamp_dV_intermediate = ((N_order_etamp+1-j))*etamp_Polyfit(j)*J.^(N_order_etamp-j);
    etamp_dV_total = etamp_dV_total + etamp_dV_intermediate;
end
d_etap_d_V = etamp_dV_total;
% d_etha_d_V = 3*V^2*etha_mp3/(n^3*D^3)+2*V*etha_mp2/(n^2*D^2)+etha_mp1/(n*D)
% d_CP_d_V = CP3*V^3/(n^3*D^3) + CP2*V^2/(n^2*D^2) + CP1*V/(n*D) + CP0;

% Estimation of Prop data
A_prop = pi*(D/2)^2;
v_i = -(1/2)*V + sqrt(1/4*V^2 + Ti_eng/(2*rho*A_prop)); % Induced velocity at prop disk
R_inf = (D/2)*sqrt((V + v_i)/(V + 2*v_i));
v_inf = 2*v_i;

%         % Estimation from Nicolai - Chapter 17
%         b = 2; % number of blades
%         a = 5.76; % lift curve slope
%         c = Prop_data.c_prop; % Prop chord at 0.75R
%         R = Prop_data.D_prop/2; % Prop Radii
%         Betar = 19.5*(pi/180); % Angle at tip
%         CD0 = 0.028;
%
%         Beta = @(r) Betar*R./r; % Eqn 17.10 Nicolai Chapter 17
%         vr = @(r) (V/2 +(b*c*a*Omega/(16*pi))) - 1 + ... % Eqn 17.19 Nicolai Chapter 17
%             (1 + (2*Omega.*r.*(Beta(r) - (V./(Omega.*r))))./(((4*pi*V^2)./(b*c*a*Omega)) + V + ((b*c*a*Omega)./(16*pi)))).^(0.5);
%         Ve = @(r) ((vr(r) + V).^2 + (Omega.*r).^2); % Eqn 17.15 Nicolai Chapter 17
%         CL = @(r) a.*(Beta(r)- ((V + vr(r))./(Omega.*r))); % Eqn 17.16 Nicolai Chapter 17
%         dT = @(r)(b/2).*rho.*(Ve(r).^2).*c.*CL(r); % Eqn 17.18 Nicolai Chapter 17
%         dP0 = @(r) ((Omega*r).^2).*Ve(r).*c*CD0; % Eqn 17.23 Nicolai Chapter 17
%         dPi = @(r) Ve(r).*c.*(Beta(r).*Omega.*r - V - vr(r)).*(V + vr(r)); % Eqn 17.23 Nicolai Chapter 17
%
%         % Limits of integration
%         ri = 0.1*R;
%         BR = 0.96*R;
%         T_quad = quad(@(r) dT(r),ri,BR);
%         P_quad = 0.5*rho*b*(quad(@(r) dP0(r),0,R) +  (quad(@(r) dPi(r),ri,BR)));

%% Variables for Prop engine Model
Propulsion.delta_T = delta_T;
% Thrust
Propulsion.Ti = Ti;
Propulsion.Ti_eng = Ti_eng;
% Thrust Derivatives
Propulsion.d_T_d_V = d_T_d_V;
Propulsion.d_Teng_d_V = d_Teng_d_V;
Propulsion.d_CT_d_V = d_CT_d_V;
Propulsion.d_T_d_M = d_T_d_V*a;
Propulsion.d_CT_d_M = d_CT_d_V;
% Power
Propulsion.Pi_eng = Pi_eng;
Propulsion.Pi = Pi;
Propulsion.Pe = Pe;
Propulsion.Pe_eng = Pe_eng;
% power Derivatives
Propulsion.d_P_d_V = d_P_d_V;
Propulsion.d_Peng_d_V = d_Peng_d_V;
Propulsion.d_CP_d_V = d_CP_d_V;
Propulsion.d_P_d_M = d_P_d_V*a;
Propulsion.d_CP_d_M = d_CP_d_V*a;
% Induced Velocity
Propulsion.v_i = v_i; % Induced velocity at prop disk
Propulsion.R_inf = R_inf; % Radius of prop wash at infinity
Propulsion.v_inf = v_inf; % induced velocity at propwash at infinity
Propulsion.delta_T = delta_T;

Propulsion.Ti = Ti;
Propulsion.RPM_max = RPM_max;
Propulsion.RPM = RPM;
Propulsion.RPS = n;
Propulsion.Omega = Omega;
Propulsion.Qi_eng = Qi_eng;

Propulsion.J = J;
Propulsion.CT = CT;
Propulsion.CP = CP;
Propulsion.CQ = CQ;

Propulsion.ethamp = ethamp;
Propulsion.etha_emp = etha_emp;

Propulsion.d_T_d_J = d_T_d_J;
Propulsion.d_CT_d_J = d_CT_d_J;
Propulsion.d_P_d_J = d_P_d_J;
Propulsion.d_Peng_d_J = d_Peng_d_J;
Propulsion.d_CP_d_J = d_CP_d_J;

Propulsion.d_etap_d_V = d_etap_d_V;
Propulsion.d_etap_d_J = d_etap_d_J;

% Dummy for Electric Engine
Propulsion.C = 0;
Propulsion.correccion = 0;