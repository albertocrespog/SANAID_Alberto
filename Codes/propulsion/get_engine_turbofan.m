function Propulsion = get_engine_turbofan(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX) % 1: TIPO DE MOTOR --> 1_TURBOFAN

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
T_SL = propul(3);
c_SL = propul(4);
civil_mil = propul(5);
eta_p = propul(6);
derivacion = propul(7);

% Coefficients correctores consumo específico
a_1=3.559957437510763;
a_2=-10.739698199171459;
a_3= 11.989635150373475;
a_4=-5.869876557884609;
a_5=2.059994459180667;

delta_T = (Fdes/n_eng)/(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
switch derivacion
    case 1
        C =  @(W) c_SL.*correccion.*(1+1.2*M).*sqrt(Temp/Temp_SL);
    case 2
        C =  @(W) c_SL.*correccion.*(1+0.33*M).*sqrt(Temp/Temp_SL);
    case 3
        C =  @(W) c_SL.*correccion.*(1+0.16875*M).*sqrt(Temp/Temp_SL);
end
Ti_eng = delta_T.*(T_SL.*(1+0.2.*M.^2).^(1.4./0.4).*(1-0.49.*sqrt(M)).*(rho./rho_SL));
Ti = Ti_eng*n_eng;
CT = Ti_eng / (qbar*S_ref);
% Propulsive derivatives as a funcion of Speed
dT_dV = 1.400000000*delta_T*T_SL*(1+.2*V^2/a^2)^2.500000000*(1-.49*sqrt(V/a))*rho*V/(rho_SL*a^2)-...
    .2450000000*delta_T*T_SL*(1+.2*V^2/a^2)^3.500000000*rho/(sqrt(V/a)*a*rho_SL);
dCT_dV = 2.800000000*delta_T*T_SL*(1+.2*V^2/a^2)^2.500000000*(1-.49*sqrt(V/a))/(rho_SL*V*S_ref*a^2)-...
    .4900000000*delta_T*T_SL*(1+.2*V^2/a^2)^3.500000000/(sqrt(V/a)*a*rho_SL*V^2*S_ref)-...
    4.000000000*delta_T*T_SL*(1+.2*V^2/a^2)^3.500000000*(1-.49*sqrt(V/a))/(rho_SL*V^3*S_ref);
% Propulsive derivatives as a funcion of Mach
dT_dM =  1.400000000*delta_T*T_SL*(1+.2*M^2)^2.500000000*(1-.49*sqrt(M))*rho*M/rho_SL-...
    .2450000000*delta_T*T_SL*(1+.2*M^2)^3.500000000*rho/(sqrt(M)*rho_SL);
dCT_dM = 2.800000000*delta_T*T_SL*(1+.2*M^2)^2.500000000*(1-.49*sqrt(M))/(rho_SL*M*a^2*S_ref)-...
    .4900000000*delta_T*T_SL*(1+.2*M^2)^3.500000000/(M^(5/2)*rho_SL*a^2*S_ref)-...
    4.000000000*delta_T*T_SL*(1+.2*M^2)^3.500000000*(1-.49*sqrt(M))/(rho_SL*M^3*a^2*S_ref);

Propulsion.CT = CT;

Propulsion.delta_T = delta_T;
Propulsion.C = C;
Propulsion.correccion = correccion;
Propulsion.Ti = Ti;
Propulsion.Ti_eng = Ti_eng;
Propulsion.d_T_d_V = dT_dV;
Propulsion.d_CT_d_V = dCT_dV;
Propulsion.d_T_d_M = dT_dM;
Propulsion.d_CT_d_M = dCT_dM;
%% Dummy variables for Jet Engine Model engine model
% Power
Propulsion.Pi_eng = 0;
Propulsion.Pi = 0;
Propulsion.Pe = 0;
Propulsion.Pe_eng = 0;
% power Derivatives
Propulsion.d_P_d_V = 0;
Propulsion.d_CP_d_V = 0;
Propulsion.d_P_d_M = 0;
Propulsion.d_CP_d_M = 0;
% Estimation of Prop data
A_prop = 0;
v_i = 0; % Induced velocity at prop disk
R_inf = 0;
v_inf = 2*v_i;
% Storing data
Propulsion.v_i = v_i; % Induced velocity at prop disk
Propulsion.R_inf = R_inf; % Radius of prop wash at infinity
Propulsion.v_inf = v_inf; % induced velocity at propwash at infinity        Propulsion.delta_T = delta_T;

% Dummy variables for propulsion
Propulsion.RPM_max = 0;
Propulsion.RPM = 0;
Propulsion.RPS = 0;
Propulsion.Omega = 0;
Propulsion.Qi_eng = 0;

Propulsion.J = 0;
Propulsion.CP = 0;
Propulsion.CQ = 0;

Propulsion.ethamp = 0;
Propulsion.etha_emp = 0;

Propulsion.d_T_d_J = 0;
Propulsion.d_CT_d_J = 0;
Propulsion.d_P_d_J = 0;
Propulsion.d_Peng_d_J = 0;
Propulsion.d_CP_d_J = 0;

Propulsion.d_etap_d_V = 0;
Propulsion.d_etap_d_J = 0;
