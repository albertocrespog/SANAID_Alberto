function Propulsion = get_engine_piston(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX) % 1: TIPO DE MOTOR -->  3_PISTON

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

%         % TARSIS Propulsive Model
%         Data_Engine = 0;
%         % Loads Data from engine model previously obtained
%         if Data_Engine == 0
%             %     PROP_var = Propulsion_Sim_VAR_v3(m_TO_1,V,conv_UNITS,S_w,h,C_D0,C_D1,C_D2,rho);
%             %     if prop_NEW==1;
%             PROP_var = Propulsion_Sim_VAR_Oct_15_v1(m_TO_1,V,conv_UNITS,S_w,h,C_D0,C_D1,C_D2,rho,prop_NEW);
%             %     end
%             P_poly = PROP_var.P_poly;
%             c_poly = PROP_var.c_poly;
%             delta_p_calc = PROP_var.delta_p_calc;
%
%             NP.PROP_var = PROP_var;
%         else
%             Data_P_poly_v3.mat
%         end
%
%         % Flag to define throtle position
%         % Posicion_Palanca = 1 - posición en calculada
%         % Posicion_Palanca = 0 - empuje 0
%         % Posicion_Palanca = 2 - posición maxima
%
%         if PROPULSION == 1
%
%             if Posicion_Palanca == 0
%                 deltaP = 0;
%                 T = 0;
%                 NP.T = T;
%                 NP.deltaP = deltaP;
%             elseif Posicion_Palanca == 1
%                 deltaP = delta_p_calc;
%
%                 if prop_NEW==1
%                     [T,c,Vnmax,delta_act] = propulsive_model_v5(V,h,deltaP);
%                 else
%                     [T,c,Vnmax] = propulsive_model_v3(V,h,deltaP);
%                 end
%
%                 NP.T = T;
%                 NP.deltaP = deltaP;
%             elseif Posicion_Palanca == 2
%                 deltaP = 1.00;
%
%                 if prop_NEW==1
%                     [T,c,Vnmax,delta_act] = propulsive_model_v5(V,h,deltaP);
%                 else
%                     [T,c,Vnmax] = propulsive_model_v3(V,h,deltaP);
%                 end
%
%                 NP.T = T;
%                 NP.deltaP = deltaP;
%             end
%
%         else
%             CL = g*m_TO_1/(q_inf*S_w);
%             C_D  = C_D0 + C_D1*CL + C_D2*CL^2;
%             T   = q_inf*S_w*C_D;
%             NP.T = T;
%         end

b_prop = OUTPUT_read_XLSX.Propulsive_flags.b_prop;
delta_T = (Fdes/n_eng)/(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
C = c_SL.*correccion.*V./eta_p;
Ti_eng = delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
Ti = Ti_eng*n_eng;
Pi_eng = Ti_eng*V;
Pi = Pi_eng*n_eng;
Pe = Pi/eta_p;
Pe_eng = Pe/n_eng;
CT = Ti_eng / (qbar*S_ref);
CP = Pi_eng / (qbar*S_ref*cmac);

% Propulsive derivatives as a funcion of Speed
dT_dV = -.1324503311*delta_T*T_SL*eta_p*(8.55*rho/rho_SL-1)/V^2;
dCT_dV = -.7947019866*delta_T*T_SL*eta_p*(8.55*rho/rho_SL-1)/(V^4*rho*S_ref);
dP_dV = 0;
dCP_dV = -.5298013244*delta_T*T_SL*(8.55*rho/rho_SL-1)/(V^3*rho*S_ref*cmac);

% Propulsive derivatives as a funcion of Mach
dT_dM = -.1324503311*delta_T*T_SL*eta_p*(8.55*rho/rho_SL-1)/(M^2*a);
dCT_dM = -.7947019866*delta_T*T_SL*eta_p*(8.55*rho/rho_SL-1)/(M^4*a^3*rho*S_ref);
dP_dM = 0;
dCP_dM = -.5298013244*delta_T*T_SL*(8.55*rho/rho_SL-1)/(M^3*a^2*rho*S_ref*cmac);

% Estimation of Prop data
A_prop = pi*(D/2)^2;
v_i = -(1/2)*V + sqrt(1/4*V^2 + Ti_eng/(2*rho*A_prop)); % Induced velocity at prop disk
R_inf = (D/2)*sqrt((V + v_i)/(V + 2*v_i));
v_inf = 2*v_i;

%% Variables for Prop engine Model
Propulsion.delta_T = delta_T;
Propulsion.C = C;
Propulsion.correccion = correccion;
% Thrust
Propulsion.Ti = Ti;
Propulsion.Ti_eng = Ti_eng;
% Thrust Derivatives
Propulsion.d_T_d_V = dT_dV;
Propulsion.d_CT_d_V = dCT_dV;
Propulsion.d_T_d_M = dT_dM;
Propulsion.d_CT_d_M = dCT_dM;
% Power
Propulsion.Pi_eng = Pi_eng;
Propulsion.Pi = Pi;
Propulsion.Pe = Pe;
Propulsion.Pe_eng = Pe_eng;
% power Derivatives
Propulsion.d_P_d_V = dP_dV;
Propulsion.d_CP_d_V = dCP_dV;
Propulsion.d_P_d_M = dP_dM;
Propulsion.d_CP_d_M = dCP_dM;
% Induced Velocity
Propulsion.v_i = v_i; % Induced velocity at prop disk
Propulsion.R_inf = R_inf; % Radius of prop wash at infinity
Propulsion.v_inf = v_inf; % induced velocity at propwash at infinity
Propulsion.delta_T = delta_T;


Propulsion_electric = get_engine_electric(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION);


% Dummy variables for propulsion
Propulsion.RPM_max = Propulsion_electric.RPM_max;
Propulsion.RPM = Propulsion_electric.RPS;
Propulsion.RPS = Propulsion_electric.RPS;
Propulsion.Omega = Propulsion_electric.Omega;
Propulsion.Qi_eng = Propulsion_electric.Qi_eng;

Propulsion.J = Propulsion_electric.J;
Propulsion.CP = Propulsion_electric.CP;
Propulsion.CQ = Propulsion_electric.CQ;

Propulsion.ethamp = Propulsion_electric.ethamp;
Propulsion.etha_emp = Propulsion_electric.etha_emp;

Propulsion.d_T_d_J = Propulsion_electric.d_T_d_J;
Propulsion.d_CT_d_J = Propulsion_electric.d_CT_d_J;
Propulsion.d_P_d_J = Propulsion_electric.d_P_d_J;
Propulsion.d_Peng_d_J = Propulsion_electric.d_Peng_d_J;
Propulsion.d_CP_d_J = Propulsion_electric.d_CP_d_J;

Propulsion.d_etap_d_V = Propulsion_electric.d_etap_d_V;
Propulsion.d_etap_d_J = Propulsion_electric.d_etap_d_J;
