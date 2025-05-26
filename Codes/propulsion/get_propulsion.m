function [Propulsion,Stab_Der] = get_propulsion(AC_CONFIGURATION,Aero_TH,Aero,Prop_data,Geo_tier, conditions, Performance, conv_UNITS, OUTPUT_read_XLSX)

% Selects the source of the polar model
C_D0 = Aero.Polar.C_D0;
C_D1 = Aero.Polar.C_D1;
C_D2 = Aero.Polar.C_D2;

%% Aerodynamic properties

h = conditions.h;
alpha_f = conditions.alpha_f;
beta_f = conditions.beta_f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
% D_prop = Prop_data.D_prop;
D_prop = OUTPUT_read_XLSX.Propulsive_flags.D_prop;
R_prop = D_prop/2;
S_heli = (pi*R_prop^2);              %superficie de la hélice
S_ref = Geo_tier.S_ref;
V = conditions.V;
rho = Performance.rho;
q_inf = 0.5*rho*V^2;
m_TOW = conditions.m_TOW;
g = conv_UNITS.g;
w_T0 = m_TOW*g;
% Estimation of Desired Thrust
CL = w_T0/(q_inf*S_ref); % assume equilibry steady state flight
CD = C_D0 + C_D1*CL + C_D2*CL^2;
Stab_Der.CD = CD;
D = CD*q_inf*S_ref;

% Propulsive Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%DERIVADA PROPULSIVA LONGITUDINAL%%%%%%%%%%%%%%%%
CD_Total_prel = C_D0 + C_D1*CL + C_D2*CL^2;  %polar aeronave
Drag = q_inf*S_ref*CD_Total_prel;
Fdes = Drag; % assume T = D

%% Propulsion information
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine
% 4: CONSUMO ESPECIFICO: % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% 5: DERIVACION(TURBOFANES) By-pass
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: AVION CIVIL = 1/MILITAR = 2 - Normativa
% 8: CAPACIDAD CALORIFICA DEL COMBUSTIBLE
% 9: DIAMETRO DEL MOTOR(TURBOHELICES)
propul = AC_CONFIGURATION.propulsion;
type_engine = propul(1);
% Type of prop: fixed pitch or variable pitch
type_prop = AC_CONFIGURATION.type_prop; 
bypass_ratio = propul(5); % By-pass Ratio


%% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
switch type_engine
    case 1 % TIPO DE MOTOR --> 1_TURBOFAN 
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX);
        
    case 2 % TIPO DE MOTOR --> 2_TURBOPROP
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX);

    case 3 % TIPO DE MOTOR --> 3_PISTON 
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX);

    case 4 % TIPO DE MOTOR --> 4_ELECTRIC_PROP
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX);

    case 5 % TIPO DE MOTOR --> 5_PISTON_CUSTOM 
        [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha_f,beta_f,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX);
end