function [Propulsion] = get_EngineProperties_v4(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX)

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

% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
switch tipo_motor
    case 1 % 1: TIPO DE MOTOR --> 1_TURBOFAN
        Propulsion = get_engine_turbofan(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX); % 1: TIPO DE MOTOR --> 1_TURBOFAN
    case 2 % 1: TIPO DE MOTOR -->  2_TURBOPROP 
        Propulsion = get_engine_turboprop(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX); % 1: TIPO DE MOTOR -->  2_TURBOPROP 
    case 3 % 1: TIPO DE MOTOR -->3_PISTON 
        Propulsion = get_engine_piston(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX); % 1: TIPO DE MOTOR -->  3_PISTON 
    case 4 % 1: TIPO DE MOTOR --> 4_ELECTRIC_PROP
        Propulsion = get_engine_electric(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX); % 1: TIPO DE MOTOR -->  3_PISTON
    case 5 % 1: TIPO DE MOTOR --> 5_PROPULSION_CUSTOM
        Propulsion = get_engine_piston_custom(h,V,rho,alpha,beta,Geo_tier,Fdes,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX); % 1: TIPO DE MOTOR -->  5_PISTON
end
