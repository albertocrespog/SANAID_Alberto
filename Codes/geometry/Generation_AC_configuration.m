function [AC_CONFIGURATION] = Generation_AC_configuration(conv_UNITS,OUTPUT_read_XLSX,case_AC)

SF = OUTPUT_read_XLSX.AC_Data_flags.SF;

%% Conditions of the Type of Aircraft and Engine
% Obtain DATA
% Data from the CATIA Original Design
% CASE1 Not making Any changes in the geometry, the airple is really stable in forward flight, trim conditions are really good but:
% -  XCG is at x_XCG: 0.8641 from the nose: 0.2699 m behind the location of the engine thrust line Esto es lo qeu mata en vuelo axial
% -  Aerodynamic center of airplane at 0.9579 m from the nose
% Case2) Removing sweep in wing
% -  XCG is at x_XCG: 0.8997 from the nose: 0.1660 m behind the location of the engine thrust line
% -  Aerodynamic center of airplane at 0.9935 m from the nose
% Case3) Putting wing backwards almos 20 cm and removing sweep
% -  XCG is at x_XCG: 1.0080 from the nose: ONLY 0.0900 m behind the location of the engine thrust line
% -  Aerodynamic center of airplane at 1.1189 m from the nose 
% Case4) Putting wing backwards almos 20 cm and removinG sweep from both wing and Vtail
% -  XCG is at x_XCG: 0.9833 from the nose: ONLY 0.0652 m behind the location of the engine thrust line
% -  Aerodynamic center of airplane at 1.0942 m from the nose 
% CASE = 4;
CASE = OUTPUT_read_XLSX.Geometry_Data_flags.CASE;
%% Determines the ammount of control surface in the w2 aileron & elevator
% For Elevons
% Control_surface 0 = ELEVON;
% Control_surface 1 = W1 AND W2: FLAPS, AILERONS + RUDDERVATOR;
% Control_surface 2 = 2-SURFACE;
% Control_surface 3 = 3-SURFACE;
% Control_surface = 1;
Control_surface = OUTPUT_read_XLSX.Geometry_Data_flags.Control_surface;
%% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% AC_type = 6 - 2 surface: cannard + wing + VTP
% AC_type = 7 - 3 surface: wing + V-tail + VTP
% AC_type = 8 - 2 surface: wing + V-tail + V-tail

AC_type = OUTPUT_read_XLSX.AC_Data_flags.AC_type;

% Available control surfaces
d_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_ail; %Definition of available control surface - aileron
d_ele = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_ele; %Definition of available control surface - elevator
d_elevon = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_elevon; %Definition of available control surface - elevon
d_flap = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_flap; %Definition of available control surface - flap
d_rudder = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudder; %Definition of available control surface - rudder
d_rudvtr = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudvtr; %Definition of available control surface - ruddervator
d_rudvtr2 = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_rudvtr2; %Definition of available control surface - ruddervator
d_can = OUTPUT_read_XLSX.InputGeometry_Data_flags.d_can; %Definition of available control surface - canard

switch AC_type
    case 1 % AC_type = 1 - flying wing + VTP
        W1 =1;
        HTP =0;
        VTP =1;
        twin_VTP = OUTPUT_read_XLSX.AC_Data_flags.twin_VTP;
        Can =0;
        Vee =0;
        Vee2 =0;
        Nac =0;

    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        W1 =1;
        HTP =1;
        VTP =1;
        twin_VTP = OUTPUT_read_XLSX.AC_Data_flags.twin_VTP;
        Can =0;
        Vee =0;
        Vee2 =0;
        Nac =0;
    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        W1 =1;
        HTP =1;
        VTP =1;
        twin_VTP = OUTPUT_read_XLSX.AC_Data_flags.twin_VTP;
        Can =1;
        Vee =0;
        Vee2 =0;
        Nac =0;
    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        W1 =1;
        HTP =0;
        VTP =0;
        twin_VTP = OUTPUT_read_XLSX.AC_Data_flags.twin_VTP;
        Can =0;
        Vee =1;
        Vee2 =0;
        Nac =0;
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        W1 =1;
        HTP =0;
        VTP =0;
        twin_VTP = OUTPUT_read_XLSX.AC_Data_flags.twin_VTP;
        Can =1;
        Vee =1;
        Vee2 =0;
        Nac =0;
    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        W1 =1;
        HTP =0;
        VTP =1;
        twin_VTP = OUTPUT_read_XLSX.AC_Data_flags.twin_VTP;
        Can =1;
        Vee =0;
        Vee2 =0;
        Nac =0;
    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        W1 =1;
        HTP =0;
        VTP =1;
        twin_VTP = OUTPUT_read_XLSX.AC_Data_flags.twin_VTP;
        Can =0;
        Vee =1;
        Vee2 =0;
        Nac =0;
    case 8 % AC_type = 8 - 2 surface: wing + V-tail + V-tail
        W1 =1;
        HTP =0;
        VTP =0;
        twin_VTP = OUTPUT_read_XLSX.AC_Data_flags.twin_VTP;
        Can =0;
        Vee =1;
        Vee2 =1;
        Nac =0;
end

%% Engine location
% Engine_loc = 1 - under wings - n-engines symetrica
% Engine_loc = 2 - fuselage - in front
% Engine_loc = 3 - fuselage - rear
% Engine_loc = 4 - wingtips - n_eng at each side
% Engine_loc = 5 - wingtips - n_eng at each side for 2 wings
Engine_loc = OUTPUT_read_XLSX.Propulsive_flags.Engine_loc;

%% Engine Configuration
% Engine_conf = 1 - pusher prop
% Engine_conf = 2 - puller prop
% Engine_conf = 3 - turbine
% Engine_conf = 2;
Engine_conf = OUTPUT_read_XLSX.Propulsive_flags.Engine_conf;

%% Scaling Factor 
% SF = 1;

%% Geometry
% Initial Geometry
% Propulsion
% Defines the Number of engines according to the type of aircraft defined
switch Engine_loc
    case 1 % Engine_loc = 1 - under wings
        n_eng = 2; %Number of engines
    case 2 % Engine_loc = 2 - fuselage front
        n_eng = 1; %Number of engines
    case 3 % Engine_loc = 2 - fuselage rear
        n_eng = 1; %Number of engines
    case 4 % wingtips - n_eng at each side
        n_eng = 2; %Number of engines
    case 5 % wingtips - n_eng at each side for 2 wings
        n_eng = 4; %Number of engines
end
% CERVERA
if case_AC == 6
    n_eng = 4; %Number of engines
end

%% Propulsion information
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine
% 4: CONSUMO ESPECIFICO: % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% 5: DERIVACION(TURBOFANES) By-pass : 1-baja, 2-media, 3-alta
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: AVION CIVIL = 1/MILITAR = 2 - Normativa
% 8: CAPACIDAD CALORIFICA DEL COMBUSTIBLE
% 9: DIAMETRO DEL MOTOR(TURBOHELICES)

% propul(1) = 4; % Type of engine
% propul(2) = n_eng; % Number of engines
% propul(3) = 5; % Thrust (lbf) or Power (shp) per engine
% propul(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% propul(5) = 2; % By-pass
% propul(6) = 0.82; % Prop efficiency
% propul(7) = 1; % Normativa
% propul(8) = 0; % Capacidad calorifica del combustible
% propul(9) = 28*2.54/100; % Diámetro de la hélice

propul = OUTPUT_read_XLSX.Propulsive_flags.propul;

% Redefines variables
empuje = propul(3);
consumo = propul(4);

% Correct variables according to type of engine
if propul(1) == 1
    propul(3) = empuje * 4.448221615255; %[N] cenversion from pounds force to Newtons
    propul(4) = consumo * 2.832546065 * 10^-5; %[kg/(N*s)] conversion from (lb/lbf*h) to [kg/(N*s)]
elseif propul(1) == 2 || propul(1) == 3
    propul(3) = empuje * 745.699872; %[W] conversion from Hp to Watts
    propul(4) = consumo * 1.68965941 * 10^-7; %[kg/(W*s)] conversion from (lb/shp*h) to [kg/(W*s)]
end
% Stores again the propulsive corrected values
OUTPUT_read_XLSX.Propulsive_flags.propul = propul;

%% Type engine
% type_prop = 1; - propeller driven airplanes with fixed pitch propellers:
% type_prop = 2; - propeller driven airplane with variable pitch (constant speed propellers):
% if propul(1) == 3 || propul(1) == 4
%     % Need to define if fixed pitch or variable pitch propeller
%     type_prop = 1; % for piston and electric assum fixed pitch
% elseif propul(1) == 2 
%     % Need to define if fixed pitch or variable pitch propeller
%     type_prop = 2; % for turboprop assume variable pitch
% elseif propul(1) == 1 
%     type_prop = 0; % for turbofan assume 0
% end
type_prop = OUTPUT_read_XLSX.Propulsive_flags.type_prop;


% Completes the data
%% Propulsive Model
% Number of Prop used
% 1 - APC 20x8
% 2 - APC 22x10 
% 3 - APC 22x12 
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W 
% 6 - APC 21x14  
% Prop_type = 4; 
Prop_type = OUTPUT_read_XLSX.Propulsive_flags.Prop_type;

% AXI 5360/24HD V2 GOLD LINE
% l_eng = SF*0.104; %length
% d_eng = SF*0.063; %diameter
% % Approximate dimenstions of the nacelles
% l_nc = l_eng*2;
% d_nc= d_eng*2; 

l_eng = OUTPUT_read_XLSX.Propulsive_flags.l_eng;
d_eng = OUTPUT_read_XLSX.Propulsive_flags.d_eng;
l_nc = OUTPUT_read_XLSX.Propulsive_flags.l_nc;
d_nc = OUTPUT_read_XLSX.Propulsive_flags.d_nc;

%% Location of the Xcg
% Defines the XCG location
% x_XCG = SF*1.1747; 
% y_XCG = SF*0;
% z_XCG = SF*0;
% 
% % Storing DATA
% XCG_data.x_XCG = x_XCG;
% XCG_data.y_XCG = y_XCG;
% XCG_data.z_XCG = z_XCG;

% Stores Aircraft Configuration
AC_CONFIGURATION.CASE = CASE;
AC_CONFIGURATION.AC_type = AC_type; 
AC_CONFIGURATION.case_AC = case_AC; 
AC_CONFIGURATION.Control_surface = Control_surface;
AC_CONFIGURATION.AC_type = AC_type;
AC_CONFIGURATION.twin_VTP = twin_VTP;
AC_CONFIGURATION.Engine_loc = Engine_loc;
AC_CONFIGURATION.Engine_conf = Engine_conf;
AC_CONFIGURATION.n_eng = n_eng;
AC_CONFIGURATION.Prop_type = Prop_type;
AC_CONFIGURATION.l_nc = l_nc;
AC_CONFIGURATION.d_nc = d_nc;
AC_CONFIGURATION.W1 = W1;
AC_CONFIGURATION.HTP = HTP;
AC_CONFIGURATION.VTP = VTP;
AC_CONFIGURATION.Can = Can;
AC_CONFIGURATION.Vee = Vee;
AC_CONFIGURATION.Vee2 = Vee2;
AC_CONFIGURATION.Nac = Nac;
% Selects the type of control surface
AC_CONFIGURATION.d_ail = d_ail;
AC_CONFIGURATION.d_ele = d_ele;
AC_CONFIGURATION.d_elevon = d_elevon;
AC_CONFIGURATION.d_flap = d_flap;
AC_CONFIGURATION.d_rudder = d_rudder;
AC_CONFIGURATION.d_rudvtr = d_rudvtr;
AC_CONFIGURATION.d_rudvtr2 = d_rudvtr2;
AC_CONFIGURATION.d_can = d_can;
% Propulsion data
AC_CONFIGURATION.propulsion = propul; 
AC_CONFIGURATION.type_prop = type_prop; 