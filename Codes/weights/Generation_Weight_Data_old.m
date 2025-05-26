function Weight_tier = Generation_Weight_Data_old(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,Weight_Estimation,SF,f_f)
% Defines Estimation of Weights according to Cefiro III densities

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
case_AC = AC_CONFIGURATION.case_AC;

g = conv_UNITS.g;

% Number of engines
n_eng = AC_CONFIGURATION.n_eng;

% Fuselage geometry
Surf_TOT = Body_Geo.Surf_TOT;
Vol_TOT = Body_Geo.Vol_TOT;
z_Area_b_1_4 = Body_Geo.z_Area_b_1_4;
z_Area_b_3_4 = Body_Geo.z_Area_b_3_4;
x_Area_b_1_4 = Body_Geo.x_Area_b_1_4;
x_Area_b_3_4 = Body_Geo.x_Area_b_3_4;
x_Area_b_max = Body_Geo.x_Area_b_max;
Area_b_max = Body_Geo.Area_b_max;
w_Area_b_max = Body_Geo.w_Area_b_max;
h_Area_b_max = Body_Geo.h_Area_b_max;
Area_top = Body_Geo.Area_top;
Area_side = Body_Geo.Area_side;
length_fus = Body_Geo.l_fus;
Area_body = Body_Geo.Area_body;
x_Area_body = Body_Geo.x_Area_body;

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        % Weight of propulsion systems
        W_eng = 0.5; % 
        m_engs = n_eng*W_eng;
        m_prop = 0.24;
        m_props = n_eng*m_prop;
        m_ESC = 0.145;
        m_ESCs = n_eng*m_ESC;
        m_prop = m_engs + m_props + m_ESCs;
        
        % subsystems
        n_servos = 10;
        m_servo = 0.05;
        m_servos = m_servo*n_servos;
        m_wiring = 0.200;
        m_metal = 0.250;
        m_subsystems = m_servos + m_wiring + m_metal;
        
        m_systems = 1;
        m_landing_gear = 0;
        
        % Weight of PAX
        n_pax = 0;
        n_crew = 0;
        m_pax = 85; % mass for passanger
        m_crew = 85; % mass for crew
        m_luggage = 20; % mass for luggage and passanger
        m_pax_crew = (n_crew + n_pax)*m_luggage + n_crew*m_crew + n_pax*m_pax;
        
        % Weight of Payload
        m_cargo = 1;
        m_payload = m_pax_crew + m_cargo;
        
        % Energy (fuel or batteries)
        m_batteries = 7;
        m_fuel = 0;
        m_energy = m_batteries + m_fuel;
        
        flag_landing_gear = 0; % determines if landing gear present (1 yes 0 no)
    case 2 % case_AC = 2 - EMERGENTIA 1:2
                % Weight of propulsion systems
        W_eng = 0.5; % 
        m_engs = n_eng*W_eng;
        m_prop = 0.24;
        m_props = n_eng*m_prop;
        m_ESC = 0.145;
        m_ESCs = n_eng*m_ESC;
        m_prop = m_engs + m_props + m_ESCs;
        
        % subsystems
        n_servos = 10;
        m_servo = 0.05;
        m_servos = m_servo*n_servos;
        m_wiring = 0.200;
        m_metal = 0.250;
        m_subsystems = m_servos + m_wiring + m_metal;
        
        m_systems = 1;
        m_landing_gear = 0;
        
        % Weight of PAX
        n_pax = 0;
        n_crew = 0;
        m_pax = 85; % mass for passanger
        m_crew = 85; % mass for crew
        m_luggage = 20; % mass for luggage and passanger
        m_pax_crew = (n_crew + n_pax)*m_luggage + n_crew*m_crew + n_pax*m_pax;
        
        % Weight of Payload
        m_cargo = 0;
        m_payload = m_pax_crew + m_cargo;
        
        % Energy (fuel or batteries)
        m_batteries = 2;
        m_fuel = 0;
        m_energy = m_batteries + m_fuel; 
        
        flag_landing_gear = 0;
    case 3 % case_AC = 3 - PEPIÑO XXL
        % Weight of propulsion systems
        W_eng = 3; % PEPIÑO XXL - DL170
        m_engs = n_eng*W_eng;
        m_prop = 0.24;
        m_props = n_eng*m_prop;
        m_ESC = 0.145;
        m_ESCs = n_eng*m_ESC;
        m_prop = m_engs + m_props + m_ESCs;
        
        % subsystems
        n_servos = 10;
        m_servo = 0.150;
        m_servos = m_servo*n_servos;
        m_wiring = 0.400;
        m_metal = 0.400;
        m_subsystems = m_servos + m_wiring + m_metal;
        
        m_systems = 10;
        m_landing_gear = 0;
        
        % Weight of PAX
        n_pax = 0;
        n_crew = 0;
        m_pax = 85; % mass for passanger
        m_crew = 85; % mass for crew
        m_luggage = 20; % mass for luggage and passanger
        m_pax_crew = (n_crew + n_pax)*m_luggage + n_crew*m_crew + n_pax*m_pax;
        
        % Weight of Payload
        m_cargo = 10;
        m_payload = m_pax_crew + m_cargo;
        
        % Energy (fuel or batteries)
        m_batteries = 1;
        m_fuel = 15;
        m_energy = m_batteries + m_fuel;
        
        flag_landing_gear = 0;
                
    case 4 % case_AC = 4 - Comercial 

        % Weight of propulsion systems
        W_eng = 61; %
        m_engs = n_eng*W_eng;
        m_prop = 5;
        m_props = n_eng*m_prop;
        m_ESC = 0;
        m_ESCs = n_eng*m_ESC;
        m_prop = m_engs + m_props + m_ESCs;

        % Subsystems
        n_servos = 0;
        m_servo = 0;
        m_servos = m_servo*n_servos;
        m_wiring = 0;
        m_metal = 0;
        m_subsystems = m_servos + m_wiring + m_metal;
        
        % Weight of PAX
        n_pax = 6;
        n_crew = 1;
        m_pax = 85; % mass for passanger
        m_crew = 85; % mass for crew
        m_luggage = 20; % mass for luggage and passanger
        m_pax_crew = (n_crew + n_pax)*m_luggage + n_crew*m_crew + n_pax*m_pax;
        
        % Weight of Payload
        m_cargo = 0;
        m_payload = m_pax_crew + m_cargo;
        
        % Energy (fuel or batteries)
        m_batteries = 0;
        m_fuel = 360;
        m_energy = m_batteries + m_fuel;
        
        flag_landing_gear = 1;
        
    case 5 % case_AC = 5 - WIG 
        % Weight of propulsion systems
        W_eng = 300; % 
        m_engs = n_eng*W_eng;
        m_prop = 20;
        m_props = n_eng*m_prop;
        m_ESC = 0;
        m_ESCs = n_eng*m_ESC;
        m_prop = m_engs + m_props + m_ESCs;
        
        % subsystems
        n_servos = 0;
        m_servo = 0;
        m_servos = m_servo*n_servos;
        m_wiring = 0;
        m_metal = 0;
        m_subsystems = m_servos + m_wiring + m_metal;
        
        % Weight of PAX
        n_pax = 0;
        n_crew = 2;
        m_pax = 85; % mass for passanger
        m_crew = 85; % mass for crew
        m_luggage = 20; % mass for luggage and passanger
        m_pax_crew = (n_crew + n_pax)*m_luggage + n_crew*m_crew + n_pax*m_pax;
        
        % Weight of Payload
        m_cargo = 20*1000;
        m_payload = m_pax_crew + m_cargo;
        
        % Energy (fuel or batteries)
        m_batteries = 0;
        m_fuel = 12700;
        m_energy = m_batteries + m_fuel;
        
        flag_landing_gear = 0;
    case 6 % case_AC = 1 - CERVERA
        % Weight of propulsion systems
        W_eng = 0; % 
        m_engs = n_eng*W_eng;
        m_prop = 0;
        m_props = n_eng*m_prop;
        m_ESC = 0;
        m_ESCs = n_eng*m_ESC;
        m_prop = m_engs + m_props + m_ESCs;
        
        % subsystems
        n_servos = 0;
        m_servo = 0.05;
        m_servos = m_servo*n_servos;
        m_wiring = 0;
        m_metal = 0;
        m_subsystems = m_servos + m_wiring + m_metal;
        
        m_systems = 0;
        m_landing_gear = 0;
        
        % Weight of PAX
        n_pax = 0;
        n_crew = 0;
        m_pax = 85; % mass for passanger
        m_crew = 85; % mass for crew
        m_luggage = 20; % mass for luggage and passanger
        m_pax_crew = (n_crew + n_pax)*m_luggage + n_crew*m_crew + n_pax*m_pax;
        
        % Weight of Payload
        m_cargo = 12;
        m_payload = m_pax_crew + m_cargo;
        
        % Energy (fuel or batteries)
        m_batteries = 5;
        m_fuel = 0;
        m_energy = m_batteries + m_fuel;
        
        flag_landing_gear = 0; % determines if landing gear present (1 yes 0 no)
end

switch Weight_Estimation
    case 1 % Composites Cefiro III
        
        [rho_f,rho_fairing,rho_w,rho_HTP,rho_VTP,rho_tb]=CFIIIdensities; % calculated from Céfiro III
        rho_Vee = (rho_HTP + rho_VTP)/2;
        m_fus = rho_f*Surf_TOT;               % Fuselage weight
        m_fairing = rho_fairing*Surf_TOT;               % Fairing weight
        rho_fus_fairing = (rho_f + rho_fairing)/2;
        rho_landing_gear = 0.04; % factor lineal landing gear
        rho_engine = 1; % factor lineal engine installation
        rho_misc = 0.10; % factor lineal msc
         
        if W1 == 1
            S_w1_s = Geo_tier.S_w1_s;
            m_w1 = f_f*rho_w*S_w1_s;                             % Wing weight.
        else
            m_w1 = 0;
        end
        if HTP == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_HTP = f_f*rho_HTP*S_w2_s;                    %Peso HTP
        else
            m_HTP = 0;
        end
        if VTP == 1
            S_VTP_s = Geo_tier.S_VTP_s;
            m_VTP = f_f*rho_VTP*S_VTP_s;                    %Peso HTP
        else
            m_VTP = 0;
        end
        if Can == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_Can = f_f*rho_w*S_w2_s;                  %Peso canard
        else
            m_Can = 0;
        end
        if Vee == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_Vee = f_f*rho_Vee*S_w2_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT;                      %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_fus_fairing;
                
        % Estimation
        m_prop = m_prop*rho_engine; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        MTOW_estimation = m_empty_estimation + m_payload + m_energy;
        m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear;                        %Peso tren de aterrizaje
        
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;
        if case_AC == 6
            m_empty = 7.9;
        end
        m_TOW = m_empty + m_payload + m_energy;

    case 2 % Wood Cefiro I
        [rho_f,rho_nose,rho_w,rho_HTP,rho_VTP,rho_tb]=CFIdensities; % calculated from Céfiro I
        rho_Vee = (rho_HTP + rho_VTP)/2;
        m_fus = rho_f*Surf_TOT;               % Fuselage weight
        m_fairing = rho_nose*Surf_TOT;               % Fairing weight
        rho_fus_fairing = (rho_f + rho_nose)/2;
        rho_landing_gear = 0.05; % factor lineal landing gear
        rho_engine = 1; % factor lineal engine installation
        rho_misc = 0.10; % factor lineal msc
    
        if W1 == 1
            S_w1_s = Geo_tier.S_w1_s;
            m_w1 = f_f*rho_w*S_w1_s;                             % Wing weight.
        else
            m_w1 = 0;
        end
        if HTP == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_HTP = f_f*rho_HTP*S_w2_s;                    %Peso HTP
        else
            m_HTP = 0;
        end
        if VTP == 1
            S_VTP_s = Geo_tier.S_VTP_s;
            m_VTP = f_f*rho_VTP*S_VTP_s;                    %Peso HTP
        else
            m_VTP = 0;
        end
        if Can == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_Can = f_f*rho_w*S_w2_s;                  %Peso canard
        else
            m_Can = 0;
        end
        if Vee == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_Vee = f_f*rho_Vee*S_w2_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT;                      %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_fus_fairing;
        
        % Estimation
        m_prop = m_prop*rho_engine; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        MTOW_estimation = m_empty_estimation + m_payload + m_energy;
        m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear;                        %Peso tren de aterrizaje 
        
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;
        if case_AC == 6
            m_empty = 7.9;
        end
        m_TOW = m_empty + m_payload + m_energy;
    case 3 % Factores Lineales
        rho_w = 49; % factor lineal wing
        rho_HTP = 27; % factor lineal htp
        rho_VTP = 27; % factor lineal wing
        rho_fus_fairing = 24; % factor lineal fuselage
        rho_landing_gear = 0.057; % factor lineal landing gear
        rho_engine = 1.4; % factor lineal engine installation
        rho_misc = 0.10; % factor lineal msc
        
        if W1 == 1
            S_w1_s = Geo_tier.S_w1_s;
            m_w1 = f_f*rho_w*S_w1_s;                             % Wing weight.
        else
            m_w1 = 0;
        end
        if HTP == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_HTP = f_f*rho_HTP*S_w2_s;                    %Peso HTP
        else
            m_HTP = 0;
        end
        if VTP == 1
            S_VTP_s = Geo_tier.S_VTP_s;
            m_VTP = f_f*rho_VTP*S_VTP_s;                    %Peso HTP
        else
            m_VTP = 0;
        end
        if Can == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_Can = f_f*rho_w*S_w2_s;                  %Peso canard
        else
            m_Can = 0;
        end
        if Vee == 1
            S_w2_s = Geo_tier.S_w2_s;
            m_Vee = f_f*rho_Vee*S_w2_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT;                      %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_fus_fairing;
        
        % Estimation
        m_prop = m_prop*rho_engine; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        
        MTOW_estimation = m_empty_estimation + m_payload + m_energy;
        m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear;                        %Peso tren de aterrizaje 
        
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;
        if case_AC == 6
            m_empty = 7.9;
        end
        
        m_TOW = m_empty + m_payload + m_energy;
end
      
Weight_tier.m_TOW = m_TOW;
Weight_tier.m_energy = m_energy;
Weight_tier.m_empty = m_empty;
Weight_tier.m_payload = m_payload;
Weight_tier.m_f_W0 = m_energy/m_TOW;
Weight_tier.m_e_W0 = m_empty/m_TOW;
Weight_tier.M_PAYLOAD.m_pax_crew = m_pax_crew;
Weight_tier.M_PAYLOAD.m_cargo = m_cargo;
Weight_tier.M_ENERGY.m_batteries = m_batteries;
Weight_tier.M_ENERGY.m_fuel = m_fuel;
Weight_tier.M_EMPTY.m_estructure = m_estructure;
Weight_tier.M_EMPTY.m_prop = m_prop;
Weight_tier.M_EMPTY.m_landing_gear = m_landing_gear;
Weight_tier.M_EMPTY.m_subsystems = m_subsystems;
Weight_tier.M_EMPTY.m_systems = m_systems;
Weight_tier.M_ESTRUCTURE.m_w1 = m_w1;
Weight_tier.M_ESTRUCTURE.m_HTP = m_HTP;
Weight_tier.M_ESTRUCTURE.m_VTP = m_VTP;
Weight_tier.M_ESTRUCTURE.m_Can = m_Can;
Weight_tier.M_ESTRUCTURE.m_Vee = m_Vee;
Weight_tier.M_ESTRUCTURE.m_fus_fairing = m_fus_fairing;
Weight_tier.M_PROP.m_engs = m_engs;
Weight_tier.M_PROP.m_ESCs = m_ESCs;
Weight_tier.M_PROP.m_props = m_props;
Weight_tier.M_SUBSYSTEMS.m_servos = m_servos;
Weight_tier.M_SUBSYSTEMS.m_wiring = m_wiring;
Weight_tier.M_SUBSYSTEMS.m_metal = m_metal;

% % Real Batteries LiFePO4
% m_cell_A123 = 76/1000;
% m_battery_S72P = m_cell_A123*14;
% SE_LiFePO4_exp = 295/m_battery_S72P;
% SE_LiFePO4_exp = 580;

Weight_tier.M_ENERGY.SE_LiFePO4 = 160; % Wh/kg Specific Energy  90–160 Wh/kg (320–580 J/g or kJ/kg)
Weight_tier.M_ENERGY.SE_LiPo = 200; % Wh/kg Specific Energy
Weight_tier.M_ENERGY.SE_FuelCells = 333; % Wh/kg Specific Energy for fuel cells
Weight_tier.M_ENERGY.SP_motor = 5; % Motor Specific Power kW/kg - TIER 1
Weight_tier.M_ENERGY.SP_ESC = 20; % Motor Specific Power kW/kg - TIER 1
Weight_tier.M_ENERGY.m_prop_din = 0.00385; % Motor Specific Power kW/kg - TIER 1

% Weight_tier.M_ENERGY.SP_motor = 8; % Motor Specific Power kW/kg - TIER 2
% Weight_tier.M_ENERGY.SP_ESC = 25; % Motor Specific Power kW/kg - TIER 2
% Weight_tier.M_ENERGY.SP_motor = 11; % Motor Specific Power kW/kg - TIER 3
% Weight_tier.M_ENERGY.SP_ESC = 30; % Motor Specific Power kW/kg - TIER 3

%% EStimation Inertia moments
L=Geo_tier.l_fus;
b=Geo_tier.b_w1;
W=m_TOW;
% Radii of Gyration for twin engine prop
Rx=0.34;
Ry=0.29;
Rz=0.44;

e = (b+L)/2;
kg2lbs = 2.204623;
m2ft = 3.28084;
slft2_2_kgm2 = 1.355817962; 
g_imp = 32.174; %[ft/s^2]

W_IMP = W*kg2lbs;
L_IMP = L*m2ft;
b_IMP = b*m2ft;
e_IMP = (b_IMP+L_IMP)/2;

% Ixx=b^2*W*Rx^2/4/1000
Ixx_IMP =(W_IMP/g_imp)*(Rx*b_IMP/2)^2;
Ixx_SI = Ixx_IMP*slft2_2_kgm2;
Ixx = Ixx_SI;

% Iyy=L^2*W*Ry^2/4/1000
Iyy_IMP = (W_IMP/g_imp)*(Ry*L_IMP/2)^2;
Iyy_SI = Iyy_IMP*slft2_2_kgm2;
Iyy = Iyy_SI;

% Izz=((b+L)/2)^2*W*Rz^2/4/1000
Izz_IMP = (W_IMP/g_imp)*(Rz*e_IMP/2)^2;
Izz_SI = Izz_IMP*slft2_2_kgm2;
Izz = Izz_SI;

%% Céfiro Inertias
Ixx=2.727;                       %inercia en eje x
Izz=9.95;                        %inercia eje z
Ixz=0.002;                     %inercia eje xz
Iyy=7.447;                      %inercia sin adimensionalizar Iy=7.447

Weight_tier.Ixx = Ixx;
Weight_tier.Iyy = Iyy;
Weight_tier.Izz = Izz;
Weight_tier.Ixz = Ixz;