function Weight_tier = Generation_Weight_Data(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,OUTPUT_read_XLSX)
% Defines Estimation of Weights according to Cefiro III densities

Weight_Estimation = OUTPUT_read_XLSX.AC_Data_flags.Weight_Estimation;

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
case_AC = AC_CONFIGURATION.case_AC;

g = conv_UNITS.g;

% Fudge factor
f_f = OUTPUT_read_XLSX.Weights_flags.f_f;
% Scaling Factor
SF = OUTPUT_read_XLSX.AC_Data_flags.SF;


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

% Fraction of weights for fuel and payload
fraction_MF = OUTPUT_read_XLSX.Weights_flags.fraction_MF; %
fraction_PL = OUTPUT_read_XLSX.Weights_flags.fraction_PL; %

% Weight of propulsion systems
W_eng = OUTPUT_read_XLSX.Weights_flags.W_eng; %
m_engs = n_eng*W_eng;
m_prop = OUTPUT_read_XLSX.Weights_flags.m_prop;
m_props = n_eng*m_prop;
m_ESC = OUTPUT_read_XLSX.Weights_flags.m_ESC;
m_ESCs = n_eng*m_ESC;
m_prop = m_engs + m_props + m_ESCs;
m_prop_fairing = OUTPUT_read_XLSX.Weights_flags.m_prop_fairing;
m_prop_fairings = n_eng*m_prop_fairing;
    
% subsystems
n_servos = OUTPUT_read_XLSX.Weights_flags.n_servos;
m_servo = OUTPUT_read_XLSX.Weights_flags.m_servo;
m_servos = m_servo*n_servos;
m_wiring = OUTPUT_read_XLSX.Weights_flags.m_wiring;
m_metal = OUTPUT_read_XLSX.Weights_flags.m_metal;
m_subsystems = m_servos + m_wiring + m_metal;
m_systems = OUTPUT_read_XLSX.Weights_flags.m_systems;
m_landing_gear = OUTPUT_read_XLSX.Weights_flags.m_landing_gear;

% Weight of structure
m_w1 = OUTPUT_read_XLSX.Weights_flags.m_w1;
m_HTP = OUTPUT_read_XLSX.Weights_flags.m_HTP;
m_VTP = OUTPUT_read_XLSX.Weights_flags.m_VTP;
m_Can = OUTPUT_read_XLSX.Weights_flags.m_Can;
m_Vee = OUTPUT_read_XLSX.Weights_flags.m_Vee;
m_Vee2 = OUTPUT_read_XLSX.Weights_flags.m_Vee2;
m_fus_fairing = OUTPUT_read_XLSX.Weights_flags.m_fus_fairing;

% Weight of PAX
n_pax = OUTPUT_read_XLSX.Weights_flags.n_pax;
n_crew = OUTPUT_read_XLSX.Weights_flags.n_crew;
m_pax = OUTPUT_read_XLSX.Weights_flags.m_pax; % mass for passanger
m_crew = OUTPUT_read_XLSX.Weights_flags.m_crew; % mass for crew
m_luggage = OUTPUT_read_XLSX.Weights_flags.m_luggage; % mass for luggage and passanger
% m_pax_crew = (n_crew + n_pax)*m_luggage + n_crew*m_crew + n_pax*m_pax;
m_crew = n_crew*m_luggage + n_crew*m_crew;
m_pax = (n_pax)*m_luggage + n_pax*m_pax;

% Weight of Payload
m_cargo = OUTPUT_read_XLSX.Weights_flags.m_cargo;
m_payload = m_cargo;

% Energy (fuel or batteries)
m_batteries = OUTPUT_read_XLSX.Weights_flags.m_batteries;
m_fuel = OUTPUT_read_XLSX.Weights_flags.m_fuel;
m_energy = (m_batteries + m_fuel);

flag_landing_gear = OUTPUT_read_XLSX.Weights_flags.flag_landing_gear; % determines if landing gear present (1 yes 0 no)
    
% [rho_f,rho_fairing,rho_w,rho_HTP,rho_VTP,rho_tb] = CFIIIdensities; % calculated from Céfiro III - Composites
% [rho_f,rho_fairing,rho_nose,rho_w,rho_HTP,rho_VTP,rho_tb]=CFIdensities; % calculated from Céfiro I
% [rho_fus_fairing,rho_w,rho_can,rho_HTP,rho_VTP,rho_Vee]=EMERGENTIadensities(AC_CONFIGURATION,OUTPUT_read_XLSX,Body_Geo); % calculated from ProVANT-EMERGENTIA 100% I

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
            S_HTP_s = Geo_tier.S_HTP_s;
            m_HTP = f_f*rho_HTP*S_HTP_s;                    %Peso HTP
        else
            m_HTP = 0;
        end

        if VTP == 1
            S_VTP_s = Geo_tier.S_VTP_s;
            m_VTP = f_f*rho_VTP*S_VTP_s;                    %Peso VTP
        else
            m_VTP = 0;
        end

        if Can == 1
            S_vee_s = Geo_tier.S_vee_s;
            m_Can = f_f*rho_w*S_vee_s;                  %Peso canard
        else
            m_Can = 0;
        end
        
        if Vee == 1
            S_vee_s = Geo_tier.S_vee_s;
            m_Vee = f_f*rho_Vee*S_vee_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        if Vee2 == 1
            S_vee2_s = Geo_tier.S_vee2_s;
            m_Vee2 = f_f*rho_Vee*S_vee2_s;                 %Peso cola en v
        else
            m_Vee2 = 0;
        end

        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT;                      %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_Vee2 + m_fus_fairing;
                
        % Estimation
        m_prop = m_prop*rho_engine; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        MTOW_estimation = m_empty_estimation + m_payload + m_energy + m_crew + m_pax ;
        m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear;                        %Peso tren de aterrizaje
        
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;
        if case_AC == 6
            m_empty = 7.9;
        end
        
        m_TOW = m_empty + m_payload + m_energy*fraction_MF + m_crew + m_pax*fraction_PL;

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
            S_HTP_s = Geo_tier.S_HTP_s;
            m_HTP = f_f*rho_HTP*S_HTP_s;                    %Peso HTP
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
            S_vee_s = Geo_tier.S_vee_s;
            m_Can = f_f*rho_w*S_vee_s;                  %Peso canard
        else
            m_Can = 0;
        end
        
        if Vee == 1
            S_vee_s = Geo_tier.S_vee_s;
            m_Vee = f_f*rho_Vee*S_vee_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        if Vee2 == 1
            S_vee2_s = Geo_tier.S_vee2_s;
            m_Vee2 = f_f*rho_Vee*S_vee2_s;                 %Peso cola en v
        else
            m_Vee2 = 0;
        end

        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT;                      %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_Vee2 + m_fus_fairing;
        
        % Estimation
        m_prop = m_prop*rho_engine; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        MTOW_estimation = m_empty_estimation + m_payload + m_energy + m_crew + m_pax;
        m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear;                        %Peso tren de aterrizaje 
        
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;
        if case_AC == 6
            m_empty = 7.9;
        end
    
        m_TOW = m_empty + m_payload + m_energy*fraction_MF + m_crew + m_pax*fraction_PL;

    case 3 % Factores Lineales
        rho_w = 49; % factor lineal wing
        rho_HTP = 27; % factor lineal htp
        rho_VTP = 27; % factor lineal wing
        rho_fus_fairing = 24; % factor lineal fuselage
        rho_landing_gear = 0.057; % factor lineal landing gear
        rho_engine = 1.4; % factor lineal engine installation
        rho_misc = 0.10; % factor lineal msc
        rho_Vee = (rho_HTP + rho_VTP)/2;

        if W1 == 1
            S_w1_s = Geo_tier.S_w1_s;
            m_w1 = f_f*rho_w*S_w1_s;                             % Wing weight.
        else
            m_w1 = 0;
        end
        
        if HTP == 1
            S_HTP_s = Geo_tier.S_HTP_s;
            m_HTP = f_f*rho_HTP*S_HTP_s;                    %Peso HTP
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
            S_vee_s = Geo_tier.S_vee_s;
            m_Can = f_f*rho_w*S_vee_s;                  %Peso canard
        else
            m_Can = 0;
        end
        
        if Vee == 1
            S_vee_s = Geo_tier.S_vee_s;
            m_Vee = f_f*rho_Vee*S_vee_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        if Vee2 == 1
            S_vee2_s = Geo_tier.S_vee2_s;
            m_Vee2 = f_f*rho_Vee*S_vee2_s;                 %Peso cola en v
        else
            m_Vee2 = 0;
        end

        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT;                      %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_Vee2 + m_fus_fairing;
        
        % Estimation
        m_prop = m_prop*rho_engine; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        
        MTOW_estimation = m_empty_estimation + m_payload + m_energy + m_crew + m_pax;
        m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear;                        %Peso tren de aterrizaje 
        
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;
        if case_AC == 6
            m_empty = 7.9;
        end
        m_TOW = m_empty + m_payload + m_energy*fraction_MF + m_crew + m_pax*fraction_PL;
    
    case 4 % Exact Estimation
        [rho_fus_fairing,rho_w,rho_can,rho_HTP,rho_VTP,rho_Vee] = EMERGENTIadensities(AC_CONFIGURATION,OUTPUT_read_XLSX,Body_Geo,Geo_tier); % calculated from ProVANT-EMERGENTIA 100% I
        rho_landing_gear = 0.05; % factor lineal landing gear
        rho_engine = 1; % factor lineal engine installation
        rho_misc = 0.10; % factor lineal msc

        if W1 == 1
            m_w1 = f_f*rho_w*Geo_tier.S_w1_s;                             % Wing weight.
        else
            m_w1 = 0;
        end

        if HTP == 1
            m_HTP = f_f*rho_HTP*Geo_tier.S_HTP_s;                    %Peso HTP
        else
            m_HTP = 0;
        end

        if VTP == 1
            m_VTP = f_f*rho_VTP*Geo_tier.S_VTP_s;                    %Peso HTP
        else
            m_VTP = 0;
        end

        if Can == 1
            m_Can = f_f*rho_can*Geo_tier.S_can_s;                  %Peso canard
        else
            m_Can = 0;
        end

        if Vee == 1
            m_Vee = f_f*rho_Vee*Geo_tier.S_vee_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        if Vee2 == 1
            m_Vee2 = f_f*rho_Vee*Geo_tier.S_vee2_s;                 %Peso cola en v
        else
            m_Vee2 = 0;
        end

        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT;                      %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_Vee2 + m_fus_fairing;
                
        % Estimation
        m_prop = m_props + m_prop_fairings; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        MTOW_estimation = m_empty_estimation + m_payload + m_energy + m_crew + m_pax;
        m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear;                        %Peso tren de aterrizaje
        
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;
        if case_AC == 6
            m_empty = 7.9;
        end
        m_TOW = m_empty + m_payload + m_energy*fraction_MF + m_crew + m_pax*fraction_PL;

    case 5 % Mix Jaleo & Cefiro III
        
        [rho_f,rho_fairing,rho_w,rho_HTP,rho_VTP,rho_tb]=CFIVdensities; % calculated from Céfiro III
        rho_Vee = (rho_HTP + rho_VTP)/2;
        m_fus = rho_f*Surf_TOT;               % Fuselage weight
        m_fairing = rho_fairing*Surf_TOT;               % Fairing weight
        rho_fus_fairing = (rho_f + rho_fairing)/2;
        rho_landing_gear = 0.04; % factor lineal landing gear
        rho_engine = 1.1; % factor lineal engine installation
        rho_misc = 0.10; % factor lineal msc
         
        if W1 == 1
            S_w1_s = Geo_tier.S_w1_s;
            m_w1 = f_f*rho_w*S_w1_s;                             % Wing weight.
        else
            m_w1 = 0;
        end

        if HTP == 1
            S_HTP_s = Geo_tier.S_HTP_s;
            m_HTP = f_f*rho_HTP*S_HTP_s;                    %Peso HTP
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
            S_can_s = Geo_tier.S_can_s;
            m_Can = f_f*rho_w*S_can_s;                  %Peso canard
        else
            m_Can = 0;
        end

        if Vee == 1
            S_vee_s = Geo_tier.S_vee_s;
            m_Vee = f_f*rho_Vee*S_vee_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        if Vee2 == 1
            S_vee2_s = Geo_tier.S_vee2_s;
            m_Vee2 = f_f*rho_Vee*S_vee2_s;                 %Peso cola en v
        else
            m_Vee2 = 0;
        end

        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT;                      %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_Vee2 + m_fus_fairing;
                
        % Estimation
        m_prop = m_prop*rho_engine; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        MTOW_estimation = m_empty_estimation + m_payload + m_energy + m_crew + m_pax ;
        % m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        % m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear;                        %Peso tren de aterrizaje
        
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;

        if case_AC == 6
            m_empty = 7.9;
        end
        
        m_TOW = m_empty + m_payload + m_energy*fraction_MF + m_crew + m_pax*fraction_PL;

end
      
if OUTPUT_read_XLSX.Weights_flags.flag_total_weights == 1
     m_empty = OUTPUT_read_XLSX.Weights_flags.ME_true;
     m_payload = OUTPUT_read_XLSX.Weights_flags.MPL_true*fraction_PL;
     m_energy = OUTPUT_read_XLSX.Weights_flags.MF_true*fraction_MF;
     m_crew = OUTPUT_read_XLSX.Weights_flags.MCREW_true;
     m_pax = 0; %Does not differentiate between cargo and pax
     m_TOW = m_empty + m_payload + m_energy + m_crew + m_pax;
end

% General Storing Weights
Weight_tier.m_TOW = m_TOW;
Weight_tier.m_empty = m_empty;
Weight_tier.m_payload = m_payload;
Weight_tier.m_energy = m_energy;
Weight_tier.m_crew = m_crew;
Weight_tier.m_pax = m_pax;

Weight_tier.m_f_W0 = m_energy/m_TOW;
Weight_tier.m_e_W0 = m_empty/m_TOW;

Weight_tier.M_PAYLOAD.m_pax = m_pax;
Weight_tier.M_PAYLOAD.m_crew = m_crew;

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
Weight_tier.M_ESTRUCTURE.m_Vee2 = m_Vee2;
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

% Weight_tier.M_ENERGY.SE_LiFePO4 = 160; % Wh/kg Specific Energy  90–160 Wh/kg (320–580 J/g or kJ/kg)
% Weight_tier.M_ENERGY.SE_LiPo = 200; % Wh/kg Specific Energy
% Weight_tier.M_ENERGY.SE_FuelCells = 333; % Wh/kg Specific Energy for fuel cells
% Weight_tier.M_ENERGY.SP_motor = 5; % Motor Specific Power kW/kg - TIER 1
% Weight_tier.M_ENERGY.SP_ESC = 20; % Motor Specific Power kW/kg - TIER 1
% Weight_tier.M_ENERGY.m_prop_din = 0.00385; % Motor Specific Power kW/kg - TIER 1

Weight_tier.M_ENERGY.SE_LiFePO4 = OUTPUT_read_XLSX.Propulsive_flags.SE_LiFePO4; % Wh/kg Specific Energy  90–160 Wh/kg (320–580 J/g or kJ/kg)
Weight_tier.M_ENERGY.SE_LiPo = OUTPUT_read_XLSX.Propulsive_flags.SE_LiPo; % Wh/kg Specific Energy
Weight_tier.M_ENERGY.SE_FuelCells = OUTPUT_read_XLSX.Propulsive_flags.SE_FuelCells; % Wh/kg Specific Energy for fuel cells
Weight_tier.M_ENERGY.SP_motor = OUTPUT_read_XLSX.Propulsive_flags.SP_motor; % Motor Specific Power kW/kg - TIER 1
Weight_tier.M_ENERGY.SP_ESC = OUTPUT_read_XLSX.Propulsive_flags.SP_ESC; % Motor Specific Power kW/kg - TIER 1
Weight_tier.M_ENERGY.m_prop_din = OUTPUT_read_XLSX.Propulsive_flags.m_prop_din; % Motor Specific Power kW/kg - TIER 1

% Weight_tier.M_ENERGY.SP_motor = 8; % Motor Specific Power kW/kg - TIER 2
% Weight_tier.M_ENERGY.SP_ESC = 25; % Motor Specific Power kW/kg - TIER 2
% Weight_tier.M_ENERGY.SP_motor = 11; % Motor Specific Power kW/kg - TIER 3
% Weight_tier.M_ENERGY.SP_ESC = 30; % Motor Specific Power kW/kg - TIER 3


% Determines if the Moments of Inertia are obtained via literature or given
% by the user
if OUTPUT_read_XLSX.Weights_flags.flag_moments_inertia == 1
    %% EStimation Inertia moments
    L=Geo_tier.l_fus;
    b=Geo_tier.b_w1;
    W=m_TOW;
    % Radii of Gyration for twin engine prop
    Rx=0.34;
    Ry=0.29;
    Rz=0.44;
    
    e = (b+L)/2;
    
    % Conversion units
    kg2lbs = conv_UNITS.kg2lbs;
    m2ft = conv_UNITS.m2ft;
    slft2_2_kgm2 = conv_UNITS.slft2_2_kgm2;
    g_imp = conv_UNITS.g_imp; %[ft/s^2] = 32.174; %[ft/s^2]
    
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

    Weight_tier.Ixx = Ixx;
    Weight_tier.Iyy = Iyy;
    Weight_tier.Izz = Izz;
    Weight_tier.Ixy = 0; % Assume zero for the estimation
    Weight_tier.Ixz = 0; % Assume zero for the estimation
    Weight_tier.Iyz = 0; % Assume zero for the estimation
    
else %% Exact Inertias
    
    Weight_tier.Ixx = OUTPUT_read_XLSX.Weights_flags.I_xx;
    Weight_tier.Iyy = OUTPUT_read_XLSX.Weights_flags.I_yy;
    Weight_tier.Izz = OUTPUT_read_XLSX.Weights_flags.I_zz;
    Weight_tier.Ixy = OUTPUT_read_XLSX.Weights_flags.I_xy;
    Weight_tier.Ixz = OUTPUT_read_XLSX.Weights_flags.I_xz;
    Weight_tier.Iyz = OUTPUT_read_XLSX.Weights_flags.I_yz;
end
