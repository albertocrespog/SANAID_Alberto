function [conditions OUTPUT_read_XLSX] = mission_changes_A320_tails_v0(OUTPUT_read_XLSX,mission_actual,conv_UNITS)
             
%Loadsa de geometry changes
Geometry_changes = ONEiRE_geometry_v0;
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
lambda_w2 = Geometry_changes.Sheet_case{10}(2);


%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;


% Identifies folder where to store Figures
filename = '../Results/ONEiRE/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

% Asymmetric weight in each wing
% right wing
conditions.conditions_TRIM_lat.sign1 = 1;
conditions.conditions_TRIM_lat.L_arm1 = 0;
conditions.conditions_TRIM_lat.Z_arm1 = 0;
conditions.conditions_TRIM_lat.L_weight1 = 0; %
% left wing
conditions.conditions_TRIM_lat.sign2 = -1;
conditions.conditions_TRIM_lat.L_arm2 = 0;
conditions.conditions_TRIM_lat.Z_arm2 = 0;
conditions.conditions_TRIM_lat.L_weight2 = 0;%

%% Rolling Moments definition
% CAse 1: NO rolling movement
% CAse 2: Rolling moment asymetric weight 1 different in each wing
% CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
% CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
conditions.conditions_TRIM_lat.Rolling_maneuver = 1;
conditions.conditions_TRIM_lat.change_dir_roll = 1;

%% Yawing Moments definition
% CAse 1: NO Yawing moment
% CAse 2: Yawing movement drag asymetry different in each wing
% CAse 3: Yawing movement drag asymetry  - only in wing 1
% CAse 4: Yawing movement drag asymetry  - only in wing 2
% CAse 5: Yawing movement engine asymetry
% CAse 6: Yawing movement engine asymetry and drag asymmetry
% change yawing moment direction
conditions.conditions_TRIM_lat.Yawing_maneuver = 1;
conditions.conditions_TRIM_lat.change_dir_yaw = 1;
conditions.conditions_TRIM_lat.direction_beta = 1; % defines the directión of sideslip

N_moment = 0; % Worse case scenarios % Actualizado
conditions.conditions_TRIM_lat.N_moment = N_moment;

switch mission_actual
    case 0 % ID_00 - Open VSP
        % Geometry                
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{1}(4);
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{1}(4)*0.2949;
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{1}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 4;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 13;     
        % Plots prefix
        prefix = 'CASO_ID_00_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        

    case 1 % ID_01
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD;

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{1}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{1}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{1}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{1}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 4;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 13;     
        % Plots prefix
        prefix = 'CASO_ID_01_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
    
    case 2 % ID_02
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{2}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{2}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{2}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{2}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 5;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 14;     
        % Plots prefix
        prefix = 'CASO_ID_02_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

    case 3 % ID_03
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{3}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{3}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{3}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{3}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 6;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 15;     
        % Plots prefix
        prefix = 'CASO_ID_03_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
    
    case 4 % ID_04
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{4}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{4}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{4}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{4}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 7;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 16;     
        % Plots prefix
        prefix = 'CASO_ID_04_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
  
    case 5 % ID_05
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{5}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{5}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{5}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{5}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 8;    
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 17;
        % Plots prefix
        prefix = 'CASO_ID_05_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

    case 6 % ID_06
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{6}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{6}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{6}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{6}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 9;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 18;     
        % Plots prefix
        prefix = 'CASO_ID_06_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

    case 7 % ID_07
        % Geometry
        Geometry_changes.Sheet_case{1}(1)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{7}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{7}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{7}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{7}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 19;     
        % Plots prefix
        prefix = 'CASO_ID_07_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

     case 8 % ID_08        
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{8}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{8}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{8}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{8}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 11;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 20;     
        % Plots prefix
        prefix = 'CASO_ID_08_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

    case 9 % ID_09
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{9}(3)/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{9}(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{9}(4)*lambda_w2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{9}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 12;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 21;     
        % Plots prefix
        prefix = 'CASO_ID_09_ONEiRE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

end