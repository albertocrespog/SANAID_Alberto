function [conditions OUTPUT_read_XLSX] = mission_changes_A320_Performance_study_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS)
             
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
        prefix = 'CASO_ID_00_A320_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
       
        % Performance
        OUTPUT_read_XLSX.Performance_pre_flags.V = 164;
        OUTPUT_read_XLSX.Performance_pre_flags.h = 1000;



        case 1 % ID_00 - Open VSP
        % Geometry                
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{1}(4);
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{1}(4)*0.2949;
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{1}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 4;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 13;     
        % Plots prefix
        prefix = 'CASO_ID_01_A320_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
        % Performance
        OUTPUT_read_XLSX.Performance_pre_flags.V = 164;
        OUTPUT_read_XLSX.Performance_pre_flags.h = 3048;

       
    case 2 % ID_01
        % Geometry                
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{1}(4);
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{1}(4)*0.2949;
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{1}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 4;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 13;     
        % Plots prefix
        prefix = 'CASO_ID_02_A320_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
        % Performance
        OUTPUT_read_XLSX.Performance_pre_flags.V = 90.0278;
        OUTPUT_read_XLSX.Performance_pre_flags.h = 3048;

 
    case 3 % ID_01
        % Geometry                
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{1}(4);
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{1}(4)*0.2949;
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{1}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 4;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 13;     
        % Plots prefix
        prefix = 'CASO_ID_03_A320_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
        % Performance
        OUTPUT_read_XLSX.Performance_pre_flags.V = 149.1889;
        OUTPUT_read_XLSX.Performance_pre_flags.h = 3048;

    case 4 % ID_01
        % Geometry                
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2 = Geometry_changes.Sheet_case{1}(4);
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2 = Geometry_changes.Sheet_case{1}(4)*0.2949;
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e = conv_UNITS.D2R*Geometry_changes.Sheet_case{1}(1);        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 4;     
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2_VLM = 13;     
        % Plots prefix
        prefix = 'CASO_ID_04_A320_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
        % Performance
        OUTPUT_read_XLSX.Performance_pre_flags.V = 231.5;
        OUTPUT_read_XLSX.Performance_pre_flags.h = 3048;

end