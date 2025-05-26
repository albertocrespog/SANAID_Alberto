function [conditions OUTPUT_read_XLSX] = mission_changes_VANTUS_24_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
             
%Loads geometry changes
Geometry_changes = EMERGENTIA_vee_tails_geometry_v0;
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

% Identifies folder where to store Figures
filename = '../Results/VANTUS/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

switch mission_actual
    case 1 % ID_01
        % Geometry
        x_loc_1R_y1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD;       
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = x_loc_1R_y1_w1_CAD + 10/100;
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1;
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = Geometry_changes.Sheet_case{1}(1); %y0
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{2}(1); %y1
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e     = Geometry_changes.Sheet_case{3}(1)*pi/180; %Lambda
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e      = Geometry_changes.Sheet_case{4}(1)*pi/180; %Gamma
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2              = Geometry_changes.Sheet_case{5}(1); %cr
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2              = Geometry_changes.Sheet_case{6}(1); %cr
        
        % KINK Geometry
        % wingspan locations
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_w1_CAD  = 0.825; % y loc of wing (w1) kink1 chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_w1_CAD  = 1.112; % y loc of wing (w1) kink2 chord LE position (distance from CAD refference point)		
        % Sweep LE
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k1_e  = 2.80000*pi/180; % Sweep of kink1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k2_e  = 12.50000*pi/180; % Sweep of kink2 (deg)	
        % Chrord
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_w1  = 0.178; % Kink1 Chord w1
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_w1  = 0.066; % Kink2 Chord w1

        % % Performance
        % OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % OUTPUT_read_XLSX.Performance_pre_flags.V = 27;

        % Plots prefix
        prefix = 'CASO_ID_01_VANTUS_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
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
    case 0 % ID_02
        % Geometry
        
        
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.00;
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1;
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = Geometry_changes.Sheet_case{1}(1); %y0
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{2}(1); %y1
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e     = Geometry_changes.Sheet_case{3}(1)*pi/180; %Lambda
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e      = Geometry_changes.Sheet_case{4}(1)*pi/180; %Gamma
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2              = Geometry_changes.Sheet_case{5}(1); %cr
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2              = Geometry_changes.Sheet_case{6}(1); %cr
        
        % KINK Geometry
        % wingspan locations
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_w1_CAD  = 0.0; % y loc of wing (w1) kink1 chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_w1_CAD  = 0.0; % y loc of wing (w1) kink2 chord LE position (distance from CAD refference point)		
        % Sweep LE
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k1_e  = 0.0; % Sweep of kink1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k2_e  = 0.0; % Sweep of kink2 (deg)	
        % Chrord
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_w1  = 0.0; % Kink1 Chord w1
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_w1  = 0.0; % Kink2 Chord w1

        % % Performance
        % OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % OUTPUT_read_XLSX.Performance_pre_flags.V = 27;

        % Plots prefix
        prefix = 'CASO_ID_02_VANTUS_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
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
end

