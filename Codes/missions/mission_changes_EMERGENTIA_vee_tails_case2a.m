function [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_vee_tails_case2a(OUTPUT_read_XLSX,mission_actual,conv_UNITS)
             
%Loads geometry changes
Geometry_changes = EMERGENTIA_vee_tails_geometry_v0;
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

% Identifies folder where to store Figures
filename = '../Results/EMERGENTIA_Manufactured/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

switch mission_actual
    case 1 % ID_01
        % Geometry
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.00;
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = Geometry_changes.Sheet_case{7}(1); %y0
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{8}(1); %y1
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e     = Geometry_changes.Sheet_case{9}(1)*pi/180; %Lambda
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e      = Geometry_changes.Sheet_case{10}(1)*pi/180; %Gamma
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2              = Geometry_changes.Sheet_case{11}(1); %cr
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2              = Geometry_changes.Sheet_case{12}(1); %cr
        
        % Aerodynamics
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        %         OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 21; %diedro 35º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 34;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_01_EMERGENTIA_';
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
    
    case 2 % ID_02
        % Geometry
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.00;
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = Geometry_changes.Sheet_case{7}(2); %y0
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{8}(2); %y1
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e     = Geometry_changes.Sheet_case{9}(2)*pi/180; %Lambda
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e      = Geometry_changes.Sheet_case{10}(2)*pi/180; %Gamma
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2              = Geometry_changes.Sheet_case{11}(2); %cr
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2              = Geometry_changes.Sheet_case{12}(2); %cr
    
        % Aerodynamics
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        %         OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 24; %diedro 40º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 37;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_02_EMERGENTIA_';
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


    case 3 % ID_03
        % Geometry
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.00;
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = Geometry_changes.Sheet_case{7}(3); %y0
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{8}(3); %y1
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e     = Geometry_changes.Sheet_case{9}(3)*pi/180; %Lambda
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e      = Geometry_changes.Sheet_case{10}(3)*pi/180; %Gamma
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2              = Geometry_changes.Sheet_case{11}(3); %cr
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2              = Geometry_changes.Sheet_case{12}(3); %cr
        
%         Aerodynamics
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        %         OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 26; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 39;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_03_EMERGENTIA_';
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

    case 4 % ID_04
        % Geometry
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.00;
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = Geometry_changes.Sheet_case{7}(4); %y0
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{8}(4); %y1
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e     = Geometry_changes.Sheet_case{9}(4)*pi/180; %Lambda
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e      = Geometry_changes.Sheet_case{10}(4)*pi/180; %Gamma
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2              = Geometry_changes.Sheet_case{11}(4); %cr
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2              = Geometry_changes.Sheet_case{12}(4); %cr
        
        % Aerodynamics
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        %         OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 28; %diedro 50º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 41;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_04_EMERGENTIA_';
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
    
    case 5 % ID_05
        % Geometry
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.00;
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = Geometry_changes.Sheet_case{7}(5); %y0
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{8}(5); %y1
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e     = Geometry_changes.Sheet_case{9}(5)*pi/180; %Lambda
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e      = Geometry_changes.Sheet_case{10}(5)*pi/180; %Gamma
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2              = Geometry_changes.Sheet_case{11}(5); %cr
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2              = Geometry_changes.Sheet_case{12}(5); %cr
        
        % Aerodynamics
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
                OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        %         OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 31; %diedro 55º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 44;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_05_EMERGENTIA_';
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

