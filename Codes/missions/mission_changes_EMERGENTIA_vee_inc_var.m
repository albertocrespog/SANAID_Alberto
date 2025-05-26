function [conditions OUTPUT_read_XLSX] = mission_changes_EMERGENTIA_vee_inc_var(OUTPUT_read_XLSX,mission_actual,conv_UNITS);

%Loadsa de geometry changes
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
    case 0 %
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1; %1.029384965115954
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453+ 0.05; %
        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 3; %
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_00_EMERGENTIA_';
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
    
    case 1 % ID_01
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.0402; %1.029384965115954
        %                  OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1; %1.029384965115954
        %                  OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453+ 0.005; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi; %1.0716
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 1; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        %         OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
        
        
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
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.0402; %1.0716
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = -0.1; %1.0716
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 + 0.005; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi; %1.0716
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 1; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 2; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V =27;
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
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.0402;
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 + 0.01; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 1; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 3; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
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
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.0402; %1.0716
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = 0; %1.0716
%                 OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 - 0.04; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
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
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG =1.0402; %1.0716
        %                 OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 - 0.04; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 2; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
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

    case 6 % ID_06
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.0402; %1.0716
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 + 0.02; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 3; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
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

    case 7 % ID_07
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.0402; %1.0716
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = 0.05; %1.0716
        %                 OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 - 0.04; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 3; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_07_EMERGENTIA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        OUTPUT_read_XLSX.PLOT_flags.fname = 'results\EMERGENTIA';

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
        
    case 8 % ID_08
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.0402; %1.0716
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = 0.1; %1.0716
        %                 OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 - 0.04; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 3; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 2; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_08_EMERGENTIA_';
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

    case 9 % ID_09
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.0402; %1.0716
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = 0.15; %1.0716
        %                 OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.871453 - 0.04; %
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e = 45/180*pi;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.561;
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD =0.79;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 3; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 3; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 10; %diedro 45º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 15;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_10_EMERGENTIA_';
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