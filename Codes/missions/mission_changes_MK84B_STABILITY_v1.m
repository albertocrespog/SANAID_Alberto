function [conditions OUTPUT_read_XLSX] = mission_changes_MK84B_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
             
%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
% Folder where to stores de Figures

%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

switch mission_actual
    case 1 % ID_01
        % Aerodynamics
        % Uso de NACA Report 823 para estimación de datos Vtail
        OUTPUT_read_XLSX.aero.datos_sin_diedro = 0;

        % Geometry

        z_loc_1R_y1_vee_CAD = 0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD;
        
        dihedral_vee_e = 45*pi/180;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = dihedral_vee_e; %

        % Plots prefix
        prefix = 'CASEID_01_AC_25_MK84b_STABILITY_';
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
    case 2 % ID_01
        % Aerodynamics
        % Uso de NACA Report 823 para estimación de datos Vtail
        OUTPUT_read_XLSX.aero.datos_sin_diedro = 0;

        % Geometry

        z_loc_1R_y1_vee_CAD = -0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD;
        
        dihedral_vee_e = -45*pi/180;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = dihedral_vee_e; %

        % Plots prefix
        prefix = 'CASEID_02_AC_25_MK84b_STABILITY_';
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
case 3 % ID_01
        % Aerodynamics
        % Uso de NACA Report 823 para estimación de datos Vtail
        OUTPUT_read_XLSX.aero.datos_sin_diedro = 0;

        % Geometry

        z_loc_1R_y1_vee_CAD = 0.5;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD;
        
        
        dihedral_vee_e = -45*pi/180;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = dihedral_vee_e; %

        % Plots prefix
        prefix = 'CASEID_03_AC_25_MK84b_STABILITY_';
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


    case 4 % ID_01
        % Aerodynamics
        % Uso de NACA Report 823 para estimación de datos Vtail
        OUTPUT_read_XLSX.aero.datos_sin_diedro = 1;

        % Geometry

        z_loc_1R_y1_vee_CAD = 0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD;
        
        dihedral_vee_e = 45*pi/180;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = dihedral_vee_e; %

        % Plots prefix
        prefix = 'CASEID_04_AC_25_MK84b_STABILITY_';
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
    case 5 % ID_01
        % Aerodynamics
        % Uso de NACA Report 823 para estimación de datos Vtail
        OUTPUT_read_XLSX.aero.datos_sin_diedro = 1;

        % Geometry

        z_loc_1R_y1_vee_CAD = -0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD;
        
        dihedral_vee_e = -45*pi/180;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = dihedral_vee_e; %

        % Plots prefix
        prefix = 'CASEID_05_AC_25_MK84b_STABILITY_';
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
case 6 % ID_01
        % Aerodynamics
        % Uso de NACA Report 823 para estimación de datos Vtail
        OUTPUT_read_XLSX.aero.datos_sin_diedro = 0;

        % Geometry

        z_loc_1R_y1_vee_CAD = 0.5;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_vee_CAD = z_loc_1R_y1_vee_CAD;
        
        
        dihedral_vee_e = -45*pi/180;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = dihedral_vee_e; %

        % Plots prefix
        prefix = 'CASEID_06_AC_25_MK84b_STABILITY_';
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
