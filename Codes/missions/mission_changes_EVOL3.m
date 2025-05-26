function [conditions OUTPUT_read_XLSX] = mission_changes_VANTUS(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
             
%Loads geometry changes
Geometry_changes = EMERGENTIA_vee_tails_geometry_v0;
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

switch mission_actual
    case 1 % ID_01
        % Geometry
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 1.00;
        %         OUTPUT_read_XLSX.InputGeometry_Data_flags.z_0XCG = 0.1;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w2_CAD = Geometry_changes.Sheet_case{1}(1); %y0
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w2_CAD = Geometry_changes.Sheet_case{2}(1); %y1
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w2_e     = Geometry_changes.Sheet_case{3}(1)*pi/180; %Lambda
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w2_e      = Geometry_changes.Sheet_case{4}(1)*pi/180; %Gamma
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w2              = Geometry_changes.Sheet_case{5}(1); %cr
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w2              = Geometry_changes.Sheet_case{6}(1); %cr
        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1 = 2; %
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w2 = 1; %
        % OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 20; %diedro 35º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w2 = 4; %diedro 35º
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 33;
        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        OUTPUT_read_XLSX.Performance_pre_flags.V = 27;
        % Plots prefix
        prefix = 'CASO_ID_01_VANTUS_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        OUTPUT_read_XLSX.PLOT_flags.fname = 'results\VANTUS';
        
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

