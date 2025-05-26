function [conditions OUTPUT_read_XLSX] = mission_changes_MILVUS_STABILITY_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
            
%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;
% Folder where to stores de Figures
m2ft = conv_UNITS.m2ft;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engine OFF
% Posicion_Palanca = 1; % Engine ON
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

%% Initialize Stability Studies
conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS);

OUTPUT_read_XLSX.Fuselage_flags.Use_Storing_DATA == 1;

switch mission_actual
    case 1 % ID_01
        
        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 2;        
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 0;
        fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %

        %% Plots prefix
        prefix = 'CASO_ID_01_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 0;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Plotting
        gamma_min = -20*D2R;
        gamma_max = -20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        gamma_max = OUTPUT_read_XLSX.Stability_flags.gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

    case 2 % ID_01
        
        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 2;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 0;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_02_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Plotting
        
    case 3 % ID_01
        
        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 2;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 0;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_03_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 0;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Plotting

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.495914669026212;

    case 4 % ID_01
        
        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 2;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 0;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_04_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;


    case 5 % ID_01
        
        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 6;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 0;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_05_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 0;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.575;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.575;

    case 6 % ID_01
        
        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 6;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 0;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_06_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.575;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.575;
    case 7 % ID_01

        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 0;
        fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_07_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 0;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.575;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.5;

    case 8 % ID_01
        
        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 0;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_08_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.575;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.5;


    case 9 % ID_01

        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 6;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 1;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_09_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 0;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.575;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.575;

    case 10 % ID_01

        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 6;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 1;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_10_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.575;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.575;
    case 11 % ID_01

        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 1;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_11_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 0;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.575;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.5;

    case 12 % ID_01

        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 4;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 12;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_12_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.575;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.5;

            case 13 % ID_01

        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 6;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 5;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 6;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 7;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 8;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 1;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %


        %% Plots prefix
        prefix = 'CASO_ID_13_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.50;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.50;
        % Extended Vtail     
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.330*cos(45*D2R) + 0.065;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = (45*D2R);
        
          case 14 % ID_01

        %% Aerodynamic
        % wing 1
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy = 2;
        % Canard
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can = 1;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy = 2;
        % Vtail
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee = 7;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy = 8;
        % Uses estimation of lateral Derivatives from FLOW
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.FLOW5_LATERAL = 1;
                fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %

        OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_vee = -2;

        %% Plots prefix
        prefix = 'CASO_ID_14_MILVUS_STABILITY_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Stability
        OUTPUT_read_XLSX.Stability_flags.V_high = 45;
        OUTPUT_read_XLSX.Stability_flags.Munk_fuselage_constribution = 1;
        % Gamma Trim conditions
        gamma_min = -20*D2R;
        gamma_max = 20*D2R;
        N_gamma_var = 12;
        OUTPUT_read_XLSX.Stability_flags.gamma_min = gamma_min;
        OUTPUT_read_XLSX.Stability_flags.gamma_max = gamma_max;
        OUTPUT_read_XLSX.Stability_flags.N_gamma_var = N_gamma_var;

        %% Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.482355966845501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.5;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_can_CAD = 0.5;
        % Extended Vtail     
        % OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 0.330*cos(45*D2R) + 0.065;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_vee_e = (45*D2R);
        



end