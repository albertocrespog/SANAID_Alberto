function [conditions OUTPUT_read_XLSX] = mission_changes_TARSIS(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

Weight_mission = T120_weightconfiguration;

conditions.m_TOW = Weight_mission.W_vec_f(1);

%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

% Percentage where bayoneta is located
pp_bayoneta = (1100.55-966.57)/469;
cR_w1_T75 = 0.421;
x_loc_1R_y1_w1_CAD_T75 = 1.45766;
x_loc_1R_y1_w1_CAD_T75_bayoneta = x_loc_1R_y1_w1_CAD_T75 + pp_bayoneta*cR_w1_T75;

x_loc_1R_y1_w2_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD;

conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

CATIAT120 = xlsread('Tabla extraida de CATIA.xlsx','D4:M25');
OUTPUT_read_XLSX.Weights_flags.CATIAT120 = CATIAT120;

% Conversion factor depending on the units of the CG_TARSIS120 (1 for m 1000 for mm)
% conv_factor = 1 old version
% conv_factor = 1000 new version
conv_factor = 1000;
OUTPUT_read_XLSX.Weights_flags.conv_factor = conv_factor;

% Extra manual modification of Xcg
Dxcg = Modifications.Dxcg;

% Weight Configuration from AERTEC
conditions.MOTOR = 1; % 1 para SP210, 2 para DA215
conditions.RACK = 1; % 1 == rack, 0 == sin rack
conditions.DEPOSITO = 1; % 1 == deposito original, 2 == depósito húmedo
conditions.FUEL = 1;  % en tanto por 1
conditions.W_PL = 3.75;

% Número de misiles
conditions.n_MSL = 4; %0,1,2,3,4

OUTPUT = CG_TARSIS120(conditions.FUEL,conditions.W_PL,conditions.RACK,conditions.n_MSL,conditions.DEPOSITO,conditions.MOTOR,CATIAT120);
        
%% Modification that include interference with RACK
y_loc_1R_y2_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD;
y_loc_1R_y1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD;
b_w1_e2 =  y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CAD;
y_offset_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD;

y_1R_y2_ail = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD;
y_1R_y1_ail =  y_1R_y2_ail - 671/1000;
K_y1_ail_w1 = (y_1R_y1_ail - y_offset_w1)/(b_w1_e2);
K_y2_ail_w1 = (y_1R_y2_ail - y_offset_w1)/(b_w1_e2);

OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_ail_w1 = K_y1_ail_w1;
OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_ail_w1 = K_y2_ail_w1;

Racks_wing_span = 0.00;
OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD + Racks_wing_span;
y_offset_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD;
b_w1_e2 = (OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD - OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD);
K_y1_ail_w1 = (y_1R_y1_ail - y_offset_w1)/(b_w1_e2);
K_y2_ail_w1 = (y_1R_y2_ail - y_offset_w1)/(b_w1_e2);
OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y1_ail_w1 = K_y1_ail_w1;
OUTPUT_read_XLSX.InputGeometry_Data_flags.K_y2_ail_w1 = K_y2_ail_w1;

% Identifies folder where to store Figures
filename = '../Results/TARSIS_75/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;


switch mission_actual
    case 0 % ID_01. c_nom, b	
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.421;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.421;
%       OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 0.421;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1; 
        
        % Plots prefix
        prefix = 'CASO_ID_00_TARSIS_';
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


    case 1 % ID_01. c_nom, b+20	 
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 20/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.421;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.421;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 4;     
        % Plots prefix
        prefix = 'CASO_ID_01_TARSIS_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix

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


    case 2 % ID_02. c_nom, b+30	 
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 30/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.421;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.421;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 5; 
        % Plots prefix
        prefix = 'CASO_ID_02_TARSIS_'
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix

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

    case 3 % ID_03. c_nom, b+50	 
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 50/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.421;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.421;
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 6; 
        % Plots prefix
        prefix = 'CASO_ID_03_TARSIS_';
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

    case 4 % ID_04. c_a, b+30	 
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 30/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.469;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.469;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 7; 
        % Aerodynamics
        DeltaCR = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 - cR_w1_T75;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD_T75_bayoneta - pp_bayoneta*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
        % Plots prefix
        prefix = 'CASO_ID_04_TARSIS_';
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

    case 5 % ID_05. c_a, b+50	 
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 50/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.469;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.469;
        
        DeltaCR = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 - cR_w1_T75;
        xCG = 1.689;
%         x_loc_1R_y1_w1_CAD_T75_bayoneta - pp_bayoneta*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1
%         xCG - OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1/4 - 0.1
        location_XAC_wing = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD + OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1/4;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD_T75_bayoneta - pp_bayoneta*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = xCG - OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1/4;
%         Delta_HTP = location_XAC_wing - OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD;
%         OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD = x_loc_1R_y1_w2_CAD + 0.00;
        
        % Correction of geometry according to CATIA - AERTEC
        % Wing Corrections
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = 1447.697/1000;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_w1_CAD = -57.955/1000;
        % HTP Corrections 
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w2_CAD = 3479.874/1000;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_w2_CAD = 99.816/1000;
        % VTP1 Corrections
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = 3392.884/1000;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = -92.192/1000;
        % VTP2 Corrections
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP2_CAD = 3392.884/1000;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP2_CAD = -92.192/1000;
        
        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8; 
        % Plots prefix      
        
        prefix = 'CASO_ID_05_TARSIS_';
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

    case 6 % ID_06. c_a, b+70	 
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 70/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.469;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.469;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 9; 
        % Aerodynamics
        DeltaCR = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 - cR_w1_T75;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD_T75_bayoneta - pp_bayoneta*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
        % Plots prefix
        prefix = 'CASO_ID_06_TARSIS_';
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

    case 7 % ID_07. c_b, b+50
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 50/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.501;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 10; 
        % Aerodynamics
        DeltaCR = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 - cR_w1_T75;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD_T75_bayoneta - pp_bayoneta*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
        % Plots prefix
        prefix = 'CASO_ID_07_TARSIS_';
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

    case 8 % ID_08. c_b, b+70	 
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 70/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.501;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 11; 
        % Aerodynamics
        DeltaCR = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 - cR_w1_T75;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD_T75_bayoneta - pp_bayoneta*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
        % Plots prefix
        prefix = 'CASO_ID_08_TARSIS_';
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

    case 9 % ID_09. c_b, b+90	 
        % Mass and Inertia: Assumed to be those of T120 with missiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);       
        % Geometry
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = 5.22/2 + 90/100;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 = 0.501;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.cT_w1 = 0.501;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 12; 
        % Aerodynamics
        DeltaCR = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1 - cR_w1_T75;
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD_T75_bayoneta - pp_bayoneta*OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
        % Plots prefix
        prefix = 'CASO_ID_09_TARSIS_';
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
