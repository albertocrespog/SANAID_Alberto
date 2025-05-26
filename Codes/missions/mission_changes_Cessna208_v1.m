function [conditions OUTPUT_read_XLSX] = mission_changes_Cessna208_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
             

D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;
m2ft = conv_UNITS.m2ft;

%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;


%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;


% Defines if the trim tabs are fixed or varied within the defined cases
fixed_trim_tab = 1; % fixed_trim_tab = 1; the fixed trim tab is defined in teh cases
conditions.conditions_TRIM_lat.fixed_trim_tab = fixed_trim_tab; % Fixed value for aileron rudder trim tab
defined_GerRatios_aileron_rudder = 1; % fixed_trim_tab = 1; the fixed trim tab is defined in teh cases
conditions.conditions_TRIM_lat.defined_GerRatios_aileron_rudder = defined_GerRatios_aileron_rudder; % Defined Gear Ratios

% Gearing 
Kra = 0; % Constant that depends on the mechanical instalation of a rudder-aileron interconection spring
Kar = 0; % Constant that depends on the mechanical instalation of a aileron-rudder interconection spring
conditions.conditions_TRIM_lat.Kra = Kra; % 
conditions.conditions_TRIM_lat.Kar = Kar; % 

% Medium Aircraft (e.g., Cessna Caravan, Beechcraft King Air): 3:1 to 4:1
Ga = 0.4*m2ft; % rad/ft -> rad/m, Stick to aileron gearing ratio
Gr = 1.5*m2ft; % rad/ft -> rad/m,
conditions.conditions_TRIM_lat.Ga = Ga; % 
conditions.conditions_TRIM_lat.Gr = Gr; % 

% Identifies folder where to store Figures
filename = '../Results/Cessna208/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

switch mission_actual
    case 1 % ID_01
        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_01_Cessna208_';
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
        conditions.conditions_TRIM_lat.LT_engine = 0; % Newton-m 

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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)

        % Condition of initial Trim Lateral Directional 
        conditions.conditions_TRIM_lat.phi_0 = 0; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 0;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0; 

        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 10*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value
        
        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

    case 2

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_02_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 2500; % Newton-m for engine un PT6A-42
                
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 0; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 0;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0; 
        
        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 10*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same

    case 3

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_03_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 0; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 0;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;  
        
        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 10*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;
        
        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same

        
    case 4

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_04_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 1*D2R; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 1;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;  
        
        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 10*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same


    case 5

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_05_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 5*D2R; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 1;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;  
        
        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 10*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same

      case 6

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_06_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 10*D2R; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 1;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;
        conditions.conditions_TRIM_lat.daT = 10*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same

         case 7

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_07_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 0; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 0;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0; 
        
        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 20*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same
        
    case 8

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_08_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 1*D2R; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 1;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;  
        
        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 20*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value
        
        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same

    case 9

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_09_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 5*D2R; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 1;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;  
        
        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 20*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same

      case 10

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_10_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 10*D2R; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 1;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0; 
        
        % Assumptions for solving Trim-Lateral Directional equations with
        % daT and daR given
        conditions.conditions_TRIM_lat.daT = 20*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

      case 11

         % Geometry

        % Plots prefix
        prefix = 'CASO_ID_11_Cessna208_';
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
        % CAse 5: Rolling moment engine
        conditions.conditions_TRIM_lat.Rolling_maneuver = 5;
        conditions.conditions_TRIM_lat.change_dir_roll = 1; 
        conditions.conditions_TRIM_lat.LT_engine = 6000; % Newton-m for electric engine
        
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
        conditions.conditions_TRIM_lat.direction_phi = 1; % defines the directión of bank angle (positive wing right down)
        
        % Condition of initial Trim Lateral Directional
        conditions.conditions_TRIM_lat.phi_0 = 10*D2R; % defines the directión of bank angle (positive wing right down)
        % For steady State & % nhat = 1; % For level turning flight
        nhat = 1;
        conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;
        conditions.conditions_TRIM_lat.daT = 20*D2R; % Fixed value
        conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

        N_moment = 0; % Worse case scenarios % Actualizado
        conditions.conditions_TRIM_lat.N_moment = N_moment;

        %%Stability
        % OUTPUT_read_XLSX.Stability_flags.SM_des = 0.25;

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot(8) = 0; % stops saving images for longitudinal trim since they are the same

        % Medium Aircraft (e.g., Cessna Caravan, Beechcraft King Air): 3:1 to 4:1
        Ga = 0.4*m2ft; % rad/ft -> rad/m, Stick to aileron gearing ratio
        Gr = 1.5*m2ft; % rad/ft -> rad/m,
        conditions.conditions_TRIM_lat.Ga = Ga*0.5; %

end
