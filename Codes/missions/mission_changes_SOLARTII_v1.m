function [conditions OUTPUT_read_XLSX] = mission_changes_SOLARTII_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);


%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

%Loadsa de geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

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

% Identifies folder where to store Figures
filename = '../Results/SOLARTII/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

switch mission_actual
    case 0 % ID_00 - Wing 11, sweep 14.25, taper = 0.6, dihedral = 0 deg
        % Geometry
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 14.25*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 0*(pi/180); % Dihedral of w1 (deg)

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)

        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= y_loc_1R_y1_VTP_CAD*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + 0.24; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Twin VTp geometry VTP2
    	y_loc_1R_y1_2VTP_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_2VTP_CAD = z_loc_1R_y1_VTP_CAD + 0.24; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_2VTP_CAD = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 10; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)

        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_01_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        

    case 1 % ID_01 - Wing 11, sweep 14.25, taper = 0.6, dihedral = 0 deg
        % Geometry
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 14.25*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 0*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 1; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)

        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_01_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 2 % ID_01 - Wing 12, sweep 14.25, taper = 0.6, dihedral = 2.5 deg
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 14.25*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 2.5*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 2; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)

        % Performance

        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_02_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 3 % ID_03 - Wing 13, sweep 14.25, taper = 0.6, dihedral = 5 deg
        % Geometry
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 14.25*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 5*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 3; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)

        % Performance

        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_03_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 4 % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 25*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 0*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 4; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)

        % Performance
        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_04_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 5 % ID_05 - Wing 22, sweep 25, taper = 0.6, dihedral = 2.5 deg
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 25*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 2.5*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 5; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)
        % Performance

        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_05_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 6 % ID_06 - Wing 23, sweep 25, taper = 0.6, dihedral = 5 deg
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 25*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 5*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 6; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)
        % Performance

        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_06_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 7 % ID_07 - Wing 31, sweep 35, taper = 0.6, dihedral = 0 deg
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 35*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 0*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 7; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)
        % Performance

        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_07_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 8 % ID_08 - Wing 32, sweep 35, taper = 0.6, dihedral = 2.5 deg
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 35*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 2.5*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 8; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)
        % Performance

        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_08_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
    case 9 % ID_09 - Wing 33, sweep 35, taper = 0.6, dihedral = 5 deg
        % Wing geometry
        y_loc_1R_y1_w1_CA = 0.5; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        y_loc_1R_y2_w1_CAD = 3; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_w1_CAD = 0.315; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        Lambda_LE_w1_e = 35*(pi/180); % Sweep of w1 (deg)
        dihedral_w1_e = 5*(pi/180); % Dihedral of w1 (deg)
        b_VTP = 0.24; % Wingspan of VTP

        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CA = y_loc_1R_y1_w1_CA; % y loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y2_w1_CAD = y_loc_1R_y2_w1_CAD; % y loc of wing (w1) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_w1_CAD = x_loc_1R_y1_w1_CAD; % x loc of wing (w1) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e = Lambda_LE_w1_e; % Sweep of w1 (deg)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e = dihedral_w1_e; % Dihedral of w1 (deg)
        
        % Twin VTp geometry VTP1
    	y_loc_1R_y1_VTP_CAD	= y_loc_1R_y2_w1_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        z_loc_1R_y1_VTP_CAD	= (y_loc_1R_y1_VTP_CAD-y_loc_1R_y1_w1_CA)*tan(dihedral_w1_e); % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	z_loc_1R_y2_VTP_CAD	= z_loc_1R_y1_VTP_CAD + b_VTP; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	x_loc_1R_y1_VTP_CAD	 = (y_loc_1R_y2_w1_CAD - y_loc_1R_y1_w1_CA)*tan(Lambda_LE_w1_e) + x_loc_1R_y1_w1_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

    	OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_2VTP_CAD = y_loc_1R_y1_VTP_CAD; % relative y loc of  (VTP) root chord LE position (distance from CAD refference point)
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y1_2VTP_CAD = z_loc_1R_y1_VTP_CAD; % relative z loc of  (VTP) root chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.z_loc_1R_y2_2VTP_CAD = z_loc_1R_y2_VTP_CAD; % z loc of wing (VTP) tip chord LE position (distance from CAD refference point)
    	OUTPUT_read_XLSX.InputGeometry_Data_flags.x_loc_1R_y1_2VTP_CAD = x_loc_1R_y1_VTP_CAD; % x loc of wing (VTP) root chord LE position (distance from CAD refference point)

        % Aerodynamics
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 9; % Selection of TXT that are used for the aerodynamic analysis (w1)
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)
        % Performance

        OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;
        % Plots prefix
        prefix = 'CASO_ID_09_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
end