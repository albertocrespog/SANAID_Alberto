function [conditions OUTPUT_read_XLSX] = mission_changes_SOLARTII_v3(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)  
           
%Loads a de geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

%% Performance for all Missiones

% Construct SEGMENT
%% CRUISE
% prompt4 = sprintf('ENTER MISSION SEGMENTS SELECTION:\n1 - TAXY\n2 - TAKE-OFF\n3 - CLIMB\n4 - VTOL CLIMB\n5 - CRUISE\n6 - LOAD DEPLOYMENT\n7 - TURN\n8 - DESCENT\n9 - VTOL DESCENT\n10 - ALTERNATIVE AIRPORT CLIMB TO 3000FT\n11 - TURN LOITTER 45MIN\n12 - LANDING\n13 - DUMMY\n\nNOTE: Introduce the mission selection as a succession of numbers separated by a blank without any commas or brackets');
% dlgtitle = 'Mission Segments Selection';
% dims = [1 70];
% answer_SEGMENTS = inputdlg(prompt4,dlgtitle,dims); % 3 5 8
% SEGMENTS = str2num(answer_SEGMENTS{1});
SEGMENTS = 5;
SEGMENTS_STYPE = zeros(1,length(SEGMENTS));
SEGMENTS_STYPE_InpDat = zeros(length(SEGMENTS),3);
OUTPUT_read_XLSX = initialilize_segments_performance_mission(SEGMENTS,OUTPUT_read_XLSX);

%% CRUISE OPTIONS
% prompt5 = sprintf('ENTER CRUISE SEGMENT SUBTYPE SELECTION:\n1 - Giving TAS and distance\n2 - Giving CL and distance\n3 - Giving Vi, Vf and throttle\n4 - Giving TAS and CD(M)\n5 - Max range giving Wf\n6 - Max endurance giving Wf\n7 - Giving TAS and m_bat');
% dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
% dims = [1 70];
% answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
% SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
SEGMENTS_STYPE = 7;
%% CRUISE OPTIONS - 7 Giving TAS and m_bat
% prompt6 = sprintf('ENTER CRUISE INPUT DATA:\n1 - V_cr [m/s]\n2 - m_bat [kg]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
% dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
% dims = [1 70];
% answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
% SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});

V_cr           = 12; % m/s
mbat           = 10; % Kg
h_inicial_cr   = 1000; % m
OUTPUT_read_XLSX.IPP_flags.V_cr           = V_cr;
OUTPUT_read_XLSX.IPP_flags.mbat           = mbat;
OUTPUT_read_XLSX.IPP_flags.h_inicial_cr   = h_inicial_cr;

% Definition of number of points
N_V_VAR_perfo = OUTPUT_read_XLSX.Performance_pre_flags.N_V_VAR_perfo; %
N_m_VAR_perfo = OUTPUT_read_XLSX.Performance_pre_flags.N_m_VAR_perfo; %


% Calculation of Mass available
mbat_low   = 2;
mbat_high = 12;
m_bat_VAR = linspace(mbat_low,mbat_high,N_m_VAR_perfo);

% Stores Variables mass for sensitivity studies
OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_low = mbat_low;
OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_high = mbat_high;
OUTPUT_read_XLSX.IPP_flags.Variable_Study1.m_bat_VAR = m_bat_VAR;

% m_empty = Weight_tier.m_empty; % Peso en vacío
% m_payload = Weight_tier.m_payload; % Carga de pago
% m_VAR_min = m_empty + m_payload + mbat_low; % total mass
% m_VAR_max = m_empty + m_payload + mbat_high; % total mass
% 
% % Variation of speeds
% V_low = V_stall_min_low;
% V_high = 1.25*V_min_ope_high;
% V_VAR = linspace(V_low,V_high,N_V_VAR_perfo);
% 
% % Stores Variables mass for sensitivity studies
% OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_low = V_low;
% OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_high = V_high;
% OUTPUT_read_XLSX.IPP_flags.Variable_Study1.m_bat_VAR = V_VAR;


% Saving Mission Configuration
OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

% case 4 % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg
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

% Centro de Gravedad
OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = 0.7424;
OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG_ac = 0.691;
OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = 0.0;

% Inercias
OUTPUT_read_XLSX.Weights_flags.I_xx = 6.64225; %
OUTPUT_read_XLSX.Weights_flags.I_yy = 0.96730; %
OUTPUT_read_XLSX.Weights_flags.I_zz = 7.60421; %
OUTPUT_read_XLSX.Weights_flags.I_xz = -0.01422; %

% Aerodynamics
OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1 = 4; % Selection of TXT that are used for the aerodynamic analysis (w1)
OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP = 11; % Selection of TXT that are used for the aerodynamic analysis (w1)

% Performance
OUTPUT_read_XLSX.IPP_flags.V_cr = OUTPUT_read_XLSX.IPP_flags.V_cr;

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
    case 1  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop01_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 14*2.54/100;

        D_prop = OUTPUT_read_XLSX.Propulsive_flags.D_prop;

    case 2  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop02_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 16*2.54/100;

    case 3  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop03_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 18*2.54/100;
    
    case 4  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop04_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 20*2.54/100;

    case 5  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop05_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 22*2.54/100;

    case 6  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop06_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 24*2.54/100;

    case 7  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop07_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 26*2.54/100;
    
    
    case 8  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop08_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 28*2.54/100;

    case 9  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop09_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 30*2.54/100;

    case 10  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop10_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 32*2.54/100;

    case 11  % ID_04 - Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_Dprop11_SOLARTII_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 34*2.54/100;

end