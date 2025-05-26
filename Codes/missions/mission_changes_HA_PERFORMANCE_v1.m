function [conditions OUTPUT_read_XLSX] = mission_changes_HA_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
m2ft = conv_UNITS.m2ft;
ft2m = conv_UNITS.ft2m;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;

% Identifies folder where to store Figures
filename = '../Results/26_HA/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

%% Initialize Stability Studies
% --------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All values here presented need to be defined at lest with this
% values which iomply that no studies are conducted, it is
% recommeded copy and past this at the beginning of the file and
% only modify in the sections if necessary to study variation of
% stability studies

%% Stability Studies
%% Asymmetric weight in each wing
% right wing
conditions.conditions_TRIM_lat.sign1 = 1; % Positive for right wing
conditions.conditions_TRIM_lat.L_arm1 = 0; % y location of mass in the wingspan
conditions.conditions_TRIM_lat.Z_arm1 = 0; % z location of mass in the wingspan
conditions.conditions_TRIM_lat.L_weight1 = 0; % Weight of mass in the wingspan
% left wing
conditions.conditions_TRIM_lat.sign2 = -1; % Negative for left wing
conditions.conditions_TRIM_lat.L_arm2 = 0; % y location of mass in the wingspan - opposite sing that right wing
conditions.conditions_TRIM_lat.Z_arm2 = 0; % z location of mass in the wingspan
conditions.conditions_TRIM_lat.L_weight2 = 0; % Weight of mass in the wingspan

%% Rolling Moments definition
% CAse 1: NO rolling movement external
% CAse 2: Rolling moment asymetric weight 1 different in each wing
% CAse 3: Rolling moment asymetric weight 2 - only in 1 wing 1
% CAse 4: Rolling moment asymetric weight 2 - only in 1 wing 2
conditions.conditions_TRIM_lat.Rolling_maneuver = 1;
conditions.conditions_TRIM_lat.change_dir_roll = 1;

%% Yawing Moments definition
% CAse 1: NO Yawing moment external
% CAse 2: Yawing movement drag asymetry different in each wing
% CAse 3: Yawing movement drag asymetry  - only in wing 1
% CAse 4: Yawing movement drag asymetry  - only in wing 2
% CAse 5: Yawing movement engine asymetry
% CAse 6: Yawing movement engine asymetry and drag asymmetry
% change yawing moment direction
conditions.conditions_TRIM_lat.Yawing_maneuver = 1;
conditions.conditions_TRIM_lat.change_dir_yaw = 1;
conditions.conditions_TRIM_lat.direction_beta = 1; % defines the directión of sideslip
% Additional Yawing moment
N_moment = 0; % Worse case scenarios % Actualizado
conditions.conditions_TRIM_lat.N_moment = N_moment;

%% Condition of initial Trim Lateral Directional
conditions.conditions_TRIM_lat.phi_0 = 0;
% defines the directión of bank angle (positive wing right down) For steady State & nhat = 1; For level turning flight
nhat = 0;
conditions.conditions_TRIM_lat.nhat = nhat; % nhat = 0;

%% From Lateral Trim conditions with Trim Tabs
% Assumptions for solving Trim-Lateral Directional equations with
% aileron trim-tab (daT) and rudder trim tab (daR) given by a
% constant value
conditions.conditions_TRIM_lat.daT = 0*D2R; % Fixed value
conditions.conditions_TRIM_lat.drT = 0*D2R; % Fixed value

% Assumptions for solving Trim-Lateral Directional equations with
% given Forces on aileron stick (Fa) and Rudder Peddals (in Newtons)
conditions.conditions_TRIM_lat.Fa = 0; % Fixed value
conditions.conditions_TRIM_lat.Fr = 0; % Fixed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------------------

%% Selects engine on or off
% Posicion_Palanca = 0; % Engine OFF
% Posicion_Palanca = 1; % Engine ON
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

% % Construct SEGMENT
% prompt4 = sprintf('ENTER MISSION SEGMENTS SELECTION:
% \n1 - TAXY
% \n2 - TAKE-OFF
% \n3 - CLIMB
% \n4 - VTOL CLIMB
% \n5 - CRUISE
% \n6 - LOAD DEPLOYMENT
% \n7 - TURN
% \n8 - DESCENT
% \n9 - VTOL DESCENT
% \n10 - ALTERNATIVE AIRPORT CLIMB TO 3000FT
% \n11 - TURN LOITTER 45MIN\n12 - LANDING
% \n13 - DUMMY

switch mission_actual
    case 1 % ID_01
        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_01_26_HA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Performance Analysis
        OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY = 1;
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP = 1; % Single Performance Study
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = 0; % Variable Study

        %% Construct SEGMENT
        %% CRUISE
        SEGMENTS = 5;
        SEGMENTS_STYPE = zeros(1,length(SEGMENTS));
        SEGMENTS_STYPE_InpDat = zeros(length(SEGMENTS),3);
        OUTPUT_read_XLSX = initialilize_segments_performance_mission(SEGMENTS,OUTPUT_read_XLSX);

        %% CRUISE OPTIONS
        % prompt5 = sprintf('ENTER CRUISE SEGMENT SUBTYPE SELECTION:
        % \n1 - Giving TAS and distance
        % \n2 - Giving CL and distance
        % \n3 - Giving Vi, Vf and throttle
        % \n4 - Giving TAS and CD(M)
        % \n5 - Max range giving Wf
        % \n6 - Max endurance giving Wf
        % \n7 - Giving TAS and m_bat');
        % dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
        % dims = [1 70];
        % answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
        % SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
        SEGMENTS_STYPE = 1;

        % OUTPUT_read_XLSX.STUDY_flags.variable_speed_weight_AP = 1;

        V_cr           = 35; % m/s
        mbat           = 10; % Kg
        h_inicial_cr   = 1000; % m
        dist_final_cr = 20000; % m
        OUTPUT_read_XLSX.IPP_flags.V_cr           = V_cr;
        OUTPUT_read_XLSX.IPP_flags.mbat           = mbat;
        OUTPUT_read_XLSX.IPP_flags.h_inicial_cr   = h_inicial_cr;
      
        OUTPUT_read_XLSX.IPP_flags.h_inicial_cr = h_inicial_cr;
        OUTPUT_read_XLSX.IPP_flags.dist_final_cr = dist_final_cr;
        % OUTPUT_read_XLSX.IPP_flags.V_cr = V_cr;
        % OUTPUT_read_XLSX.IPP_flags.delta_T_cr = delta_T_cr;
        % OUTPUT_read_XLSX.IPP_flags.V_ini_cr = V_ini_cr;
        % OUTPUT_read_XLSX.IPP_flags.V_fin_cr = V_fin_cr;
        % OUTPUT_read_XLSX.IPP_flags.fuel_cr = fuel_cr;
        % OUTPUT_read_XLSX.IPP_flags.Cd0_cr = Cd0_cr;
        % OUTPUT_read_XLSX.IPP_flags.k1_cr = k1_cr;
        % OUTPUT_read_XLSX.IPP_flags.k2_cr = k2_cr;
        % OUTPUT_read_XLSX.IPP_flags.mbat = mbat;

        %% Saving Mission Configuration
        OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
        OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
        OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

        % OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
        % OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
        % OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

    case 2 % ID_01
        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_02_26_HA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Performance Analysis
        OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY = 1;
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP = 0; % Single Performance Study
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = 1; % Variable Study

        %% Construct SEGMENT
        %% CRUISE
        SEGMENTS = 5;
        SEGMENTS_STYPE = zeros(1,length(SEGMENTS));
        SEGMENTS_STYPE_InpDat = zeros(length(SEGMENTS),3);
        OUTPUT_read_XLSX = initialilize_segments_performance_mission(SEGMENTS,OUTPUT_read_XLSX);

        %% CRUISE OPTIONS
        % prompt5 = sprintf('ENTER CRUISE SEGMENT SUBTYPE SELECTION:
        % \n1 - Giving TAS and distance
        % \n2 - Giving CL and distance
        % \n3 - Giving Vi, Vf and throttle
        % \n4 - Giving TAS and CD(M)
        % \n5 - Max range giving Wf
        % \n6 - Max endurance giving Wf
        % \n7 - Giving TAS and m_bat');
        % dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
        % dims = [1 70];
        % answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
        % SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
        SEGMENTS_STYPE = 7;

        OUTPUT_read_XLSX.STUDY_flags.variable_speed_weight_AP = 1;

        V_cr           = 35; % m/s
        mbat           = 10; % Kg
        h_inicial_cr   = 1000; % m
        OUTPUT_read_XLSX.IPP_flags.V_cr           = V_cr;
        OUTPUT_read_XLSX.IPP_flags.mbat           = mbat;
        OUTPUT_read_XLSX.IPP_flags.h_inicial_cr   = h_inicial_cr;

        % Definition of number of points
        N_V_VAR_perfo = 15; %
        N_m_VAR_perfo = 15; %
        OUTPUT_read_XLSX.Performance_pre_flags.N_V_VAR_perfo = N_V_VAR_perfo; %
        OUTPUT_read_XLSX.Performance_pre_flags.N_m_VAR_perfo = N_m_VAR_perfo; %

        % Calculation of Mass available
        mbat_low   = 1;
        mbat_high = 10;
        m_bat_VAR = linspace(mbat_low,mbat_high,N_m_VAR_perfo);

        % Stores Variables mass for sensitivity studies
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_low = mbat_low;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_high = mbat_high;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.m_bat_VAR = m_bat_VAR;

        %% Variation of speeds
        V_low = 15;
        V_high = 40;
        V_VAR = linspace(V_low,V_high,N_V_VAR_perfo);

        h_inicial_cr   = 5000*ft2m; % m

        % Stores Variables mass for sensitivity studies
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.V_low = V_low;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.V_high = V_high;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.V_VAR = V_VAR;

        OUTPUT_read_XLSX.IPP_flags.h_inicial_cr = h_inicial_cr;
        % OUTPUT_read_XLSX.IPP_flags.dist_final_cr = dist_final_cr;
        % OUTPUT_read_XLSX.IPP_flags.V_cr = V_cr;
        % OUTPUT_read_XLSX.IPP_flags.delta_T_cr = delta_T_cr;
        % OUTPUT_read_XLSX.IPP_flags.V_ini_cr = V_ini_cr;
        % OUTPUT_read_XLSX.IPP_flags.V_fin_cr = V_fin_cr;
        % OUTPUT_read_XLSX.IPP_flags.fuel_cr = fuel_cr;
        % OUTPUT_read_XLSX.IPP_flags.Cd0_cr = Cd0_cr;
        % OUTPUT_read_XLSX.IPP_flags.k1_cr = k1_cr;
        % OUTPUT_read_XLSX.IPP_flags.k2_cr = k2_cr;
        % OUTPUT_read_XLSX.IPP_flags.mbat = mbat;



        %% Saving Mission Configuration
        OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
        OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
        OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

        % OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
        % OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
        % OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

end
