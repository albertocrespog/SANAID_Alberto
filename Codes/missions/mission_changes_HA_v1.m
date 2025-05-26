function [conditions OUTPUT_read_XLSX] = mission_changes_HA_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
             
%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
% Folder where to stores de Figures

m2ft = conv_UNITS.m2ft;

% Identifies folder where to store Figures
filename = '../Results/26_HA/Figs';
OUTPUT_read_XLSX.PLOT_flags.fname = filename;

%% Selects engine on or off
% Posicion_Palanca = 0; % Engien off
Posicion_Palanca = 0; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

switch mission_actual
    case 1 % ID_01
        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_01_HA_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;
        
        % Performance
                %% DESCENT

        %% descent SEGMENT ALONE

        %% Select the Number of Segments
        % prompt4 = sprintf('ENTER MISSION SEGMENTS SELECTION:\n1 - TAXY\n2 - TAKE-OFF\n3 - CLIMB\n4 - VTOL CLIMB\n5 - CRUISE\n6 - LOAD DEPLOYMENT\n7 - TURN\n8 - DESCENT\n9 - VTOL DESCENT\n10 - ALTERNATIVE AIRPORT CLIMB TO 3000FT\n11 - TURN LOITTER 45MIN\n12 - LANDING\n13 - DUMMY\n\nNOTE: Introduce the mission selection as a succession of numbers separated by a blank without any commas or brackets');
        % dlgtitle = 'Mission Segments Selection';
        % dims = [1 70];
        % answer_SEGMENTS = inputdlg(prompt4,dlgtitle,dims); % 3 5 8
        % SEGMENTS = str2num(answer_SEGMENTS{1});
        SEGMENTS = 8;
        SEGMENTS_STYPE = zeros(1,length(SEGMENTS));
        SEGMENTS_STYPE_InpDat = zeros(length(SEGMENTS),3);
        OUTPUT_read_XLSX = initialilize_segments_performance_mission(SEGMENTS,OUTPUT_read_XLSX);

        %% descent OPTIONS

        % prompt5 = sprintf('ENTER DESCENT SEGMENT SUBTYPE SELECTION:\n1 - Giving M and gamma\n2 - Giving EAS and gamma\n3 - Giving TAS and gamma\n4 - Giving M and throttle\n5 - Giving EAS and throttle\n6 - Giving TAS and throttle\n7 - Giving Vi, Vf and gamma\n8 - Minimum gamma\n9 - Slowest sink\n10 - Giving Vi, Vf and throttle');
        % dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
        % dims = [1 70];
        % answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
        % SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
        SEGMENTS_STYPE = 4;
        V_cr           = 12; % m/s
        delta_T_d = 0.01;
        h_inicial_d = 10000/m2ft;
        h_final_d = 0;
        [Temp_init,rho_init,p_init,a_init]=atmos_inter_mio(h_inicial_d);
        Mach_d = V_cr/a_init;

        OUTPUT_read_XLSX.IPP_flags.Mach_d         = Mach_d;
        OUTPUT_read_XLSX.IPP_flags.delta_T_d      = delta_T_d;
        OUTPUT_read_XLSX.IPP_flags.h_inicial_d    = h_inicial_d;
        OUTPUT_read_XLSX.IPP_flags.h_final_d      = h_final_d;

        % Saving Mission Configuration
        OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
        OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
        OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

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
