function [conditions OUTPUT_read_XLSX] = mission_changes_FR_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)
             
%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
% Folder where to stores de Figures

m2ft = conv_UNITS.m2ft;

% Identifies folder where to store Figures
filename = '../Results/27_FR/Figs';
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

V_cr           = 200; % m/s
mbat           = 1; % Kg
h_inicial_cr   = 1000; % m
OUTPUT_read_XLSX.IPP_flags.V_cr           = V_cr;
OUTPUT_read_XLSX.IPP_flags.mbat           = mbat;
OUTPUT_read_XLSX.IPP_flags.h_inicial_cr   = h_inicial_cr;

% Definition of number of points
N_V_VAR_perfo = OUTPUT_read_XLSX.Performance_pre_flags.N_V_VAR_perfo; %
N_m_VAR_perfo = OUTPUT_read_XLSX.Performance_pre_flags.N_m_VAR_perfo; %

% Calculation of Mass available
mbat_low   = 0.1;
mbat_high = 0.2;
m_bat_VAR = linspace(mbat_low,mbat_high,N_m_VAR_perfo);

% Stores Variables mass for sensitivity studies
OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_low = mbat_low;
OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_high = mbat_high;
OUTPUT_read_XLSX.IPP_flags.Variable_Study1.m_bat_VAR = m_bat_VAR;

% Saving Mission Configuration
OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

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

