function [conditions OUTPUT_read_XLSX] = mission_changes_MILVUS2_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
m2ft = conv_UNITS.m2ft;
ft2m = conv_UNITS.ft2m;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;

% % Identifies folder where to sAirtore Figures
% filename = '../Results/29_MILVUS/Figs';
% OUTPUT_read_XLSX.PLOT_flags.fname = filename;

%% Initialize Stability Studies
conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS);

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
        
        %% Aerodynamic
        fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_01_29_MILVUS2_PERFORMANCE_';
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
        SEGMENTS_STYPE = 7;


        % Type of battery
        % SE_LiFePO4
        % SE_LiPo
        % SE_FuelCells
       
        % Battery
        nominal_voltage = 22.2;
        total_capacity = 22000/1000; % Ah
        mass_battery = 1.963;
        SE_battery = total_capacity*nominal_voltage/mass_battery;
        OUTPUT_read_XLSX.Propulsive_flags.type_battery = 2; %
        OUTPUT_read_XLSX.Propulsive_flags.SE_LiPo = SE_battery; %
        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 17*2.54/100;

        OUTPUT_read_XLSX.STUDY_flags.variable_speed_weight_AP = 1;

        V_cr           = 32; % m/s
        mbat_max           = 3.872; % Kg
        mbat_max = OUTPUT_read_XLSX.Weights_flags.MF_true;
        h_inicial_cr   = 1200; % m
        OUTPUT_read_XLSX.IPP_flags.V_cr           = V_cr;
        OUTPUT_read_XLSX.IPP_flags.mbat           = mbat_max;
        OUTPUT_read_XLSX.IPP_flags.h_inicial_cr   = h_inicial_cr;

        % Definition of number of points
        N_V_VAR_perfo = 15; %
        N_m_VAR_perfo = 15; %
        OUTPUT_read_XLSX.Performance_pre_flags.N_V_VAR_perfo = N_V_VAR_perfo; %
        OUTPUT_read_XLSX.Performance_pre_flags.N_m_VAR_perfo = N_m_VAR_perfo; %

        % Calculation of Mass available
        mbat_low   = 1;
        mbat_high = mbat_max;
        m_bat_VAR = linspace(mbat_low,mbat_high,N_m_VAR_perfo);

        % Stores Variables mass for sensitivity studies
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_low = mbat_low;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_high = mbat_high;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.m_bat_VAR = m_bat_VAR;

        %% Variation of speeds
        V_low = 20;
        V_high = 40;
        V_VAR = linspace(V_low,V_high,N_V_VAR_perfo);

        h_inicial_cr   = 1200; % m

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

        case 2 % ID_01
        
        %% Aerodynamic
        fuse_aero_FLOW_and_CBM = 2;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.fuse_aero_FLOW_and_CBM = fuse_aero_FLOW_and_CBM; %
        
        polar_model = 3;
        OUTPUT_read_XLSX.Aerodynamic_Data_flags.polar_model = polar_model;   

        % Geometry

        % Plots prefix
        prefix = 'CASO_ID_02_6_MILVUS_PERFORMANCE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Performance Analysis
        OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY = 1;
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP = 0; % Single Performance Study
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = 1; % Variable Study
        OUTPUT_read_XLSX.STUDY_flags.variable_speed_weight_AP_fuel = 1;


        %% Construct SEGMENT
        %% CRUISE
        SEGMENTS = 5;
        % SEGMENTS_STYPE = zeros(1,length(SEGMENTS));
        % SEGMENTS_STYPE_InpDat = zeros(length(SEGMENTS),3);
        % OUTPUT_read_XLSX = initialilize_segments_performance_mission(SEGMENTS,OUTPUT_read_XLSX);

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

        
        % Type of battery
        % SE_LiFePO4
        % SE_LiPo
        % SE_FuelCells
       
        % Battery
        nominal_voltage = 22.2;
        total_capacity = 22000/1000; % Ah
        mass_battery = 1.963;

                


        SE_battery = total_capacity*nominal_voltage/mass_battery;
        OUTPUT_read_XLSX.Propulsive_flags.type_battery = 2; %
        OUTPUT_read_XLSX.Propulsive_flags.SE_LiPo = SE_battery; %
        OUTPUT_read_XLSX.Propulsive_flags.D_prop = 15*2.54/100;
        
        OUTPUT_read_XLSX.STUDY_flags.variable_speed_weight_AP = 1;

        V_cr           = 32; % m/s
        mbat_max           = 3.872; % Kg
        mbat_max = OUTPUT_read_XLSX.Weights_flags.MF_true;
        h_inicial_cr   = 1200; % m
        OUTPUT_read_XLSX.IPP_flags.V_cr           = V_cr;
        OUTPUT_read_XLSX.IPP_flags.mbat           = mbat_max;
        OUTPUT_read_XLSX.IPP_flags.h_inicial_cr   = h_inicial_cr;

        % Definition of number of points
        N_V_VAR_perfo = 15; %
        N_m_VAR_perfo = 15; %
        OUTPUT_read_XLSX.Performance_pre_flags.N_V_VAR_perfo = N_V_VAR_perfo; %
        OUTPUT_read_XLSX.Performance_pre_flags.N_m_VAR_perfo = N_m_VAR_perfo; %

        % Calculation of Mass available
        mbat_low   = 1;
        mbat_high = mbat_max;
        m_bat_VAR = linspace(mbat_low,mbat_high,N_m_VAR_perfo);

        % Stores Variables mass for sensitivity studies
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_low = mbat_low;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.mbat_high = mbat_high;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study1.m_bat_VAR = m_bat_VAR;

        %% Variation of speeds
        V_low = 20;
        V_high = 40;
        V_VAR = linspace(V_low,V_high,N_V_VAR_perfo);

        h_inicial_cr   = 1200; % m

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
