function [conditions OUTPUT_read_XLSX] = mission_changes_PEPINO_PERFORMANCE_v1(OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications)

%Loads geometry changes
dummy = 1;
conditions.dummy = dummy;
conditions.m_TOW = OUTPUT_read_XLSX.Weights_flags.MTOW_true;
m2ft = conv_UNITS.m2ft;
ft2m = conv_UNITS.ft2m;
D2R = conv_UNITS.D2R;
R2D = conv_UNITS.R2D;

%% Initialize Stability Studies
conditions = initialization_STABILITY_Analysis(conditions,conv_UNITS);

%% Selects engine on or off
% Posicion_Palanca = 0; % Engine OFF
% Posicion_Palanca = 1; % Engine ON
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

%% Construct SEGMENT
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

%% Initialization PLOTS in Missions
OUTPUT_read_XLSX = initialization_PLOTS_in_MISSIONS(OUTPUT_read_XLSX);
% plot_perfo1 = Prints PLOTS PERFORMANCE STUDY
% plot_perfo2 = Prints plots of Performance for Variable V and mass - Electric
% plot_perfo3 = Prints plots of Performance for Variable h and V
% plot_perfo4 = Prints plots of Performance Glide for Variable h and V
% plot_perfo5 = Prints plots of Performance Glide for Variable h Max

% Variable Mass and Speed Performance Studies - Electric - Perfo1
% Variable Mass and Speed Performance Studies - Fuel - Perfo2
% Glid Performance Variable Speed - Altitude - Perfo3
% Glid Performance Variable Speed - Altitude - Glide - Perfo4
% Glid Performance Variable Speed - Altitude - Glide Max - Perfo5
% Variable Speed Performance Studies - Perfo6
% Variable Mass Performance Studies - Perfo7

%% Selects engine on or off
% Posicion_Palanca = 0; % Engine OFF
% Posicion_Palanca = 1; % Engine ON
Posicion_Palanca = 1; % Engien on
conditions.Posicion_Palanca = Posicion_Palanca;

switch mission_actual
    case 1

        %% Plots prefix
        prefix = 'CASEID_01_AC_03_PEPINO_PERFORMANCE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Geometry

        %% Porpulsion
        Thrust = 176*1.5;
        OUTPUT_read_XLSX.Propulsive_flags.propul(3) = Thrust; % EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine

        %% Performance Analysis
        OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY = 1; % Activates Performance
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP = 1; % Single Performance Study
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = 0; % Variable Study

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot_perfo1 = 1;

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

        OUTPUT_read_XLSX.STUDY_flags.variable_speed_weight_AP = 1;

        V_cr           = 200; % m/s
        h_inicial_cr   = 10000*ft2m; % m
        dist_final_cr = 300*1000;
        OUTPUT_read_XLSX.IPP_flags.V_cr           = V_cr;
        OUTPUT_read_XLSX.IPP_flags.h_inicial_cr   = h_inicial_cr;
        OUTPUT_read_XLSX.IPP_flags.dist_final_cr = dist_final_cr;


        %% Saving Mission Configuration
        OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
        OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
        OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

    case 2

        % Plots prefix
        prefix = 'CASEID_02_AC_03_PEPINO_PERFORMANCE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Geometry

        %% Porpulsion
        Thrust = 176*1.5;
        OUTPUT_read_XLSX.Propulsive_flags.propul(3) = Thrust; % EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine

        %% Performance Analysis
        OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY = 1; % Activates Performance
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP = 0; % Single Performance Study
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = 1; % Variable Study
        OUTPUT_read_XLSX.STUDY_flags.Perfo3 = 1; % Variable Study
        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot_perfo3 = 1;


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

        dist_final_cr = 300*1000;
        OUTPUT_read_XLSX.IPP_flags.dist_final_cr = dist_final_cr;

        % Stores Variables altitude for sensitivity studies
        h_low = 5000*ft2m;
        h_high = 10000*ft2m;

        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_low = h_low;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_high = h_high;

        % Stores Variables altitude for sensitivity studies
        V_low = 100;
        V_high = 240;

        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.V_low = V_low;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.V_high = V_high;

        % Number of Vector elements for both altitude and velocity
        % sensitivity study
        N_h_VAR_perfo = 3 ;
        N_V_VAR_perfo = 10;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.N_h_VAR_perfo = N_h_VAR_perfo; %
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.N_V_VAR_perfo = N_V_VAR_perfo; %
        
        % OUTPUT_read_XLSX.IPP_flags.h_inicial_cr = h_inicial_cr;
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

    case 3 % ID_01
        
        % Plots prefix
        prefix = 'CASEID_03_AC_03_PEPINO_PERFORMANCE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        % Geometry

        % Porpulsion
        Thrust = 176*0.5;
        OUTPUT_read_XLSX.Propulsive_flags.propul(3) = Thrust; % EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine

       
         %% Performance Analysis
        OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY = 1; % Activates Performance
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP = 1; % Single Performance Study
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = 0; % Variable Study

        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot_perfo1 = 1;
        
        %% Construct SEGMENT
        %% CRUISE
        SEGMENTS = 8;
        SEGMENTS_STYPE = zeros(1,length(SEGMENTS));
        SEGMENTS_STYPE_InpDat = zeros(length(SEGMENTS),3);
        OUTPUT_read_XLSX = initialilize_segments_performance_mission(SEGMENTS,OUTPUT_read_XLSX);

        %% DESCENT
        %% descent OPTIONS
        % prompt5 = sprintf('ENTER DESCENT SEGMENT SUBTYPE SELECTION:
        % \n1 - Giving M and gamma
        % \n2 - Giving EAS and gamma
        % \n3 - Giving TAS and gamma
        % \n4 - Giving M and throttle
        % \n5 - Giving EAS and throttle
        % \n6 - Giving TAS and throttle
        % \n7 - Giving Vi, Vf and gamma
        % \n8 - Minimum gamma
        % \n9 - Slowest sink
        % \n10 - Giving Vi, Vf and throttle
        % dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
        % dims = [1 70];
        % answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
        % SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});

        SEGMENTS_STYPE = 6;

        TAS_d = 200; % (m/s)
        delta_T_d = 0.0;
        h_inicial_d = 10000/m2ft;
        h_final_d = 0;
        [Temp_init,rho_init,p_init,a_init]=atmos_inter_mio(h_inicial_d);
        Mach_d = TAS_d/a_init;

        OUTPUT_read_XLSX.IPP_flags.Mach_d         = Mach_d;
        OUTPUT_read_XLSX.IPP_flags.delta_T_d      = delta_T_d;
        OUTPUT_read_XLSX.IPP_flags.h_inicial_d    = h_inicial_d;
        OUTPUT_read_XLSX.IPP_flags.h_final_d      = h_final_d;
        OUTPUT_read_XLSX.IPP_flags.TAS_d         = TAS_d;
        OUTPUT_read_XLSX.IPP_flags.h_final_d      = h_final_d;

        % gamma_d = OUTPUT_read_XLSX.IPP_flags.gamma_d = gamma_d;
        % EAS_d = OUTPUT_read_XLSX.IPP_flags.EAS_d = EAS_d;
        % TAS_d = OUTPUT_read_XLSX.IPP_flags.TAS_d = TAS_d;
        % V_ini_d = OUTPUT_read_XLSX.IPP_flags.V_ini_d = V_ini_d;
        % V_fin_d = OUTPUT_read_XLSX.IPP_flags.V_fin_d = V_fin_d;

        % Saving Mission Configuration
        OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
        OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
        OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;


    case 4 % ID_01
        % Geometry

        % Plots prefix
        prefix = 'CASEID_04_AC_03_PEPINO_PERFORMANCE_';
        OUTPUT_read_XLSX.PLOT_flags.prefix = prefix;

        %% Performance Analysis
        OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY = 1; % Activates Performance
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP = 0; % Single Performance Study
        OUTPUT_read_XLSX.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var = 1; % Variable Study
        OUTPUT_read_XLSX.STUDY_flags.Perfo3 = 1; % Variable Study
        %% PLOTS
        OUTPUT_read_XLSX.PLOT_flags.plot_perfo3 = 1;

        %% Construct SEGMENT
        %% DESCENT
        SEGMENTS = 8;
        SEGMENTS_STYPE = zeros(1,length(SEGMENTS));
        SEGMENTS_STYPE_InpDat = zeros(length(SEGMENTS),3);
        OUTPUT_read_XLSX = initialilize_segments_performance_mission(SEGMENTS,OUTPUT_read_XLSX);

        %% DESCENT
        %% descent OPTIONS
        % prompt5 = sprintf('ENTER DESCENT SEGMENT SUBTYPE SELECTION:
        % \n1 - Giving M and gamma
        % \n2 - Giving EAS and gamma
        % \n3 - Giving TAS and gamma
        % \n4 - Giving M and throttle
        % \n5 - Giving EAS and throttle
        % \n6 - Giving TAS and throttle
        % \n7 - Giving Vi, Vf and gamma
        % \n8 - Minimum gamma
        % \n9 - Slowest sink
        % \n10 - Giving Vi, Vf and throttle
        % dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
        % dims = [1 70];
        % answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
        % SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});

        % Stores Variables altitude for sensitivity studies
        h_low = 10000*ft2m;
        h_high = 30000*ft2m;
        h_terminal = 0*ft2m;

        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_low = h_low;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_high = h_high;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.h_terminal = h_terminal;

        % Stores Variables altitude for sensitivity studies
        V_low = 180;
        V_high = 240;

        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.V_low = V_low;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.V_high = V_high;

        % Number of Vector elements for both altitude and velocity
        % sensitivity study
        N_h_VAR_perfo = 5;
        N_V_VAR_perfo = 6;
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.N_h_VAR_perfo = N_h_VAR_perfo; %
        OUTPUT_read_XLSX.IPP_flags.Variable_Study5.N_V_VAR_perfo = N_V_VAR_perfo; %

        SEGMENTS_STYPE = 4;

        % Saving Mission Configuration
        OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
        OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
        OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;

end
