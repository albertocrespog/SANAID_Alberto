function [Plot_Options,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_AERO_DATA,Storing_PERFORMANCE_DATA_1,Storing_PERFORMANCE_DATA_2,Storing_STABILITY_DATA1,Storing_STABILITY_DATA2,Storing_STABILITY_DATA2B,...
    Storing_STABILITY_DATA3,Storing_STABILITY_DATA4A,Storing_STABILITY_DATA4B,Storing_STABILITY_DATA4C,Storing_STABILITY_DATA4D,Storing_STABILITY_DATA5,Storing_PROPULSION_DATA,OUTPUT_read_XLSX_mod] = conduct_Analysis(conv_UNITS,AC_CONFIGURATION,...
    OUTPUT_read_XLSX,case_AC,Fig,Plot_Options,mission_actual,Propulsion,Modifications)

MISSIONS_STUDY = OUTPUT_read_XLSX.STUDY_flags.MISSIONS_STUDY; % Conducts Mission Studies

if MISSIONS_STUDY == 1
    [conditions OUTPUT_read_XLSX_mod] = select_mission_DATA(case_AC,OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
%     [conditions OUTPUT_read_XLSX] = select_mission_DATA(case_AC,OUTPUT_read_XLSX,mission_actual,conv_UNITS);
    Saving_data_OUTPUT_read_XLSX(conditions,OUTPUT_read_XLSX_mod);
else
    dummy = 1;
    conditions.dummy = dummy; 
    OUTPUT_read_XLSX_mod = OUTPUT_read_XLSX;
end

Saving_AC_CONFIGURATION(AC_CONFIGURATION,OUTPUT_read_XLSX_mod);

Dxcg = Modifications.Dxcg; 
% Geometric Data
% Geo_input_tier = Generation_Input_Geometric_Data_v2(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod);
Geo_input_tier = Generation_Input_Geometric_Data_v3(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod);
% Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod,case_AC);
% Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod,case_AC);
Geo_tier = Generation_Geometric_Data_v6(Geo_input_tier,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod,case_AC);
%% Generation of Geometry
[Fig,Body_Geo,meshData] = Generation_Fuselage_Data_old(Geo_tier,OUTPUT_read_XLSX_mod,Fig);% Defines Propulsion DATA 
Storing_GEO_DATA = Saving_data_Geo(Geo_tier,Body_Geo,meshData,OUTPUT_read_XLSX_mod);

%% Defines Estimation of Weights according different methods densities
Weight_tier = Generation_Weight_Data(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,OUTPUT_read_XLSX_mod);
if MISSIONS_STUDY == 1
    Weight_tier.m_TOW = conditions.m_TOW;
    Storing_WEIGHT_DATA = Saving_data_Weight(Weight_tier,OUTPUT_read_XLSX_mod);
else
    Storing_WEIGHT_DATA = Saving_data_Weight(Weight_tier,OUTPUT_read_XLSX_mod);
end

%% Generates the file for Aerodynamic Data
% Indentifies preliminary performance results
Performance_preliminar = Generate_Performance_preliminary(OUTPUT_read_XLSX_mod);
 
%% Defines the Reading files    
if OUTPUT_read_XLSX_mod.STUDY_flags.AERODYNAMIC_STUDY == 1
    [Storing_AERO_DATA] = Conduct_Aerodynamic_Study(OUTPUT_read_XLSX_mod,AC_CONFIGURATION,conv_UNITS,Performance_preliminar,...
        case_AC,Storing_WEIGHT_DATA,Storing_GEO_DATA,conditions);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_AERO_DATA = dummy;
    % Warning
    warning_areo = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1 - Code in PAUSE';
    disp(warning_areo)
    pause
end

%% Data for analysis of Prop DATA
VECTOR_Prop = Get_Comparing_vectors_Propulsion(case_AC,OUTPUT_read_XLSX_mod);

%% Generates the file for Prop Data to be used in the codes
% actualizes the Prop geometry
% Select:
% Propulsion model: SEE read_prop_files_May2020.m FOR DETAILS OF PROP MODEL!!
Flag.APC_model = 1;
Flag.Wind_tunnel_model = 1;
Flag.compare_prop_models = 1;

%% Propulsive Model
[Prop_data] = Generation_Propulsion_Data(AC_CONFIGURATION,OUTPUT_read_XLSX_mod);
Storing_PROPULSION_DATA.Prop_data = Prop_data;

% Defines Propulsion DATA
% [Fig] = plot_prop_APC(Data_P,Plot_Options,Fig,prefix,VECTOR_Prop)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Conducts Propulsion optimization %
if OUTPUT_read_XLSX_mod.STUDY_flags.PROP_STUDY == 1
    Conduct_Prop_Optimization
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if OUTPUT_read_XLSX_mod.STUDY_flags.PERFORMANCE_STUDY == 1
    [OUTPUT_read_XLSX_mod, Storing_PERFORMANCE_DATA_1, Storing_PERFORMANCE_DATA_2, Storing_PERFORMANCE_DATA_3] = Conduct_Performance_Study_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
        OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA)% Conducts Performance optimization %
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA_1 = dummy;
    Storing_PERFORMANCE_DATA_2 = dummy;
    Storing_PERFORMANCE_DATA_3 = dummy;
end


% if OUTPUT_read_XLSX_mod.STUDY_flags.PERFORMANCE_STUDY == 1
%     % Performance nominal
%     Storing_PERFORMANCE_DATA_1 = Conduct_Performance_Optimization_v2(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
%         OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA);
%     % Performance Analysis integrated with AP Codes varying Conditions
%     if OUTPUT_read_XLSX_mod.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var == 1
% 
%         if OUTPUT_read_XLSX_mod.STUDY_flags.variable_speed_weight_AP == 1
%             Storing_PERFORMANCE_DATA_2 = Conduct_Performance_Optimization_v4(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
%                 OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA);
%         else
%             % Stores dummy variable for continuity
%             dummy = 1;
%             Storing_PERFORMANCE_DATA_2 = dummy;
%         end
% 
%      % Performance Analysis integrated with AP Codes varying Conditions
%         if OUTPUT_read_XLSX_mod.STUDY_flags.variable_speed_AP == 1
%             Storing_PERFORMANCE_DATA_3 = Conduct_Performance_Optimization_v5(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
%                 OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA);
%         else
%             % Stores dummy variable for continuity
%             dummy = 1;
%             Storing_PERFORMANCE_DATA_3 = dummy;
%         end
%     end
% else
%     % Stores dummy variable for continuity
%     dummy = 1;
%     Storing_PERFORMANCE_DATA_1 = dummy;
%     Storing_PERFORMANCE_DATA_2 = dummy;
%     Storing_PERFORMANCE_DATA_3 = dummy;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CONFIGURACIÓN PUSHER. MODELO DE HELICÓPTERO%%%%%%%%%%%%%%%
Posicion_Palanca = 1;
% Will be used in future version of propulsion
%% Propulsion Generation
alpha_f = 0*conv_UNITS.D2R;
beta_f = 0*conv_UNITS.D2R;
conditions.alpha_f = alpha_f;
conditions.beta_f = beta_f;
conditions.h = Storing_AERO_DATA.Performance.h;
conditions.V = Storing_AERO_DATA.Performance.V;
conditions.study_var_xcg = 0;
conditions.x_XCG = OUTPUT_read_XLSX_mod.InputGeometry_Data_flags.x_XCG;
conditions.z_XCG = OUTPUT_read_XLSX_mod.InputGeometry_Data_flags.z_XCG;

% Selects XCG depending if it's from desired stability conditions in Forward Flight (XCG_FF = 1) or
% AXIAL Flight XCG_FF = 0)
XCG_FF = OUTPUT_read_XLSX_mod.Stability_flags.XCG_FF;

if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dynamic Stability Formulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Trim Studies
    %     if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Trim == 1
    %         conditions.m_TOW = m_TOW;
    %         % only_trim = OUTPUT_read_XLSX.Stability_flags.only_trim;
    %         % Forces only trim so that obtains preliminary resutls
    %         only_trim = 1;
    %         StabilityModel = OUTPUT_read_XLSX.Stability_flags;
    %         [TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Propulsion_values] = ...
    %             Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    %             Body_Geo,Design_criteria,Posicion_Palanca,SM_des,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,...
    %             Performance_preliminar);
    %
    % %         %% XCG range
    % %         x_XCG_fwd = TRIM_RESULTS.x_XCG_fwd;
    % %         x_XCG_rwd = TRIM_RESULTS.x_XCG_rwd;
    % %         N_x_XCG_VAR = 100;
    % %         x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
    % %         Plot_Options.x_XCG_VAR = x_XCG_VAR;
    %     end
    
%     rho = 1.225; % At sea level
%     % rho = 0.96296; % At 8000 ft
%     % rho = 0.74628; % At 16000 ft
%     conditions.rho = rho;

% % case 5 - conf_missiles = 5 weight variation with 0 missiles;
% % case 4 - conf_missiles = 4 weight variation with 1 missiles;
% % case 3 - conf_missiles = 3 weight variation with 2 missiles;
% % case 2 - conf_missiles = 2 weight variation with 3 missiles;
% % case 1 - conf_missiles = 1 weight variation with 4 missiles;
% conf_missiles = 1;
% conditions.conf_missiles = conf_missiles;

    %% Regular Stability Study
    if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Regular == 1
        only_trim = OUTPUT_read_XLSX_mod.Stability_flags.only_trim;
        Storing_STABILITY_DATA1 = Conduct_Stability_Study1(conditions,Prop_data,conv_UNITS,Posicion_Palanca,XCG_FF,...
            OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA,Storing_WEIGHT_DATA,...
            Storing_AERO_DATA,case_AC,Dxcg);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA1 = dummy;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Sensitivity Study for Velocity and Mass
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
        only_trim = OUTPUT_read_XLSX_mod.Stability_flags.only_trim;
        Storing_STABILITY_DATA2 = Conduct_Stability_Study2(conditions,Prop_data,conv_UNITS,...
            Posicion_Palanca,XCG_FF,OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,case_AC,Plot_Options,...
            Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_AERO_DATA,Modifications);
        % 
        Storing_STABILITY_DATA2B = Conduct_Stability_Study2B(conditions,Prop_data,conv_UNITS,...
            Posicion_Palanca,XCG_FF,OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,case_AC,Plot_Options,...
            Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_AERO_DATA,Modifications);

        Plot_Options = Storing_STABILITY_DATA2.Plot_Options;
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA2 = dummy;
        Storing_STABILITY_DATA2B = dummy;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Sensitivity for varying Xcg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Study of variation of Velocity with XCG selected for trim conditions at desired Static Margin
    if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_var_XCG
        only_trim = OUTPUT_read_XLSX_mod.Stability_flags.only_trim;
        Storing_STABILITY_DATA3 = Conduct_Stability_Study3(conditions,Prop_data,conv_UNITS,...
    Posicion_Palanca,XCG_FF,OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,case_AC,Storing_STABILITY_DATA1,Plot_Options,...
    Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_AERO_DATA);
        Plot_Options = Storing_STABILITY_DATA3.Plot_Options;
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA3 = dummy;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Lateral Directional Trim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_lat == 1
        Storing_STABILITY_DATA4A = Conduct_Stability_Study4A(conditions,conv_UNITS,Storing_STABILITY_DATA1,OUTPUT_read_XLSX_mod,...
            Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,case_AC);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA4A = dummy;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Lateral Directional Trim - Asymmetries
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_lat_asymmetries == 1
        Storing_STABILITY_DATA4B = Conduct_Stability_Study4B(conditions,conv_UNITS,Storing_STABILITY_DATA1,OUTPUT_read_XLSX_mod,...
            Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,case_AC);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA4B = dummy;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Lateral Directional Trim - Accelerations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_lat_accelerations == 1
        Storing_STABILITY_DATA4C = Conduct_Stability_Study4C(conditions,conv_UNITS,Storing_STABILITY_DATA1,OUTPUT_read_XLSX_mod,...
            Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,case_AC);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA4C = dummy;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Lateral Directional Trim - Including Tabs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_lat_Trim_TAB == 1
        Storing_STABILITY_DATA4D = Conduct_Stability_Study4D(conditions,conv_UNITS,Storing_STABILITY_DATA1,OUTPUT_read_XLSX_mod,...
            Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,case_AC,AC_CONFIGURATION);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA4D = dummy;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stability Turning 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Turning == 1
        Storing_STABILITY_DATA5 = Conduct_Stability_Study5(conv_UNITS,Storing_STABILITY_DATA1,OUTPUT_read_XLSX_mod,...
            Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA);
    %   Plot_Options = Storing_STABILITY_DATA5.Plot_Options;
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_STABILITY_DATA5 = dummy;
    end

else
    Plot_Options.dummy = 1;
    dummy = 1;
    Storing_STABILITY_DATA1 = dummy;
    Storing_STABILITY_DATA2 = dummy;
    Storing_STABILITY_DATA2B = dummy;
    Storing_STABILITY_DATA3 = dummy;
    Storing_STABILITY_DATA4A = dummy;
    Storing_STABILITY_DATA4B = dummy;
    Storing_STABILITY_DATA4C = dummy;
    Storing_STABILITY_DATA4D = dummy;
    Storing_STABILITY_DATA5 = dummy;
end



