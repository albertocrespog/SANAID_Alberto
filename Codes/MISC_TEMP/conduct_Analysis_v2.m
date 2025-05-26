function [Plot_Options,Storing_DATA,OUTPUT_read_XLSX_mod] = conduct_Analysis_v2(conv_UNITS,AC_CONFIGURATION,...
    OUTPUT_read_XLSX,case_AC,Fig,Plot_Options,mission_actual,Propulsion,Modifications,filenameS)

%% DEFINES MISSION STUDY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
MISSIONS_STUDY = OUTPUT_read_XLSX.STUDY_flags.MISSIONS_STUDY; % Conducts Mission Studies

if MISSIONS_STUDY == 1
    [conditions, OUTPUT_read_XLSX_mod] = select_mission_DATA(case_AC,OUTPUT_read_XLSX,mission_actual,conv_UNITS,Modifications);
%     [conditions OUTPUT_read_XLSX] = select_mission_DATA(case_AC,OUTPUT_read_XLSX,mission_actual,conv_UNITS);
    Saving_data_OUTPUT_read_XLSX(conditions,OUTPUT_read_XLSX_mod,filenameS);
else
    dummy = 1;
    conditions.dummy = dummy; 
    OUTPUT_read_XLSX_mod = OUTPUT_read_XLSX;
end

Saving_AC_CONFIGURATION(AC_CONFIGURATION,OUTPUT_read_XLSX_mod,filenameS);
% Dxcg = Modifications.Dxcg; 

%% Geometric Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Message Display
st1 = 'Start Generate Geometry - Mission: ';
st1B = 'End Generate Geometry - Mission: ';
st2 = num2str(mission_actual');
Message0 = '-------------------------------------';
Message4 = strcat(st1,st2);
Message4B = strcat(st1B,st2);
disp(Message0)
disp(Message4)
disp(Message0)

% Geo_input_tier = Generation_Input_Geometric_Data_v2(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod);
Geo_input_tier = Generation_Input_Geometric_Data_v3(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod);
% Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod,case_AC);
% Geo_tier = Generation_Geometric_Data_v4(Geo_input_tier,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod,case_AC);
Geo_tier = Generation_Geometric_Data_v6(Geo_input_tier,conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_mod,case_AC);

% Message Display
disp(Message0)
disp(Message4B)
disp(Message0)

%% BODY GEOMETRY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Message Display
st1 = 'Start Read Fuselage DATA - Mission: ';
st1B = 'End Read Fuselage DATA - Mission: ';
st2 = num2str(mission_actual');
Message4 = strcat(st1,st2);
Message4B = strcat(st1B,st2);
disp(Message0)
disp(Message4)
disp(Message0)

% Read data already stored to accelerate
if OUTPUT_read_XLSX.Fuselage_flags.Use_Storing_DATA == 0
    % Generation of Geometry
    [Fig,Body_Geo,meshData] = Generation_Fuselage_Data_old(Geo_tier,OUTPUT_read_XLSX_mod,Fig);% Defines Propulsion DATA
    Storing_GEO_DATA_1 = Saving_data_Geo(Geo_tier,Body_Geo,meshData,OUTPUT_read_XLSX_mod,filenameS);
else
    [Storing_GEO_DATA_1,Body_Geo,meshData] = Read_data_Geo(OUTPUT_read_XLSX_mod,filenameS,Geo_tier);
end

% Message Display
disp(Message0)
disp(Message4B)
disp(Message0)

%% Defines Estimation of Weights according different methods densities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Message Display
st1 = 'Start Weight Estimation - Mission: ';
st2 = num2str(mission_actual');
Message4 = strcat(st1,st2);
disp(Message0)
disp(Message4)
disp(Message0)

Weight_tier = Generation_Weight_Data(Geo_tier,Body_Geo,AC_CONFIGURATION,conv_UNITS,OUTPUT_read_XLSX_mod);
% Read data already stored to accelerate
if MISSIONS_STUDY == 1
    Weight_tier.m_TOW = conditions.m_TOW;
    Storing_WEIGHT_DATA_1 = Saving_data_Weight(Weight_tier,OUTPUT_read_XLSX_mod,filenameS);
else
    Storing_WEIGHT_DATA_1 = Saving_data_Weight(Weight_tier,OUTPUT_read_XLSX_mod,filenameS);
end

%% Generates the file for Aerodynamic Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Indentifies preliminary performance results
Performance_preliminar = Generate_Performance_preliminary(OUTPUT_read_XLSX_mod);
 
%% AERODYNAMIC STUDY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Message Display
st1 = 'Start Aerodynamic Study - Mission: ';
st1B = 'End Aerodynamic Study - Mission: ';
st2 = num2str(mission_actual');
Message4 = strcat(st1,st2);
Message4B = strcat(st1B,st2);
disp(Message0)
disp(Message4)
disp(Message0)

if OUTPUT_read_XLSX_mod.STUDY_flags.AERODYNAMIC_STUDY == 1
    [Storing_AERO_DATA_1] = Conduct_Aerodynamic_Study(OUTPUT_read_XLSX_mod,AC_CONFIGURATION,conv_UNITS,Performance_preliminar,...
        case_AC,Storing_WEIGHT_DATA_1,Storing_GEO_DATA_1,conditions,filenameS);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_AERO_DATA_1 = dummy;
    % Warning
    warning_areo = 'WARNING!!! No results generated, please Make sure that AERODYNAMIC_STUDY == 1 - Code in PAUSE';
    disp(warning_areo)
    pause
end
% Message Display
disp(Message0)
disp(Message4B)
disp(Message0)

%% PROPULSION STUDY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Propulsive Model
[Prop_data] = Generation_Propulsion_Data(AC_CONFIGURATION,OUTPUT_read_XLSX_mod,filenameS);
Storing_PROPULSION_DATA_1.Prop_data = Prop_data;

% Data for analysis of Prop DATA
VECTOR_Prop = Get_Comparing_vectors_Propulsion(case_AC,OUTPUT_read_XLSX_mod);

% Generates the file for Prop Data to be used in the codes
% actualizes the Prop geometry
% Select:
% Propulsion model: SEE read_prop_files_May2020.m FOR DETAILS OF PROP MODEL!!
Flag.APC_model = 1;
Flag.Wind_tunnel_model = 1;
Flag.compare_prop_models = 1;

%% Conducts Propulsion optimization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if OUTPUT_read_XLSX_mod.STUDY_flags.PROP_STUDY == 1
    Conduct_Prop_Optimization
end

%% PERFORMANCE STUDY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if OUTPUT_read_XLSX_mod.STUDY_flags.PERFORMANCE_STUDY == 1
    % Message Display
    st1 = 'Start Performance Study - Mission: ';
    st1B = 'End Performance Study - Mission: ';
    st2 = num2str(mission_actual');
    Message4 = strcat(st1,st2);
    Message4B = strcat(st1B,st2);
    disp(Message0)
    disp(Message4)
    disp(Message0)

    % PROPULSION STUDY
    [OUTPUT_read_XLSX_mod, Storing_PERFORMANCE_DATA] = Conduct_Performance_Study_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
        OUTPUT_read_XLSX_mod,Storing_AERO_DATA_1,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_PROPULSION_DATA_1,filenameS); 

    % Message Display
    disp(Message0)
    disp(Message4B)
    disp(Message0)
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_1 = dummy;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_21 = dummy;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_22 = dummy;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_23 = dummy;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_24 = dummy;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_25 = dummy;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_26 = dummy;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_27 = dummy;

end

%% STABILITY STUDY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY == 1
    %% EVOLUTION DATA
    % Will be used in future version of propulsion
    % Propulsion Generation
    alpha_f = 0*conv_UNITS.D2R;
    beta_f = 0*conv_UNITS.D2R;
    conditions.alpha_f = alpha_f;
    conditions.beta_f = beta_f;
    conditions.h = Storing_AERO_DATA_1.Performance.h;
    conditions.V = Storing_AERO_DATA_1.Performance.V;
    conditions.study_var_xcg = 0;
    conditions.x_XCG = OUTPUT_read_XLSX_mod.InputGeometry_Data_flags.x_XCG;
    conditions.z_XCG = OUTPUT_read_XLSX_mod.InputGeometry_Data_flags.z_XCG;


    %% Selects engione on or off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    % Posicion_Palanca = 0; % Engien off
    % Posicion_Palanca = 1; % Engien on
    Posicion_Palanca = conditions.Posicion_Palanca;

    % Dynamic Stability Formulation
    % Message Display
    st1 = 'Start Stability Study - Mission: ';
    st1B = 'End Stability Study - Mission: ';
    st2 = num2str(mission_actual');
    Message4 = strcat(st1,st2);
    Message4B = strcat(st1B,st2);
    disp(Message0)
    disp(Message4)
    disp(Message0)

    % STABILTY STUDY
    [Storing_STABILITY_DATA, Plot_Options] = Conduct_Stability_MainStudy_v1(AC_CONFIGURATION,case_AC,OUTPUT_read_XLSX_mod,conditions,conv_UNITS,...
        Posicion_Palanca,Modifications,Plot_Options,Propulsion,Prop_data,...
        Performance_preliminar,Storing_AERO_DATA_1,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_PROPULSION_DATA_1,filenameS);

    %% Message Display
    disp(Message0)
    disp(Message4B)
    disp(Message0)

else
    % Stores dummy variable for continuity
    Plot_Options.dummy = 1;
    dummy = 1;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_1 = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_2 = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_2B = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_2C = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_3 = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_4A = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_4B = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_4C = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_4D = dummy;
    Storing_STABILITY_DATA.Storing_STABILITY_DATA_5 = dummy;
end

%% Stores in single Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
Storing_DATA.Storing_GEO_DATA_1 = Storing_GEO_DATA_1;
Storing_DATA.Storing_WEIGHT_DATA_1 = Storing_WEIGHT_DATA_1;
Storing_DATA.Storing_AERO_DATA_1 = Storing_AERO_DATA_1;
% Performance
Storing_DATA.Storing_PROPULSION_DATA_1 = Storing_PROPULSION_DATA_1;
Storing_DATA.Storing_PERFORMANCE_DATA_1 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_1;
Storing_DATA.Storing_PERFORMANCE_DATA_21 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_21;
Storing_DATA.Storing_PERFORMANCE_DATA_22 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_22;
Storing_DATA.Storing_PERFORMANCE_DATA_23 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_23;
Storing_DATA.Storing_PERFORMANCE_DATA_24 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_24;
Storing_DATA.Storing_PERFORMANCE_DATA_25 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_25;
Storing_DATA.Storing_PERFORMANCE_DATA_26 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_26;
Storing_DATA.Storing_PERFORMANCE_DATA_27 = Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_27;
% Stability
Storing_DATA.Storing_STABILITY_DATA_1 = Storing_STABILITY_DATA.Storing_STABILITY_DATA_1;
Storing_DATA.Storing_STABILITY_DATA_2 = Storing_STABILITY_DATA.Storing_STABILITY_DATA_2;
Storing_DATA.Storing_STABILITY_DATA_2B = Storing_STABILITY_DATA.Storing_STABILITY_DATA_2B;
Storing_DATA.Storing_STABILITY_DATA_2C = Storing_STABILITY_DATA.Storing_STABILITY_DATA_2C;
Storing_DATA.Storing_STABILITY_DATA_3 = Storing_STABILITY_DATA.Storing_STABILITY_DATA_3;
Storing_DATA.Storing_STABILITY_DATA_4A = Storing_STABILITY_DATA.Storing_STABILITY_DATA_4A;
Storing_DATA.Storing_STABILITY_DATA_4B = Storing_STABILITY_DATA.Storing_STABILITY_DATA_4B;
Storing_DATA.Storing_STABILITY_DATA_4C = Storing_STABILITY_DATA.Storing_STABILITY_DATA_4C;
Storing_DATA.Storing_STABILITY_DATA_4D = Storing_STABILITY_DATA.Storing_STABILITY_DATA_4D;
Storing_DATA.Storing_STABILITY_DATA_5 = Storing_STABILITY_DATA.Storing_STABILITY_DATA_5;
