%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
format long
path(pathdef)
% clc
% Call function that Defines Generic Path so that files can be organized in folders
get_add_path

% Units conversion
conv_UNITS = conversion_UNITS;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
g = conv_UNITS.g;

%% Initializes figures
Fig = 0;

%% Prompts User for reading process
[case_AC OUTPUT_read_XLSX] = Initial_prompt_advanced;

%% MATLAB Compatibility
% Flag that determines MATLAB incompatibility (For example functions that do not work in MATLAB 2017 )
MATLAB_in = OUTPUT_read_XLSX.MATLAB_flags.MATLAB_in;
CHECK_Efficiency = OUTPUT_read_XLSX.MATLAB_flags.CHECK_Efficiency;
Detailed_Profile = OUTPUT_read_XLSX.MATLAB_flags.Detailed_Profile;

%% Defines options for the plots
[mark_Type COLOR_scheme, Plot_Options] = Generate_Plot_Options_init(OUTPUT_read_XLSX);
Plot_Options.MATLAB_in = MATLAB_in;

%% Estimation of prop diameter, just preliminary for geometric conditions
% Stores the flags
% Propulsive_flags.type_battery = type_battery; %
% Propulsive_flags.Engine_loc = Engine_loc; %
% Propulsive_flags.Engine_conf = Engine_conf; %
prop_known = OUTPUT_read_XLSX.Propulsive_flags.prop_known; %
prop_properties = OUTPUT_read_XLSX.Propulsive_flags.prop_properties; %

% type_battery used
% case 1 LiFePO4
% case 2 LiPo
% case 3 FuelCells
type_battery = OUTPUT_read_XLSX.Propulsive_flags.type_battery; %
AC_CONFIGURATION.type_battery = type_battery;
alpha = 0;
beta = 0;

if prop_known == 1
    number = prop_properties/100;
    integ = floor(number);
    fract = (number-integ)*10;
    D_prop = integ*2.54/100;
    OUTPUT_read_XLSX.Propulsive_flags.D_prop = D_prop;
else
    D_prop = OUTPUT_read_XLSX.Propulsive_flags.D_prop; %
end

%% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% AC_type = 6 - 2 surface: cannard + wing + VTP
% AC_type = OUTPUT_read_XLSX.AC_Data_flags.AC_type;
[AC_CONFIGURATION] = Generation_AC_configuration(conv_UNITS,OUTPUT_read_XLSX,case_AC);
AC_CONFIGURATION.type_battery = OUTPUT_read_XLSX.Propulsive_flags.type_battery;

% Propulsion data for combustion Engine
Propulsion = Propulsion_Data_CE(AC_CONFIGURATION);

Modifications.Dxcg = -0.00;
MISSIONS_STUDY = OUTPUT_read_XLSX.STUDY_flags.MISSIONS_STUDY; % Conducts Mission Studies
if MISSIONS_STUDY == 1
    missions = define_mission_DATA(case_AC,OUTPUT_read_XLSX);
    N_missions = length(missions);
    for n=1:N_missions
        mission_actual = missions(n);
%         [conditions OUTPUT_read_XLSX] = select_mission_DATA(case_AC,OUTPUT_read_XLSX,mission_actual,conv_UNITS);
%         Saving_data_OUTPUT_read_XLSX(conditions,OUTPUT_read_XLSX);
        % Stores the mission data
        OUTPUT_read_XLSX_VEC{n} = OUTPUT_read_XLSX;
        [Plot_Options,Storing_GEO_DATA_1{n},Storing_WEIGHT_DATA_1{n},Storing_AERO_DATA_1{n},Storing_PERFORMANCE_DATA_1{n},Storing_STABILITY_DATA_1{n},...
            Storing_STABILITY_DATA_2{n},Storing_STABILITY_DATA_2B{n},Storing_STABILITY_DATA_3{n},Storing_STABILITY_DATA_4{n},Storing_STABILITY_DATA_5{n},...
            Storing_PROPULSION_DATA{n},OUTPUT_read_XLSX_VEC{n}] = conduct_Analysis(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_VEC{n},...
            case_AC,Fig,Plot_Options,mission_actual,Propulsion,Modifications);

        %% Defines options for the plots
        %% Generate Plots
        [M_alpha_cero{n},V_alpha_cero{n}] = GENERATE_PLOTS(OUTPUT_read_XLSX_VEC{n},Plot_Options,Fig,Storing_GEO_DATA_1{n},Storing_WEIGHT_DATA_1{n},Storing_AERO_DATA_1{n},...
            Storing_STABILITY_DATA_1{n},Storing_STABILITY_DATA_2{n},Storing_STABILITY_DATA_2B{n},Storing_STABILITY_DATA_3{n},Storing_STABILITY_DATA_4{n},Storing_STABILITY_DATA_5{n},...
            Storing_PERFORMANCE_DATA_1{n},Storing_PROPULSION_DATA{n},conv_UNITS,AC_CONFIGURATION,case_AC);
    end
    % Writting the data
    for n=1:N_missions
        mission_actual = missions(n);
        if OUTPUT_read_XLSX.MATLAB_flags.write_data == 1
            Write_DATA_complete(Storing_GEO_DATA_1{n},Storing_WEIGHT_DATA_1{n},Storing_AERO_DATA_1{n},Storing_PERFORMANCE_DATA_1{n},...
                Storing_STABILITY_DATA_1{n},Storing_STABILITY_DATA_2{n},Storing_STABILITY_DATA_2B{n},Storing_STABILITY_DATA_3{n},Storing_STABILITY_DATA_4{n},...
                Storing_STABILITY_DATA_5{n},Storing_PROPULSION_DATA{n},M_alpha_cero{n},V_alpha_cero{n},conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX);
        end
    end
else
    % Stores the mission data
    n = 1;
    mission_actual = 1;
    conditions = 1;
    Saving_data_OUTPUT_read_XLSX(conditions,OUTPUT_read_XLSX);
%     conditions.m_TOW = Storing_WEIGHT_DATA.Weight_tier.m_TOW;
    [Plot_Options,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_AERO_DATA_1,Storing_PERFORMANCE_DATA_1,Storing_STABILITY_DATA_1,...
        Storing_STABILITY_DATA_2,Storing_STABILITY_DATA_2B,Storing_STABILITY_DATA_3,Storing_STABILITY_DATA_4,Storing_STABILITY_DATA_5,...
        Storing_PROPULSION_DATA,OUTPUT_read_XLSX] = conduct_Analysis(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX,case_AC,Fig,Plot_Options,...
        mission_actual,Propulsion,Modifications);  

    %% Defines options for the plots
    %% Generate Plots
    [M_alpha_cero,V_alpha_cero] = GENERATE_PLOTS(OUTPUT_read_XLSX,Plot_Options,Fig,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_AERO_DATA_1,...
        Storing_STABILITY_DATA_1,Storing_STABILITY_DATA_2,Storing_STABILITY_DATA_2B,Storing_STABILITY_DATA_3,Storing_STABILITY_DATA_4,Storing_STABILITY_DATA_5,...
        Storing_PERFORMANCE_DATA_1,Storing_PROPULSION_DATA,conv_UNITS,AC_CONFIGURATION,case_AC);

    if OUTPUT_read_XLSX.MATLAB_flags.write_data == 1
        Write_DATA_complete(Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_AERO_DATA_1,Storing_PERFORMANCE_DATA_1,Storing_STABILITY_DATA_1,...
            Storing_STABILITY_DATA_2,Storing_STABILITY_DATA_2B,Storing_STABILITY_DATA_3,Storing_STABILITY_DATA_4,Storing_STABILITY_DATA_5,...
            Storing_PROPULSION_DATA,M_alpha_cero,V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX);
    end

end

if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Saving_DATA_complete(Plot_Options,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_AERO_DATA_1,Storing_PERFORMANCE_DATA_1,Storing_STABILITY_DATA_1,...
        Storing_STABILITY_DATA_2,Storing_STABILITY_DATA_2B,Storing_STABILITY_DATA_3,Storing_STABILITY_DATA_4,Storing_STABILITY_DATA_5,...
        Storing_PROPULSION_DATA,M_alpha_cero,V_alpha_cero,case_AC);
end

%% Checks Efficiency of code
Check_efficiency(CHECK_Efficiency,Detailed_Profile)