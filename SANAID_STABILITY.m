%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SANAID - Stability ANAlysis Interactive Design Tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
format long
tic

path(pathdef)
% Call function that Defines Generic Path so that files can be organized in folders
get_add_path_v2

% Units conversion
conv_UNITS = conversion_UNITS;

%% Initializes figures
Fig = 0;

%% Prompts User for reading process
[case_AC OUTPUT_read_XLSX] = Initial_prompt_advanced_v3;

%% Message Display
Message0 = '-------------------------------------';
Message1 = 'Finished Reading DATA';
disp(Message0)
disp(Message1)
disp(Message0)

%% MATLAB Compatibility
% Flag that determines MATLAB incompatibility (For example functions that do not work in MATLAB 2017 )
MATLAB_in = OUTPUT_read_XLSX.MATLAB_flags.MATLAB_in;
CHECK_Efficiency = OUTPUT_read_XLSX.MATLAB_flags.CHECK_Efficiency;
Detailed_Profile = OUTPUT_read_XLSX.MATLAB_flags.Detailed_Profile;
if CHECK_Efficiency == 1
    profile on
end

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
[AC_CONFIGURATION] = Generation_AC_configuration(conv_UNITS,OUTPUT_read_XLSX,case_AC);
AC_CONFIGURATION.type_battery = OUTPUT_read_XLSX.Propulsive_flags.type_battery;

% Propulsion data for combustion Engine
Propulsion = Propulsion_Data_CE(AC_CONFIGURATION);

%% Selects the FOLDER file location (Unification of address)
[filenameS,Sheet_AC_var] = selection_filename_location(case_AC);

Modifications.Dxcg = 0.00;
MISSIONS_STUDY = OUTPUT_read_XLSX.STUDY_flags.MISSIONS_STUDY; % Conducts Mission Studies

if MISSIONS_STUDY == 1
    missions = define_mission_DATA(case_AC,OUTPUT_read_XLSX);
    N_missions = length(missions);
    for n=1:N_missions
    
        %% DEFINES MISSION BEIN ANALYZED
        mission_actual = missions(n);

        %% CONDUCT ANALYSIS
        % Message Display
        st1 = 'Start Studies - Mission: ';
        st1B = 'Ends Studies - Mission: ';
        st2 = num2str(mission_actual');
        Message2 = strcat(st1,st2);
        Message2B = strcat(st1B,st2);
        disp(Message0)
        disp(Message2)
        disp(Message0)

        % Stores the mission data
        OUTPUT_read_XLSX_VEC{n} = OUTPUT_read_XLSX;
        [Plot_Options,Storing_DATA{n},OUTPUT_read_XLSX_VEC{n}] = conduct_Analysis_v2(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX_VEC{n},...
            case_AC,Fig,Plot_Options,mission_actual,Propulsion,Modifications,filenameS);

        % Message Display
        disp(Message0)
        disp(Message2B)
        disp(Message0)

        %% Generate Plots       
        % Message Display
        st1 = 'Start PLOTS - Mission: ';
        st1B = 'Ends PLOTS - Mission: ';
        st2 = num2str(mission_actual');
        Message4 = strcat(st1,st2);
        Message4B = strcat(st1B,st2);
        disp(Message0)
        disp(Message4)
        disp(Message0)

        [M_alpha_cero{n},V_alpha_cero{n}] = GENERATE_PLOTS_v2(OUTPUT_read_XLSX_VEC{n},Storing_DATA{n},Plot_Options,Fig,conv_UNITS,AC_CONFIGURATION,case_AC,filenameS);

        % Message Display
        disp(Message0)
        disp(Message4B)
        disp(Message0)

    end
    %% Writting the data
    for n=1:N_missions
        mission_actual = missions(n);
        if OUTPUT_read_XLSX.MATLAB_flags.write_data == 1
            
            %% WRITE DATA
            % Message Display
            st1 = 'Start Writing Data - Mission: ';
            st1B = 'Ends Writing Data - Mission: ';
            st2 = num2str(mission_actual');
            Message3 = strcat(st1,st2);
            Message3B = strcat(st1B,st2);
            disp(Message0)
            disp(Message3)
            disp(Message0)

            Write_DATA_complete_v2(Storing_DATA{n},M_alpha_cero{n},V_alpha_cero{n},conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filenameS,Sheet_AC_var);

            %% Message Display
            disp(Message0)
            disp(Message3B)
            disp(Message0)
        end
    end

    %% SAVING DATAFOR  multiple iterations
    if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
        %% Message Display
        st1 = 'Start Saving Data - Mission: ';
        st1B = 'Ends Saving Data - Mission: ';
        st2 = num2str(mission_actual');
        Message4 = strcat(st1,st2);
        Message4B = strcat(st1B,st2);
        disp(Message0)
        disp(Message4)
        disp(Message0)

        %% SAVING DATA
        Saving_DATA_complete(Plot_Options,Storing_DATA{n},M_alpha_cero,V_alpha_cero,case_AC,OUTPUT_read_XLSX,filenameS,Sheet_AC_var);
        
        % Message Display
        disp(Message0)
        disp(Message4B)
        disp(Message0)
    end

else % SINGLE ANALYSIS
    %% Stores the mission data
    n = 1;
    mission_actual = 1;
    conditions = 1;
    Saving_data_OUTPUT_read_XLSX(conditions,OUTPUT_read_XLSX,filenameS);

    % Message Display
    st1 = 'Start Studies - Mission: ';
    st1B = 'Start Studies - Mission: ';
    st2 = num2str(mission_actual');
    Message2 = strcat(st1,st2);
    Message2B = strcat(st1B,st2);
    disp(Message0)
    disp(Message2)
    disp(Message0)

    %% CONDUCT ANALYSIS
    [Plot_Options,Storing_DATA,OUTPUT_read_XLSX] = conduct_Analysis_v2(conv_UNITS,AC_CONFIGURATION,OUTPUT_read_XLSX,...
        case_AC,Fig,Plot_Options,mission_actual,Propulsion,Modifications,filenameS);
    %% Message Display
    disp(Message0)
    disp(Message2B)
    disp(Message0)

    %% Generate Plots
    % Message Display
    st1 = 'Start PLOTS - Mission: ';
    st1B = 'Ends PLOTS - Mission: ';
    st2 = num2str(mission_actual');
    Message4 = strcat(st1,st2);
    Message4B = strcat(st1B,st2);
    disp(Message0)
    disp(Message4)
    disp(Message0)
    
    %% Generate Plots
    [M_alpha_cero,V_alpha_cero]  = GENERATE_PLOTS_v2(OUTPUT_read_XLSX,Storing_DATA,Plot_Options,Fig,conv_UNITS,AC_CONFIGURATION,case_AC,filenameS);
    % Message Display
    disp(Message0)
    disp(Message4B)
    disp(Message0)

    %% Writting the data
    if OUTPUT_read_XLSX.MATLAB_flags.write_data == 1

        %% Message Display
        st1 = 'Start Writing Data - Mission: ';
        st1 = 'Ends Writing Data - Mission: ';
        st2 = num2str(mission_actual');
        Message3 = strcat(st1,st2);
        Message3B = strcat(st1B,st2);
        disp(Message0)
        disp(Message3)
        disp(Message0)

        %% WRITE DATA
        Write_DATA_complete_v2(Storing_DATA,M_alpha_cero{n},V_alpha_cero,conv_UNITS,mission_actual,n,case_AC,OUTPUT_read_XLSX,filenameS,Sheet_AC_var);

        %% Message Display
        disp(Message0)
        disp(Message3B)
        disp(Message0)
    end

    %% SAVING DATAFOR For only one iteration
    if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
        %% Message Display
        st1 = 'Start Saving Data - Mission: ';
        st1B = 'Ends Saving Data - Mission: ';
        st2 = num2str(mission_actual');
        Message4 = strcat(st1,st2);
        Message4B = strcat(st1B,st2);
        disp(Message0)
        disp(Message4)
        disp(Message0)
        
        %% SAVING DATA
        Saving_DATA_complete(Plot_Options,Storing_DATA,M_alpha_cero,V_alpha_cero,case_AC,OUTPUT_read_XLSX,filenameS,Sheet_AC_var);
    
        % Message Display
        disp(Message0)
        disp(Message4B)
        disp(Message0)
    end
end

%% Checks Efficiency of code
Check_efficiency(CHECK_Efficiency,Detailed_Profile)
toc
