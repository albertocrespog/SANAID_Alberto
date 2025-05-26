function [Storing_STABILITY_DATA,Plot_Options] = Conduct_Stability_MainStudy_v1(AC_CONFIGURATION,case_AC,OUTPUT_read_XLSX_mod,conditions,conv_UNITS,...
        Posicion_Palanca,Modifications,Plot_Options,Propulsion,Prop_data,...
            Performance_preliminar,Storing_AERO_DATA_1,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,Storing_PROPULSION_DATA_1,filenameS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Conducts Stability Study

% It is only activated for the Vn study
conditions.V_n = 0;        

%% Regular Stability Study
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Regular == 1
    only_trim = OUTPUT_read_XLSX_mod.Stability_flags.only_trim;
    Storing_STABILITY_DATA_1 = Conduct_Stability_Study1_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,filenameS);

    % Plot_Options = Storing_STABILITY_DATA_1.Plot_Options;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_1 = dummy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability Sensitivity Study for Velocity and Mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_varV_m0 == 1
    only_trim = OUTPUT_read_XLSX_mod.Stability_flags.only_trim;
    Storing_STABILITY_DATA_2 = Conduct_Stability_Study2_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);
    % Variation of V and n
    Storing_STABILITY_DATA_2B = Conduct_Stability_Study2B_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);

    Plot_Options = Storing_STABILITY_DATA_2.Plot_Options;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_2 = dummy;
    Storing_STABILITY_DATA_2B = dummy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability Sensitivity Study for Velocity and Gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_varV_gamma == 1
    only_trim = OUTPUT_read_XLSX_mod.Stability_flags.only_trim;

    % Variation of V and gamma
    Storing_STABILITY_DATA_2C = Conduct_Stability_Study2C_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);

    Plot_Options = Storing_STABILITY_DATA_2C.Plot_Options;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_2C = dummy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability Sensitivity for varying Xcg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Study of variation of Velocity with XCG selected for trim conditions at desired Static Margin
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_var_XCG
    only_trim = OUTPUT_read_XLSX_mod.Stability_flags.only_trim;
    Storing_STABILITY_DATA_3 = Conduct_Stability_Study3_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);
    
    Plot_Options = Storing_STABILITY_DATA_3.Plot_Options;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_3 = dummy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability Lateral Directional Trim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_lat == 1
    Storing_STABILITY_DATA_4A = Conduct_Stability_Study4A_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);
    % Plot_Options = Storing_STABILITY_DATA_4A.Plot_Options;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_4A = dummy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability Lateral Directional Trim - Asymmetries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_lat_asymmetries == 1
    Storing_STABILITY_DATA_4B = Conduct_Stability_Study4B_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);
    % Plot_Options = Storing_STABILITY_DATA_4B.Plot_Options;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_4B = dummy;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability Lateral Directional Trim - Accelerations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_lat_accelerations == 1
    Storing_STABILITY_DATA_4C = Conduct_Stability_Study4C(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_4C = dummy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability Lateral Directional Trim - Including Tabs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Trim_lat_Trim_TAB == 1
    Storing_STABILITY_DATA_4D = Conduct_Stability_Study4D(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_4D = dummy;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability Turning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX_mod.STUDY_flags.STABILITY_STUDY_Turning == 1
    Storing_STABILITY_DATA_5 = Conduct_Stability_Study5_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
        OUTPUT_read_XLSX_mod,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA_1,Storing_WEIGHT_DATA_1,...
        Storing_AERO_DATA_1,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS);
    %   Plot_Options = Storing_STABILITY_DATA5.Plot_Options;
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_5 = dummy;
end

Storing_STABILITY_DATA.Storing_STABILITY_DATA_1 = Storing_STABILITY_DATA_1;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_2 = Storing_STABILITY_DATA_2;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_2B = Storing_STABILITY_DATA_2B;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_2C = Storing_STABILITY_DATA_2C;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_3 = Storing_STABILITY_DATA_3;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_4A = Storing_STABILITY_DATA_4A;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_4B = Storing_STABILITY_DATA_4B;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_4C = Storing_STABILITY_DATA_4C;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_4D = Storing_STABILITY_DATA_4D;
Storing_STABILITY_DATA.Storing_STABILITY_DATA_5 = Storing_STABILITY_DATA_5;

