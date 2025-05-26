function [OUTPUT_read_XLSX_mod, Storing_PERFORMANCE_DATA] = Conduct_Performance_Study_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
    OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Conducts Performance optimization %
% Performance nominal
if OUTPUT_read_XLSX_mod.STUDY_flags.ANALYSIS_PERFORMANCE_AP == 1
    Storing_PERFORMANCE_DATA_1 = Conduct_Performance_Optimization_v2(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
        OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA_1 = dummy;
end

% Performance Analysis integrated with AP Codes varying Conditions
if OUTPUT_read_XLSX_mod.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var == 1
    % Perfo1 Variable Mass & Speed Performance Studies - Electric
    % Perfo2 Variable Mass & Speed Performance Studies - Fuel
    % Perfo3 Performance Variable Speed - Altitude
    % Perfo4 Glide Performance Variable Speed - Altitude
    % Perfo5 Glide Performance Variable - Altitude for MAx Glide
    % Perfo6 Variable Speed Performance Studies
    % Perfo7 Variable Mass Performance Studies

    if OUTPUT_read_XLSX_mod.STUDY_flags.Perfo1 == 1
        Storing_PERFORMANCE_DATA_21 = Conduct_Performance_Optimization_Perfo1_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_PERFORMANCE_DATA_21 = dummy;
    end

    if OUTPUT_read_XLSX_mod.STUDY_flags.Perfo2 == 1
        Storing_PERFORMANCE_DATA_22 = Conduct_Performance_Optimization_Perfo2_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_PERFORMANCE_DATA_22 = dummy;
    end

        % Performance Analysis integrated with AP Codes varying Conditions
    %% Glide Performancef
    if OUTPUT_read_XLSX_mod.STUDY_flags.Perfo3 == 1
        Storing_PERFORMANCE_DATA_23 = Conduct_Performance_Optimization_Perfo3_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_PERFORMANCE_DATA_23 = dummy;
    end

    % Performance Analysis integrated with AP Codes varying Conditions
    %% Glide Performance
    if OUTPUT_read_XLSX_mod.STUDY_flags.Perfo4 == 1
        Storing_PERFORMANCE_DATA_24 = Conduct_Performance_Optimization_Perfo4_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_PERFORMANCE_DATA_24 = dummy;
    end

        % Performance Analysis integrated with AP Codes varying Conditions
    %% MAx Glide Performance
    if OUTPUT_read_XLSX_mod.STUDY_flags.Perfo5 == 1
        Storing_PERFORMANCE_DATA_25 = Conduct_Performance_Optimization_Perfo5_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_PERFORMANCE_DATA_25 = dummy;
    end

    % Performance Analysis integrated with AP Codes varying Conditions
    if OUTPUT_read_XLSX_mod.STUDY_flags.Perfo6 == 1
        Storing_PERFORMANCE_DATA_26 = Conduct_Performance_Optimization_Perfo6_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_PERFORMANCE_DATA_26 = dummy;
    end

    % Performance Analysis integrated with AP Codes varying Conditions
    if OUTPUT_read_XLSX_mod.STUDY_flags.Perfo7 == 1
        Storing_PERFORMANCE_DATA_27 = Conduct_Performance_Optimization_Perfo7_v1(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
            OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA,filenameS);
    else
        % Stores dummy variable for continuity
        dummy = 1;
        Storing_PERFORMANCE_DATA_27 = dummy;
    end

    % Storing DATA
    % Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_1 = Storing_PERFORMANCE_DATA_1;
    % Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_2 = Storing_PERFORMANCE_DATA_21;
    % Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_2B = Storing_PERFORMANCE_DATA_22;
    % Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_3 = Storing_PERFORMANCE_DATA_26;
    % Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_4 = Storing_PERFORMANCE_DATA_24;
    % Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_4C = Storing_PERFORMANCE_DATA_23;
    % Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_5 = Storing_PERFORMANCE_DATA_25;
    % Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_27 = Storing_PERFORMANCE_DATA_27;

    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_1 = Storing_PERFORMANCE_DATA_1;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_21 = Storing_PERFORMANCE_DATA_21;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_22 = Storing_PERFORMANCE_DATA_22;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_23 = Storing_PERFORMANCE_DATA_23;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_24 = Storing_PERFORMANCE_DATA_24;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_25 = Storing_PERFORMANCE_DATA_25;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_26 = Storing_PERFORMANCE_DATA_26;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_27 = Storing_PERFORMANCE_DATA_27;

else
    dummy = 1;
    Storing_PERFORMANCE_DATA_21 = dummy;
    Storing_PERFORMANCE_DATA_22 = dummy;
    Storing_PERFORMANCE_DATA_26 = dummy;
    Storing_PERFORMANCE_DATA_24 = dummy;
    Storing_PERFORMANCE_DATA_23 = dummy;
    Storing_PERFORMANCE_DATA_25 = dummy;
    Storing_PERFORMANCE_DATA_27 = dummy;

    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_1 = Storing_PERFORMANCE_DATA_1;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_21 = Storing_PERFORMANCE_DATA_21;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_22 = Storing_PERFORMANCE_DATA_22;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_23 = Storing_PERFORMANCE_DATA_23;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_24 = Storing_PERFORMANCE_DATA_24;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_25 = Storing_PERFORMANCE_DATA_25;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_26 = Storing_PERFORMANCE_DATA_26;
    Storing_PERFORMANCE_DATA.Storing_PERFORMANCE_DATA_27 = Storing_PERFORMANCE_DATA_27;
end