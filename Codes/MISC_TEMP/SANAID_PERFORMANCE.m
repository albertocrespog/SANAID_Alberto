function [OUTPUT_read_XLSX_mod, Storing_PERFORMANCE_DATA_1, Storing_PERFORMANCE_DATA_2, Storing_PERFORMANCE_DATA_3] = Conduct_Performance_Study(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
        OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Conducts Performance optimization %
if OUTPUT_read_XLSX_mod.STUDY_flags.PERFORMANCE_STUDY == 1
    % Performance nominal
    Storing_PERFORMANCE_DATA_1 = Conduct_Performance_Optimization_v2(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
        OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA);
    % Performance Analysis integrated with AP Codes varying Conditions
    if OUTPUT_read_XLSX_mod.STUDY_flags.ANALYSIS_PERFORMANCE_AP_var == 1

        if OUTPUT_read_XLSX_mod.STUDY_flags.variable_speed_weight_AP == 1
            Storing_PERFORMANCE_DATA_2 = Conduct_Performance_Optimization_v4(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
                OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA);
        else
            % Stores dummy variable for continuity
            dummy = 1;
            Storing_PERFORMANCE_DATA_2 = dummy;
        end

     % Performance Analysis integrated with AP Codes varying Conditions
        if OUTPUT_read_XLSX_mod.STUDY_flags.variable_speed_AP == 1
            Storing_PERFORMANCE_DATA_3 = Conduct_Performance_Optimization_v5(AC_CONFIGURATION,case_AC,Propulsion,Prop_data,conv_UNITS,...
                OUTPUT_read_XLSX_mod,Storing_AERO_DATA,Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_PROPULSION_DATA);
        else
            % Stores dummy variable for continuity
            dummy = 1;
            Storing_PERFORMANCE_DATA_3 = dummy;
        end
    end
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_PERFORMANCE_DATA_1 = dummy;
    Storing_PERFORMANCE_DATA_2 = dummy;
    Storing_PERFORMANCE_DATA_3 = dummy;
end