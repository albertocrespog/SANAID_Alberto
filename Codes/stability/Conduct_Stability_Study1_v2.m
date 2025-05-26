function Storing_STABILITY_DATA_1 = Conduct_Stability_Study1_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
            OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA,Storing_WEIGHT_DATA,...
            Storing_AERO_DATA,case_AC,Plot_Options,Modifications,filenameS)

XCG_FF = OUTPUT_read_XLSX.Stability_flags.XCG_FF;
Dxcg = Modifications.Dxcg;

Geo_tier = Storing_GEO_DATA.Geo_tier;
Body_Geo = Storing_GEO_DATA.Body_Geo;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Aero = Storing_AERO_DATA.Aero;
Aero_TH = Storing_AERO_DATA.Aero_TH;
Design_criteria = Storing_AERO_DATA.Design_criteria;
Performance = Storing_AERO_DATA.Performance;
DATA_Ae = Storing_AERO_DATA.DATA_Ae;

% Wing loading for nominal flight
n = 1;
conditions.n = n;
gamma = 0;
conditions.gamma = gamma;
% conf_missiles = conditions.conf_missiles;
% % Calculates the variable xCG as a functin of weight
% if OUTPUT_read_XLSX.Weights_flags.get_shift_XCG_variation == 1
%     conditions = get_Weights_Range(case_AC,conf_missiles);
%     m_TOW = max(Weight_Range);
%     x_XCG = get_Variable_x_XCG(m_TOW,case_AC,conf_missiles) + Dxcg;% - 0.078381885757103;
%     conditions.m_TOW = m_TOW;
%     conditions.x_XCG = x_XCG;    
% else
%     m_TOW = Weight_tier.m_TOW;
%     conditions.m_TOW = m_TOW;
% end
m_TOW = conditions.m_TOW;


%     conditions.m_TOW = m_TOW;

%% Updates the Velocity Vector according to the minimum speeds
CL_max_ac = Aero.CL_max_ac;
V_stall= sqrt(2*m_TOW*conv_UNITS.g/(Geo_tier.S_ref*Performance_preliminar.rho*CL_max_ac));
V_min_ope = 1.3*V_stall;
% conditions.V = V_min_ope;

StabilityModel = OUTPUT_read_XLSX.Stability_flags;
[TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Propulsion_values, Effects] = ...
    Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
    Body_Geo,Design_criteria,Posicion_Palanca,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar);

TRIM_RESULTS
Trim_ITER

StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
% Defines that only one case is conducted
Variable_Study = 0;

%ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
    show_messages_screen = OUTPUT_read_XLSX.MATLAB_flags.show_messages_screen;
    Stab_Dyn_Long = longitudinal_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
        Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen,filenameS);
else
    Dummy = 0;
    Stab_Dyn_Long = Dummy;
end
%ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
    show_messages_screen = OUTPUT_read_XLSX.MATLAB_flags.show_messages_screen;
    Stab_Dyn_LatDir = lateral_directional_analysis(Performance,Stab_Der,conv_UNITS,StabilityModel,...
        Weight_tier,TRIM_RESULTS,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen,filenameS);
else
    Dummy = 0;
    Stab_Dyn_LatDir = Dummy;
end

% Saving variables
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_STABILITY_DATA_1 = Saving_data_Stability_TRIM(TRIM_RESULTS,Trim_ITER,Stab_Der,Stab_Der_parts,Stab_Dyn_Long,Stab_Dyn_LatDir,Propulsion_values,OUTPUT_read_XLSX,Effects,filenameS);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_1 = dummy;
end