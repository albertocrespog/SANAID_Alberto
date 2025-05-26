function Storing_STABILITY_DATA_3 = Conduct_Stability_Study3_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
            OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar,Storing_GEO_DATA,Storing_WEIGHT_DATA,...
            Storing_AERO_DATA,case_AC,Plot_Options,Modifications,Storing_STABILITY_DATA_1,filenameS)

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

%% XCG range
x_XCG_fwd = Storing_STABILITY_DATA_1.TRIM_RESULTS.x_XCG_fwd;
x_XCG_rwd = Storing_STABILITY_DATA_1.TRIM_RESULTS.x_XCG_rwd;
N_x_XCG_VAR = 100;
x_XCG_VAR = linspace(x_XCG_fwd,x_XCG_rwd*1.1,N_x_XCG_VAR);
Plot_Options.x_XCG_VAR = x_XCG_VAR;

conditions.study_var_xcg = 1;
% Wing loading for nominal flight
n = 1;
conditions.n = n;

conditions.V_n = 0;        
for i=1:N_x_XCG_VAR
    % Modification of the propulsion arms in X direction
    x_eng_xbar = Geo_tier.x_eng_xbar;
    % Determine the propulsive arms - Positive for engine behind Xcg
    x_d_T = x_eng_xbar - conditions.x_XCG; % Positive for engine behind Xcg
    % Storing DATA
    Geo_tier.x_d_T = x_d_T;
%     conditions.m_TOW = m_TOW;
    conditions.V = Performance.V;
    conditions.x_XCG = x_XCG_VAR(i);
    gamma = 0;
    conditions.gamma = gamma;
    
    [TRIM_RESULTS_calc_var_xcg,Trim_ITER_calc_var_xcg,Stab_Der_calc_var_xcg,Stab_Der_parts_calc_var_xcg,Propulsion_calc_var_xcg] = ...
        Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
        Body_Geo,Design_criteria,Posicion_Palanca,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar);
    TRIM_RESULTS_var_XCG{i} = TRIM_RESULTS_calc_var_xcg;
    Trim_ITER_var_XCG{i} = Trim_ITER_calc_var_xcg;
    Stab_Der_var_XCG{i} = Stab_Der_calc_var_xcg;
    Stab_Der_parts_var_XCG{i} = Stab_Der_parts_calc_var_xcg;
    Propulsion_var_xcg{i} = Propulsion_calc_var_xcg;
    
    %% Dynamic Stability Analysis
    StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
    Variable_Study = 1;
    %ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
        show_messages_screen = 0;
        Stab_Dyn_Long_var_XCG{i} = longitudinal_analysis(Performance,Stab_Der_calc_var_xcg,conv_UNITS,StabilityModel,...
            Weight_tier,TRIM_RESULTS_calc_var_xcg,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
    else
        Dummy = 0;
        Stab_Dyn_Long_var_XCG = Dummy;
    end
    %ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
    if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
        show_messages_screen = 0;
        Stab_Dyn_LatDir_var_XCG{i} = lateral_directional_analysis(Performance,Stab_Der_calc_var_xcg,conv_UNITS,StabilityModel,...
            Weight_tier,TRIM_RESULTS_calc_var_xcg,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
    else
        Dummy = 0;
        Stab_Dyn_LatDir_var_XCG = Dummy;
    end
end
%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_STABILITY_DATA_3 = Saving_data_Stability_varXCG(TRIM_RESULTS_var_XCG,Trim_ITER_var_XCG,Stab_Der_var_XCG,Stab_Der_parts_var_XCG,...
        Stab_Dyn_Long_var_XCG,Stab_Dyn_LatDir_var_XCG,case_AC,Propulsion_var_xcg,Plot_Options,OUTPUT_read_XLSX,filenameS);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_3 = dummy;
end
