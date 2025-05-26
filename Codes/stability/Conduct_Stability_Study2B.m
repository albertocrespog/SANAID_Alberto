function Storing_STABILITY_DATA_2B = Conduct_Stability_Study2B(conditions,Prop_data,conv_UNITS,...
    Posicion_Palanca,XCG_FF,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar,case_AC,Plot_Options,...
    Storing_GEO_DATA,Storing_WEIGHT_DATA,Storing_AERO_DATA,Modifications)

Geo_tier = Storing_GEO_DATA.Geo_tier;
Body_Geo = Storing_GEO_DATA.Body_Geo;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;
Aero = Storing_AERO_DATA.Aero;
Aero_TH = Storing_AERO_DATA.Aero_TH;
Design_criteria = Storing_AERO_DATA.Design_criteria;
Performance = Storing_AERO_DATA.Performance;
DATA_Ae = Storing_AERO_DATA.DATA_Ae;

%% Weight range
% m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_systems + Weight_tier.m_batteries/2 + Weight_tier.m_fuel/2);
m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy);
m_mid = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy)/2;
%     m_low = Weight_tier.m_TOW - (Weight_tier.m_energy);
%     m_mid = Weight_tier.m_TOW - (Weight_tier.m_energy)/2;
m_high = Weight_tier.m_TOW;
N_m_VAR = OUTPUT_read_XLSX.Stability_flags.N_m_VAR;
m_VAR = linspace(m_low,m_high,N_m_VAR);

% case 5 - conf_missiles = 5 weight variation with 0 missiles;
% case 4 - conf_missiles = 4 weight variation with 1 missiles;
% case 3 - conf_missiles = 3 weight variation with 2 missiles;
% case 2 - conf_missiles = 2 weight variation with 3 missiles;
% case 1 - conf_missiles = 1 weight variation with 4 missiles;
% conf_missiles = 1;
% conf_missiles = conditions.conf_missiles;

if OUTPUT_read_XLSX.Weights_flags.get_shift_XCG_variation == 1
    Weight_Range = get_Weights_Range(case_AC,conditions,OUTPUT_read_XLSX);
    m_low = min(Weight_Range);
    m_high = max(Weight_Range);
    N_m_VAR = OUTPUT_read_XLSX.Stability_flags.N_m_VAR;
    m_VAR = linspace(m_low,m_high,N_m_VAR);
else
    m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy);
    m_mid = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy)/2;
    m_high = Weight_tier.m_TOW;
    N_m_VAR = OUTPUT_read_XLSX.Stability_flags.N_m_VAR;
    m_VAR = linspace(m_low,m_high,N_m_VAR);
end

%% Speed range
V_low = OUTPUT_read_XLSX.Stability_flags.V_low;
V_high = OUTPUT_read_XLSX.Stability_flags.V_high;
N_V_VAR = OUTPUT_read_XLSX.Stability_flags.N_V_VAR;

n_min = OUTPUT_read_XLSX.Stability_flags.n_min;
n_max = OUTPUT_read_XLSX.Stability_flags.n_max;
N_n_VAR = OUTPUT_read_XLSX.Stability_flags.N_n_VAR;
n_VAR = linspace(n_min,n_max,N_n_VAR);
Plot_Options.n_VAR = n_VAR;

%% Updates the Velocity Vector according to the minimum speeds
CL_max_ac = Aero.CL_max_ac;

V_stall_min_low = sqrt(2*m_low*conv_UNITS.g/(Geo_tier.S_w1*Performance.rho*CL_max_ac));
V_stall_min_high = sqrt(2*m_high*conv_UNITS.g/(Geo_tier.S_w1*Performance.rho*CL_max_ac));
V_min_ope_low = 1.3*V_stall_min_low;
V_min_ope_high = 1.3*V_stall_min_high;

% Correct the minimum speed  from excel
V_low = V_stall_min_low;

% Correct the maximum speed from excel
if V_high < V_min_ope_high
    V_high = V_min_ope_high;
end

V_VAR = linspace(V_low,V_high,N_V_VAR);
Plot_Options.V_VAR = V_VAR;

% reduces the calculations since only for trim conditions
%% Study of variation of Velocity and mass with XCG selected for trim conditions at desired Static Margin
% Calculates the variable xCG as a functin of weight
if OUTPUT_read_XLSX.Weights_flags.get_shift_XCG_variation == 1
    [conditions,OUTPUT_read_XLSX] = get_Variable_x_XCG(m_high,case_AC,conditions,OUTPUT_read_XLSX,Modifications);
    x_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
    z_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;
else
    x_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
    z_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;
end
        
for i=1:N_n_VAR
    n = n_VAR(i);
    % actual values into de variable study
    conditions.x_XCG = conditions.x_XCG;
    conditions.x_XCG_var = conditions.x_XCG;
    % loop
    for j=1:N_V_VAR
        x_XCG_VAR(i,j) = x_XCG;
        z_XCG_VAR(i,j) = z_XCG;
        
        % Determines stall conditions
        conditions.V = V_VAR(j);
        conditions.m_TOW = m_high;
        conditions.n = n;

        V_stall = sqrt(2*m_high*conv_UNITS.g/(Geo_tier.S_w1*Performance_preliminar.rho*CL_max_ac));
        V_TO = 1.2*V_stall;
        V_min_ope = 1.3*V_stall;
        % Stores values into de restrictions so that can be plotted
        Restrictions.V_stall = V_stall;
        Restrictions.V_TO = V_TO;
        Restrictions.V_min_ope = V_min_ope;
        Restrictions.V_min_ope = V_min_ope;
        % Calculations
        [TRIM_RESULTS_calc_var_V_m_v,Trim_ITER_calc_var_V_m_v,Stab_Der_calc_var_V_m_v,Stab_Der_parts_calc_var_V_m_v,Propulsion_calc_var_V_m_v] = ...
            Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
            Body_Geo,Design_criteria,Posicion_Palanca,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar);
        % Storing
        Restrictions_var_V_m{i,j} = Restrictions;
        TRIM_RESULTS_var_V_m{i,j} = TRIM_RESULTS_calc_var_V_m_v;
        Trim_ITER_var_V_m{i,j} = Trim_ITER_calc_var_V_m_v;
        Stab_Der_var_V_m{i,j} = Stab_Der_calc_var_V_m_v;
        Stab_Der_parts_V_m{i,j} = Stab_Der_parts_calc_var_V_m_v;
        Propulsion_var_V_m{i,j} =  Propulsion_calc_var_V_m_v;
        
        StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
        Variable_Study = 1;
        
        %ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
            show_messages_screen = 0;
            Stab_Dyn_Long_var_V_m{i,j} = longitudinal_analysis(Performance,Stab_Der_calc_var_V_m_v,conv_UNITS,StabilityModel,...
                Weight_tier,TRIM_RESULTS_calc_var_V_m_v,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
        else
            Dummy = 0;
            Stab_Dyn_Long_var_V_m = Dummy;
        end
        %ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
            show_messages_screen = 0;
            Stab_Dyn_LatDir_var_V_m{i,j} = lateral_directional_analysis(Performance,Stab_Der_calc_var_V_m_v,conv_UNITS,StabilityModel,...
                Weight_tier,TRIM_RESULTS_calc_var_V_m_v,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
        else
            Dummy = 0;
            Stab_Dyn_LatDir_var_V_m = Dummy;
        end
    end
end

%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_STABILITY_DATA_2B = Saving_data_Stability_varV(TRIM_RESULTS_var_V_m,Trim_ITER_var_V_m,Stab_Der_var_V_m,Stab_Der_parts_V_m,...
        Stab_Dyn_Long_var_V_m,Stab_Dyn_LatDir_var_V_m,case_AC,Propulsion_var_V_m,Restrictions_var_V_m,Plot_Options,OUTPUT_read_XLSX,x_XCG_VAR);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_2B = dummy;
end