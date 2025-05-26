function Storing_STABILITY_DATA_2 = Conduct_Stability_Study2_v2(conditions,Prop_data,Propulsion,conv_UNITS,Posicion_Palanca,...
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

%% Weight range
% m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_systems + Weight_tier.m_batteries/2 + Weight_tier.m_fuel/2);
% m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy);
% m_mid = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy)/2;
%     m_low = Weight_tier.m_TOW - (Weight_tier.m_energy);
%     m_mid = Weight_tier.m_TOW - (Weight_tier.m_energy)/2;
% m_high = Weight_tier.m_TOW;
% N_m_VAR = OUTPUT_read_XLSX.Stability_flags.N_m_VAR;
% m_VAR = linspace(m_low,m_high,N_m_VAR);

%ding Wing load
n = 1;
conditions.n = n;

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

Plot_Options.W_VAR = m_VAR;
Plot_Options.n_VAR = n;

%% Speed range
V_low = OUTPUT_read_XLSX.Stability_flags.V_low;
V_high = OUTPUT_read_XLSX.Stability_flags.V_high;
N_V_VAR = OUTPUT_read_XLSX.Stability_flags.N_V_VAR;

%% Updates the Velocity Vector according to the minimum speeds
CL_max_ac = Aero.CL_max_ac;
% V_stall_ac = sqrt((2*m_TOW*g)/(rho*S_ref*CL_max_ac_ope));
% V_min_ac = Flight_SF*V_stall_ac;
% Performance.V_stall_ac = V_stall_ac;
% Performance.V_min_ac = V_min_ac;

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

Dxcg = Modifications.Dxcg;

% tic
% % reduces the calculations since only for trim conditions
% % Study of variation of Velocity and mass with XCG selected for trim conditions at desired Static Margin
% for i=1:N_V_VAR
%     % actual values into de variable study    
%     conditions.V = V_VAR(i);
%     conditions.x_XCG = conditions.x_XCG;
%     conditions.x_XCG_var = conditions.x_XCG;
%     % loop
%     for j=1:N_m_VAR
%         % Determines stall conditions
%         conditions.m_TOW = m_VAR(j);
% 
%         % Calculates the variable xCG as a functin of weight
%         if OUTPUT_read_XLSX.Weights_flags.get_shift_XCG_variation == 1
% %             OUTPUT = get_Variable_x_XCG(m_VAR(j),case_AC,conditions,OUTPUT_read_XLSX) + Dxcg;
%             [conditions, OUTPUT_read_XLSX] = get_Variable_x_XCG(m_VAR(j),case_AC,conditions,OUTPUT_read_XLSX,Modifications);
% %             OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG + Dxcg;
%             x_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
%             z_XCG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;
%             x_XCG_VAR(i,j) = x_XCG;
%             z_XCG_VAR(i,j) = z_XCG;
%         else
%             dummy = 1;
%             x_XCG_VAR(i,j) = dummy;
%         end
% 
%         V_stall = sqrt(2*m_VAR(j)*conv_UNITS.g/(Geo_tier.S_w1*Performance_preliminar.rho*CL_max_ac));
%         V_TO = 1.2*V_stall;
%         V_min_ope = 1.3*V_stall;
%         % Stores values into de restrictions so that can be plotted
%         Restrictions.V_stall = V_stall;
%         Restrictions.V_TO = V_TO;
%         Restrictions.V_min_ope = V_min_ope;
%         % Calculations
%         [TRIM_RESULTS_calc_var_V_m_v,Trim_ITER_calc_var_V_m_v,Stab_Der_calc_var_V_m_v,Stab_Der_parts_calc_var_V_m_v,Propulsion_calc_var_V_m_v] = ...
%             Calculo_Stability_Derivatives_April2022_v1(conditions,Aero,Aero_TH,Geo_tier,Weight_tier,Prop_data,conv_UNITS,...
%             Body_Geo,Design_criteria,Posicion_Palanca,XCG_FF,DATA_Ae,Performance,OUTPUT_read_XLSX,AC_CONFIGURATION,only_trim,Performance_preliminar);
%         % Storing
%         Restrictions_var_V_m{i,j} = Restrictions;
%         TRIM_RESULTS_var_V_m{i,j} = TRIM_RESULTS_calc_var_V_m_v;
%         Trim_ITER_var_V_m{i,j} = Trim_ITER_calc_var_V_m_v;
%         Stab_Der_var_V_m{i,j} = Stab_Der_calc_var_V_m_v;
%         Stab_Der_parts_V_m{i,j} = Stab_Der_parts_calc_var_V_m_v;
%         Propulsion_var_V_m{i,j} =  Propulsion_calc_var_V_m_v;
% 
%         StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
%         Variable_Study = 1;
% 
%         %ESTABILIDAD DINAMICA LONGITUDINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn == 1
%             show_messages_screen = 0;
%             Stab_Dyn_Long_var_V_m{i,j} = longitudinal_analysis(Performance,Stab_Der_calc_var_V_m_v,conv_UNITS,StabilityModel,...
%                 Weight_tier,TRIM_RESULTS_calc_var_V_m_v,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
%         else
%             Dummy = 0;
%             Stab_Dyn_Long_var_V_m = Dummy;
%         end
%         %ESTABILIDAD DINAMICA LATERAL-DIRECCIONAL%%%%%%%%%%%%%%%%%%%%%%
%         if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn == 1
%             show_messages_screen = 0;
%             Stab_Dyn_LatDir_var_V_m{i,j} = lateral_directional_analysis(Performance,Stab_Der_calc_var_V_m_v,conv_UNITS,StabilityModel,...
%                 Weight_tier,TRIM_RESULTS_calc_var_V_m_v,Geo_tier,Variable_Study,conditions,OUTPUT_read_XLSX,show_messages_screen);
%         else
%             Dummy = 0;
%             Stab_Dyn_LatDir_var_V_m = Dummy;
%         end
%     end
% end

% Preallocate variables
x_XCG_VAR = zeros(N_V_VAR, N_m_VAR); % Preallocate numeric arrays
z_XCG_VAR = zeros(N_V_VAR, N_m_VAR);

Restrictions_var_V_m = cell(N_V_VAR, N_m_VAR); % Preallocate cell arrays
TRIM_RESULTS_var_V_m = cell(N_V_VAR, N_m_VAR);
Trim_ITER_var_V_m = cell(N_V_VAR, N_m_VAR);
Stab_Der_var_V_m = cell(N_V_VAR, N_m_VAR);
Stab_Der_parts_V_m = cell(N_V_VAR, N_m_VAR);
Propulsion_var_V_m = cell(N_V_VAR, N_m_VAR);

Stab_Dyn_Long_var_V_m = cell(N_V_VAR, N_m_VAR);
Stab_Dyn_LatDir_var_V_m = cell(N_V_VAR, N_m_VAR);

% Precompute constants if independent of inner loop
g = conv_UNITS.g;
rho = Performance_preliminar.rho;
S_w1 = Geo_tier.S_w1;

StabilityModel = OUTPUT_read_XLSX.Stability_flags.StabilityModel;
Variable_Study = 1;

% Outer loop: Iterate over velocity variables
for i = 1:N_V_VAR
    conditions.V = V_VAR(i);
    
    % Inner loop: Iterate over mass variables
    for j = 1:N_m_VAR
        % Set up conditions
        conditions.m_TOW = m_VAR(j);
        gamma = 0;
        conditions.gamma = gamma;
        
        % Calculate XCG variations
        if OUTPUT_read_XLSX.Weights_flags.get_shift_XCG_variation == 1
            [conditions, OUTPUT_read_XLSX] = get_Variable_x_XCG(m_VAR(j), case_AC, conditions, OUTPUT_read_XLSX, Modifications);
            x_XCG_VAR(i, j) = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
            z_XCG_VAR(i, j) = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;
        else
            x_XCG_VAR(i, j) = 1; % Default value
            z_XCG_VAR(i, j) = 1; % Default value
        end
        
        % Calculate stall speed and operating speeds
        V_stall = sqrt(2 * m_VAR(j) * g / (S_w1 * rho * CL_max_ac));
        V_TO = 1.2 * V_stall;
        V_min_ope = 1.3 * V_stall;
        
        % Store restrictions
        Restrictions.V_stall = V_stall;
        Restrictions.V_TO = V_TO;
        Restrictions.V_min_ope = V_min_ope;
        Restrictions_var_V_m{i, j} = Restrictions;
        
        % Calculate stability derivatives and store results
        [TRIM_RESULTS_calc_var, Trim_ITER_calc_var, Stab_Der_calc_var, Stab_Der_parts_calc_var, Propulsion_calc_var] = ...
            Calculo_Stability_Derivatives_April2022_v1(conditions, Aero, Aero_TH, Geo_tier, Weight_tier, Prop_data, conv_UNITS, ...
            Body_Geo, Design_criteria, Posicion_Palanca, XCG_FF, DATA_Ae, Performance, OUTPUT_read_XLSX, AC_CONFIGURATION, only_trim, Performance_preliminar);
        
        TRIM_RESULTS_var_V_m{i, j} = TRIM_RESULTS_calc_var;
        Trim_ITER_var_V_m{i, j} = Trim_ITER_calc_var;
        Stab_Der_var_V_m{i, j} = Stab_Der_calc_var;
        Stab_Der_parts_V_m{i, j} = Stab_Der_parts_calc_var;
        Propulsion_var_V_m{i, j} = Propulsion_calc_var;
        
        % Longitudinal dynamic stability
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn
            Stab_Dyn_Long_var_V_m{i, j} = longitudinal_analysis(Performance, Stab_Der_calc_var, conv_UNITS, ...
                StabilityModel, Weight_tier, TRIM_RESULTS_calc_var, Geo_tier, 1, conditions, OUTPUT_read_XLSX, 0);
        else
            Stab_Dyn_Long_var_V_m{i, j} = 0; % Default value
        end
        
        % Lateral-directional dynamic stability
        if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn
            Stab_Dyn_LatDir_var_V_m{i, j} = lateral_directional_analysis(Performance, Stab_Der_calc_var, conv_UNITS, ...
                StabilityModel, Weight_tier, TRIM_RESULTS_calc_var, Geo_tier, 1, conditions, OUTPUT_read_XLSX, 0);
        else
            Stab_Dyn_LatDir_var_V_m{i, j} = 0; % Default value
        end
    end
end

% pause
% tic
% % Preallocate the variables
% x_XCG_VAR = zeros(N_V_VAR, N_m_VAR);
% z_XCG_VAR = zeros(N_V_VAR, N_m_VAR);
% 
% Restrictions_var_V_m = cell(N_V_VAR, N_m_VAR);
% TRIM_RESULTS_var_V_m = cell(N_V_VAR, N_m_VAR);
% Trim_ITER_var_V_m = cell(N_V_VAR, N_m_VAR);
% Stab_Der_var_V_m = cell(N_V_VAR, N_m_VAR);
% Stab_Der_parts_V_m = cell(N_V_VAR, N_m_VAR);
% Propulsion_var_V_m = cell(N_V_VAR, N_m_VAR);
% 
% Stab_Dyn_Long_var_V_m = cell(N_V_VAR, N_m_VAR);
% Stab_Dyn_LatDir_var_V_m = cell(N_V_VAR, N_m_VAR);
% 
% % Start parfor loop
% parfor i = 1:N_V_VAR
%     % Temporary storage for outer loop iteration
%     temp_Restrictions_var = cell(1, N_m_VAR);
%     temp_TRIM_RESULTS_var = cell(1, N_m_VAR);
%     temp_Trim_ITER_var = cell(1, N_m_VAR);
%     temp_Stab_Der_var = cell(1, N_m_VAR);
%     temp_Stab_Der_parts_var = cell(1, N_m_VAR);
%     temp_Propulsion_var = cell(1, N_m_VAR);
%     temp_Stab_Dyn_Long_var = cell(1, N_m_VAR);
%     temp_Stab_Dyn_LatDir_var = cell(1, N_m_VAR);
% 
%     for j = 1:N_m_VAR
%         % Set up conditions
%         conditions.V = V_VAR(i);
%         conditions.m_TOW = m_VAR(j);
% 
%         % Calculate XCG variations
%         if OUTPUT_read_XLSX.Weights_flags.get_shift_XCG_variation == 1
%             [conditions, OUTPUT_read_XLSX] = get_Variable_x_XCG(m_VAR(j), case_AC, conditions, OUTPUT_read_XLSX, Modifications);
%             x_XCG_VAR(i, j) = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG;
%             z_XCG_VAR(i, j) = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG;
%         else
%             x_XCG_VAR(i, j) = 1; % Default value
%             z_XCG_VAR(i, j) = 1; % Default value
%         end
% 
%         % Calculate stall speed and operating speeds
%         V_stall = sqrt(2 * m_VAR(j) * conv_UNITS.g / (Geo_tier.S_w1 * Performance_preliminar.rho * CL_max_ac));
%         V_TO = 1.2 * V_stall;
%         V_min_ope = 1.3 * V_stall;
% 
%         % Store restrictions directly in temporary cell
%         temp_Restrictions_var{j} = struct( ...
%             'V_stall', V_stall, ...
%             'V_TO', V_TO, ...
%             'V_min_ope', V_min_ope ...
%         );
% 
%         % Calculate stability derivatives and store results
%         [TRIM_RESULTS_calc_var, Trim_ITER_calc_var, Stab_Der_calc_var, Stab_Der_parts_calc_var, Propulsion_calc_var] = ...
%             Calculo_Stability_Derivatives_April2022_v1(conditions, Aero, Aero_TH, Geo_tier, Weight_tier, Prop_data, conv_UNITS, ...
%             Body_Geo, Design_criteria, Posicion_Palanca, XCG_FF, DATA_Ae, Performance, OUTPUT_read_XLSX, AC_CONFIGURATION, only_trim, Performance_preliminar);
% 
%         temp_TRIM_RESULTS_var{j} = TRIM_RESULTS_calc_var;
%         temp_Trim_ITER_var{j} = Trim_ITER_calc_var;
%         temp_Stab_Der_var{j} = Stab_Der_calc_var;
%         temp_Stab_Der_parts_var{j} = Stab_Der_parts_calc_var;
%         temp_Propulsion_var{j} = Propulsion_calc_var;
% 
%         % Longitudinal dynamic stability
%         if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_Long_dyn
%             temp_Stab_Dyn_Long_var{j} = longitudinal_analysis(Performance, Stab_Der_calc_var, conv_UNITS, ...
%                 StabilityModel, Weight_tier, TRIM_RESULTS_calc_var, Geo_tier, 1, conditions, OUTPUT_read_XLSX, 0);
%         else
%             temp_Stab_Dyn_Long_var{j} = 0; % Default value
%         end
% 
%         % Lateral-directional dynamic stability
%         if OUTPUT_read_XLSX.STUDY_flags.STABILITY_STUDY_LatDir_dyn
%             temp_Stab_Dyn_LatDir_var{j} = lateral_directional_analysis(Performance, Stab_Der_calc_var, conv_UNITS, ...
%                 StabilityModel, Weight_tier, TRIM_RESULTS_calc_var, Geo_tier, 1, conditions, OUTPUT_read_XLSX, 0);
%         else
%             temp_Stab_Dyn_LatDir_var{j} = 0; % Default value
%         end
%     end
% 
%     % Assign temporary variables back to main storage
%     Restrictions_var_V_m(i, :) = temp_Restrictions_var;
%     TRIM_RESULTS_var_V_m(i, :) = temp_TRIM_RESULTS_var;
%     Trim_ITER_var_V_m(i, :) = temp_Trim_ITER_var;
%     Stab_Der_var_V_m(i, :) = temp_Stab_Der_var;
%     Stab_Der_parts_V_m(i, :) = temp_Stab_Der_parts_var;
%     Propulsion_var_V_m(i, :) = temp_Propulsion_var;
%     Stab_Dyn_Long_var_V_m(i, :) = temp_Stab_Dyn_Long_var;
%     Stab_Dyn_LatDir_var_V_m(i, :) = temp_Stab_Dyn_LatDir_var;
% end
% toc
% pause

%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1
    Storing_STABILITY_DATA_2 = Saving_data_Stability_varV(TRIM_RESULTS_var_V_m,Trim_ITER_var_V_m,Stab_Der_var_V_m,Stab_Der_parts_V_m,...
        Stab_Dyn_Long_var_V_m,Stab_Dyn_LatDir_var_V_m,case_AC,Propulsion_var_V_m,Restrictions_var_V_m,Plot_Options,OUTPUT_read_XLSX,x_XCG_VAR,filenameS);
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_2 = dummy;
end