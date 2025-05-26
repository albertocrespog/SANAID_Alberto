function Storing_STABILITY_DATA_4C = Conduct_Stability_Study4C(conditions,conv_UNITS,Storing_STABILITY_DATA_1,OUTPUT_read_XLSX,...
    Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,case_AC,filenameS)

Geo_tier = Storing_GEO_DATA.Geo_tier;
Body_Geo = Storing_GEO_DATA.Body_Geo;
Performance = Storing_AERO_DATA.Performance;
Aero = Storing_AERO_DATA.Aero;
Weight_tier = Storing_WEIGHT_DATA.Weight_tier;

%%%%%%%%%%%%%%%%%%%% TRIM LATERAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions.study_var_xcg = 0;
% Cosntant beta
beta = OUTPUT_read_XLSX.Stability_flags.beta;
% variable beta
beta_i = OUTPUT_read_XLSX.Stability_flags.beta_i;
beta_f = OUTPUT_read_XLSX.Stability_flags.beta_f;
N_Delta_beta = OUTPUT_read_XLSX.Stability_flags.N_Delta_beta;
beta_vec = linspace(beta_i,beta_f,N_Delta_beta);
% Storing study conditions
conditions_TRIM_lat.beta = beta;
conditions_TRIM_lat.beta_vec = beta_vec;

%% Weight range
% m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_systems + Weight_tier.m_batteries/2 + Weight_tier.m_fuel/2);
m_low = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy);
m_mid = Weight_tier.m_TOW - (Weight_tier.m_payload + Weight_tier.m_energy)/2;
m_high = Weight_tier.m_TOW;

% Determines stall conditions
m_TOW = conditions.m_TOW;
%% Updates the Velocity Vector according to the minimum speeds
V_stall_min_low = sqrt(2*m_low*conv_UNITS.g/(Geo_tier.S_w1*Performance.rho*Aero.CL_max_w1_CR));
V_stall_min_high = sqrt(2*m_high*conv_UNITS.g/(Geo_tier.S_w1*Performance.rho*Aero.CL_max_w1_CR));
V_TO = 1.2*V_stall_min_low;
V_min_ope = 1.3*V_stall_min_low;

conditions_TRIM_lat.V_TO = V_TO;
conditions_TRIM_lat.m_TOW = m_TOW;

%% Speed range
V_low = OUTPUT_read_XLSX.Stability_flags.V_low;
V_high = OUTPUT_read_XLSX.Stability_flags.V_high;
N_V_VAR = OUTPUT_read_XLSX.Stability_flags.N_V_VAR;
V_min_ope_low = 1.3*V_stall_min_low;
V_min_ope_high = 1.3*V_stall_min_high;

% Correct the minimum speed  from excel
V_low = V_stall_min_low;

% Correct the maximum speed from excel
if V_high < V_min_ope_high
    V_high = V_min_ope_high;
end

% variable speed conditions
conditions_TRIM_lat.V_VAR = linspace(V_low,V_high,N_V_VAR);

% rho = 1.225; % At sea level
% rho = 0.96296; % At 8000 ft
% rho = 0.74628; % At 16000 ft
conditions_TRIM_lat.rho = Performance.rho;

% Trim lateral acelertions for deflections
[Trim_ITER_LAT3] = Calculo_Trim_ITER_LAT_accelerations(conv_UNITS,conditions_TRIM_lat,...
    Storing_GEO_DATA,Storing_AERO_DATA,Storing_WEIGHT_DATA,Storing_STABILITY_DATA_1,conditions,case_AC,OUTPUT_read_XLSX);


%% Saves the data to a mat file
if OUTPUT_read_XLSX.MATLAB_flags.save_MAT == 1   
    Storing_STABILITY_DATA_4C = Saving_Trim_LatC(Trim_ITER_LAT3,conditions_TRIM_lat,OUTPUT_read_XLSX,filenameS);   
else
    % Stores dummy variable for continuity
    dummy = 1;
    Storing_STABILITY_DATA_4C = dummy;
end